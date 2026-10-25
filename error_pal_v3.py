#To evaluate the preinjector errors.
#This is pal version, that can run the multiple parmela at once. Limited by the # of CPU
#author: Erdong Wang
#Usage: python error_ana_pal.py <input_filename.inp> <error_config.yaml>
#Version 1.5 PDF summary added


import os
import sys
import shutil
import yaml
import glob
import uuid
import subprocess
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import matplotlib.pyplot as plt


# 1. Find element indices and main frequency
def find_ele_ind(lines):
    elements = {
        "quad": [],
        "solenoid": [],
        "cell1300": [],
        "cell197": [],
        "trwave": [],
        "bend": [],
        "steerer": [],
        "laser": [],
    }
    mainfreq = None
    for i, raw in enumerate(lines):
        line = raw.strip().lower()
        if line.startswith("run"):
            parts = raw.split()
            if len(parts) > 3:
                mainfreq = float(parts[3])

        if line.startswith("quad"):
            elements["quad"].append(i)
        elif line.startswith("solenoid"):
            elements["solenoid"].append(i)
        elif line.startswith("bend"):
            elements["bend"].append(i)
        elif line.startswith("steerer"):
            elements["steerer"].append(i)
        elif line.startswith("input"):
            elements["laser"].append(i)
        elif line.startswith("cell"):
            # A CELL followed immediately by TRWAVE is the header of a traveling-wave structure.
            if i + 1 < len(lines) and lines[i + 1].strip().lower().startswith("trwave"):
                elements["trwave"].append(i)
            else:
                parts = raw.split()
                if len(parts) > 9:
                    try:
                        rf_freq = int(float(parts[9]))
                    except ValueError:
                        continue
                    if rf_freq == 1300:
                        elements["cell1300"].append(i)
                    elif rf_freq == 197:
                        elements["cell197"].append(i)
    return elements, mainfreq

# 2. Truncated normal helper
def truncated_normal(mean, sigma, bound, size):
    from scipy.stats import truncnorm
    import numpy as np
    if sigma == 0:
        return np.full(size, mean)
    a, b = (-bound - mean) / sigma, (bound - mean) / sigma
    return truncnorm.rvs(a, b, loc=mean, scale=sigma, size=size)


def _param(params, name, fallback_name=None):
    """Read a YAML parameter, with optional backward-compatible fallback."""
    if name in params:
        return params[name]
    if fallback_name is not None and fallback_name in params:
        return params[fallback_name]
    raise KeyError(f"Missing required YAML parameter: {name}")


def _rf_distribution(elements, params, element_key, yaml_prefix, fallback_prefix=None):
    phase = truncated_normal(
        _param(params, f"{yaml_prefix}_rf_phase_mean", f"{fallback_prefix}_rf_phase_mean" if fallback_prefix else None),
        _param(params, f"{yaml_prefix}_rf_phase_sig", f"{fallback_prefix}_rf_phase_sig" if fallback_prefix else None),
        _param(params, f"{yaml_prefix}_rf_phase_bound", f"{fallback_prefix}_rf_phase_bound" if fallback_prefix else None),
        len(elements[element_key]),
    )
    amp = truncated_normal(
        _param(params, f"{yaml_prefix}_rf_amp_mean", f"{fallback_prefix}_rf_amp_mean" if fallback_prefix else None),
        _param(params, f"{yaml_prefix}_rf_amp_sig", f"{fallback_prefix}_rf_amp_sig" if fallback_prefix else None),
        _param(params, f"{yaml_prefix}_rf_amp_bound", f"{fallback_prefix}_rf_amp_bound" if fallback_prefix else None),
        len(elements[element_key]),
    )
    return list(zip(phase, amp))

# 3. Generate distributions
def randseed(elements, params):
    dist = {}
    # quadrupole, solenoid and steerer use power supply amp
    dist["quad"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["quad"]))
    dist["solenoid"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["solenoid"]))
    dist["steerer"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["steerer"]))

    # standing-wave cells, split by RF frequency in part 10 of the CELL line.
    # Fallback to old cell_* YAML names so older YAML files still work.
    dist["cell1300"] = _rf_distribution(elements, params, "cell1300", "cell1300", fallback_prefix="cell")
    dist["cell197"] = _rf_distribution(elements, params, "cell197", "cell197", fallback_prefix="cell")

    # trwave
    dist["trwave"] = _rf_distribution(elements, params, "trwave", "trwave")

    # laser INPUT jitter
    # laser_longi_jitter is read in ps and converted to degrees by ps*2856*360*1e-6.
    laser_longi_ps = truncated_normal(
        params["laser_longi_jitter_mean"],
        params["laser_longi_jitter_sig"],
        params["laser_longi_jitter_bound"],
        len(elements["laser"]),
    )
    laser_longi_deg = laser_longi_ps * 2856.0 * 360.0 * 1.0e-6
    laser_x = truncated_normal(
        params["laser_xjitter_mean"],
        params["laser_xjitter_sig"],
        params["laser_xjitter_bound"],
        len(elements["laser"]),
    )
    laser_y = truncated_normal(
        params["laser_yjitter_mean"],
        params["laser_yjitter_sig"],
        params["laser_yjitter_bound"],
        len(elements["laser"]),
    )
    dist["laser"] = list(zip(laser_longi_deg, laser_x, laser_y))

    # bend
    dist["bend"] = truncated_normal(params["bend_amp_mean"], params["bend_amp_sig"], params["bend_amp_bound"], len(elements["bend"]))

    return dist


def _split_data_and_comment(raw):
    """Split a PARMELA command into editable fields plus the inline comment, if present."""
    newline = "\n" if raw.endswith("\n") else ""
    body = raw[:-1] if newline else raw
    comment_start = len(body)
    for marker in ("!", ";"):
        pos = body.find(marker)
        if pos != -1:
            comment_start = min(comment_start, pos)
    data = body[:comment_start].strip()
    comment = body[comment_start:].strip() if comment_start < len(body) else ""
    return data.split(), comment, newline


def _join_data_and_comment(parts, comment, newline):
    data = " ".join(parts)
    if comment:
        data = f"{data} {comment}"
    return data + newline


def _add_to_field(parts, field_index, delta):
    """Add delta to a 0-based split-field index, padding optional missing fields with 0."""
    while len(parts) <= field_index:
        parts.append("0")
    parts[field_index] = str(float(parts[field_index]) + delta)


def _apply_cell_perturbation(lines, indices, distribution):
    for idx, (p, a) in zip(indices, distribution):
        parts = lines[idx].split()
        if len(parts) >= 6:
            if p != 0:
                parts[4] = str(float(parts[4]) + p)
            if a != 0:
                parts[5] = str(float(parts[5]) * (1 + a))
            lines[idx] = ' '.join(parts) + '\n'

# 4. Apply element perturbations
def apply_perturbations(input_path, params, folder):
    # read lines
    with open(input_path, 'r') as f:
        lines = f.readlines()
    elements, mainfreq = find_ele_ind(lines)
    dist = randseed(elements, params)

    # modify trwave
    for idx, (p, a) in zip(elements["trwave"], dist["trwave"]):
        for j in range(idx, idx + 86):
            if j >= len(lines): break
            parts = lines[j].split()
            if len(parts) >= 6:
                if p != 0:
                    parts[4] = str(float(parts[4]) + p)
                if a != 0:
                    parts[5] = str(float(parts[5]) * (1 + a))
                lines[j] = ' '.join(parts) + '\n'

    # modify standing-wave cells
    _apply_cell_perturbation(lines, elements["cell1300"], dist["cell1300"])
    _apply_cell_perturbation(lines, elements["cell197"], dist["cell197"])

    # laser INPUT jitter
    # parts[6]: longitudinal laser jitter, degrees
    # parts[9]: x laser jitter, cm
    # parts[10]: y laser jitter, cm
    for idx, (longi_deg, x_jitter, y_jitter) in zip(elements["laser"], dist["laser"]):
        if longi_deg == 0 and x_jitter == 0 and y_jitter == 0:
            continue
        parts, comment, newline = _split_data_and_comment(lines[idx])
        if longi_deg != 0:
            _add_to_field(parts, 7, longi_deg)
        if x_jitter != 0:
            _add_to_field(parts, 10, x_jitter)
        if y_jitter != 0:
            _add_to_field(parts, 11, y_jitter)
        lines[idx] = _join_data_and_comment(parts, comment, newline)

    # solenoid
    for idx, a in zip(elements["solenoid"], dist["solenoid"]):
        parts = lines[idx].split()
        if len(parts) >= 5 and a != 0:
            parts[4] = str(float(parts[4]) * (1 + a))
            lines[idx] = ' '.join(parts) + '\n'
    # quad
    for idx, a in zip(elements["quad"], dist["quad"]):
        parts = lines[idx].split()
        if len(parts) >= 5 and a != 0:
            parts[4] = str(float(parts[4]) * (1 + a))
            lines[idx] = ' '.join(parts) + '\n'
 
    # steerer
    for idx, a in zip(elements["steerer"], dist["steerer"]):
        parts = lines[idx].split()
        if len(parts) >= 6 and a != 0:
            parts[4] = str(float(parts[4]) * (1 + a))
            parts[5] = str(float(parts[5]) * (1 + a))
            lines[idx] = ' '.join(parts) + '\n'
 
     # bend
    for idx, a in zip(elements["bend"], dist["bend"]):
        parts = lines[idx].split()
        if len(parts) >= 5 and a != 0:
            parts[4] = str(float(parts[4]) * (1 + a))
            lines[idx] = ' '.join(parts) + '\n'

    # write modified file
    pert_file = os.path.join(folder, os.path.basename(input_path).replace('.inp', '_erranaly.inp'))
    with open(pert_file, 'w') as f:
        f.writelines(lines)
    return pert_file

# 5. Run a single case
def run_case(args):
    folder, pert_file = args
    cwd = os.getcwd()
    os.chdir(folder)
    try:
        subprocess.run(["parmela", os.path.basename(pert_file)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # Find all files ending with .T2 or .T3
        files_to_delete = glob.glob('*.T2') + glob.glob('*.T3')
        # Loop through the list and delete each file
        for f in files_to_delete:
            try:
                os.remove(f)
            except OSError as e:
                print(f"Error deleting file {f}: {e}")
        # -------------------------------------------
    finally:
        os.chdir(cwd)
    # return path to output TBL
    tbl_files = glob.glob(os.path.join(folder, 'TIMESTEPEMITTANCE.TBL'))
    if not tbl_files:
        tbl_files = glob.glob(os.path.join(folder, '*.tbl'))
        if tbl_files:
            print(f"Warning: TIMESTEPEMITTANCE.TBL not found. Using {os.path.basename(tbl_files[0])} instead.")

    if tbl_files:
        return tbl_files[0]
    return None


# 6. Aggregate results
def aggregate_results(tbl_paths, params, run_id):
    # prepare header in main analysis file
    main_file = "error_analysis_dat.txt"
    header = "T(deg) Z(cm) Xun(mm-mrad) Yun(mm-mrad) Zun(mm-mrad) Xn(mm-mrad) Yn(mm-mrad) Zn(mm-mrad) Xrms(mm) Yrms(mm) Zrmz(mm) <kE>(MeV) Del-Erms <X>(mm) <Xpn>(mrad) <Y>(mm) <Ypn>(mrad) <Z>(cm) <Zpn>(rad) EZref(MV/m)"
    with open(main_file, 'w') as out:
        out.write(header + '\n')
    # append each -3 line
    for tbl in tbl_paths:
        if tbl is None: continue
        with open(tbl, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 3:
                line = lines[-3].strip()
                with open(main_file, 'a') as out:
                    out.write(line + '\n')
    # rename
    out_file = f"error_analysis_dat_{run_id}.txt"
    shutil.move(main_file, out_file)
    return out_file


BEAM_PARAMETER_EXPLANATIONS = {
    "T(deg)": "Time",
    "Z(cm)": "Beam position",
    "Xun(mm-mrad)": "RMS un-normalized x emittance",
    "Yun(mm-mrad)": "RMS un-normalized y emittance",
    "Zun(mm-mrad)": "RMS un-normalized z emittance",
    "Xn(mm-mrad)": "RMS normalized x emittance",
    "Yn(mm-mrad)": "RMS normalized y emittance",
    "Zn(mm-mrad)": "RMS normalized z emittance",
    "Xrms(mm)": "RMS x radius",
    "Yrms(mm)": "RMS y radius",
    "Zrmz(mm)": "RMS z radius",
    "<kE>(MeV)": "Average kinetic energy",
    "Del-Erms": "dE/E",
    "<X>(mm)": "x orbit",
    "<Xpn>(mrad)": "Average normalized x angular error",
    "<Y>(mm)": "y orbit",
    "<Ypn>(mrad)": "Average normalized y angular error",
    "<Z>(cm)": "z orbit",
    "<Zpn>(rad)": "Average normalized z angular error",
    "EZref(MV/m)": "RF field",
}


def _fmt_float(value):
    """Compact float formatting for PDF tables."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value == 0:
        return "0"
    abs_value = abs(value)
    if abs_value >= 1.0e5 or abs_value < 1.0e-3:
        return f"{value:.6e}"
    return f"{value:.6g}"


def _split_yaml_error_key(key):
    for suffix in ("_mean", "_sig", "_bound"):
        if key.endswith(suffix):
            return key[:-len(suffix)], suffix[1:]
    return None, None


def _infer_error_unit(param_name):
    name = param_name.lower()
    if "laser_longi_jitter" in name:
        return "ps"
    if "laser_xjitter" in name or "laser_yjitter" in name:
        return "cm"
    if "phase" in name:
        return "deg"
    if "amp" in name:
        return "fraction"
    return ""


def _yaml_error_setup_rows(params):
    """Return rows [parameter, unit, mean, sigma, bound] from YAML settings."""
    grouped = {}
    order = []
    for key, value in params.items():
        base, field = _split_yaml_error_key(str(key))
        if base is None:
            continue
        if base not in grouped:
            grouped[base] = {"mean": "", "sig": "", "bound": ""}
            order.append(base)
        grouped[base][field] = value

    rows = []
    for base in order:
        values = grouped[base]
        rows.append([
            base,
            _infer_error_unit(base),
            _fmt_float(values.get("mean", "")),
            _fmt_float(values.get("sig", "")),
            _fmt_float(values.get("bound", "")),
        ])
    return rows


def _read_error_analysis_dat(error_dat_path):
    """Read error_analysis_dat_[run_id].txt and force numeric columns."""
    df = pd.read_csv(error_dat_path, sep=r"\s+", engine="python")
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _beam_statistics_rows(df):
    """Return rows [parameter, description, average, std, avg-3std, avg+3std]."""
    stats = []
    for col in df.columns:
        series = df[col].dropna()
        if series.empty:
            continue
        avg = series.mean()
        # For Monte Carlo runs, use sample standard deviation. If only one run exists, use 0.
        std = series.std(ddof=1) if len(series) > 1 else 0.0
        stats.append([
            col,
            BEAM_PARAMETER_EXPLANATIONS.get(col, ""),
            _fmt_float(avg),
            _fmt_float(std),
            _fmt_float(avg - 3.0 * std),
            _fmt_float(avg + 3.0 * std),
        ])
    return stats


def _add_table(story, table_data, widths, repeat_rows=1, font_size=7):
    from reportlab.lib import colors
    from reportlab.platypus import LongTable

    table = LongTable(table_data, colWidths=widths, repeatRows=repeat_rows)
    table.setStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8EEF7")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#111111")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("ALIGN", (2, 1), (-1, -1), "RIGHT"),
        ("ALIGN", (1, 1), (1, -1), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#B8B8B8")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7F7F7")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ])
    story.append(table)


def generate_error_summary_pdf(params, error_dat_path, run_id, output_pdf=None):
    """Generate error_Summary_[run_id].pdf from YAML setup and aggregate beam results."""
    try:
        from reportlab.lib.pagesizes import letter, landscape
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
    except ImportError as exc:
        print("Warning: reportlab is required to generate the PDF summary. Install it with: pip install reportlab")
        print(f"PDF summary skipped: {exc}")
        return None

    if output_pdf is None:
        output_pdf = f"error_Summary_{run_id}.pdf"

    df = _read_error_analysis_dat(error_dat_path)
    setup_rows = _yaml_error_setup_rows(params)
    stats_rows = _beam_statistics_rows(df)

    doc = SimpleDocTemplate(
        output_pdf,
        pagesize=landscape(letter),
        rightMargin=24,
        leftMargin=24,
        topMargin=24,
        bottomMargin=24,
    )
    styles = getSampleStyleSheet()
    story = []

    story.append(Paragraph(f"Error Summary - Run ID: {run_id}", styles["Title"]))
    story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles["Normal"]))
    story.append(Paragraph(f"Runs in YAML: {params.get('runs', '')}; rows in error_analysis_dat: {len(df)}", styles["Normal"]))
    story.append(Spacer(1, 10))

    story.append(Paragraph("Initial error setup from YAML", styles["Heading2"]))
    story.append(Paragraph("The table lists the truncated-normal mean+-bound values with sigma values read from the YAML file. ", styles["Normal"]))
    story.append(Spacer(1, 6))
    setup_table = [["Error parameter", "Unit", "Mean", "Sigma", "Bound"]] + setup_rows
    _add_table(story, setup_table, [230, 70, 100, 100, 100], repeat_rows=1, font_size=8)

    story.append(PageBreak())
    story.append(Paragraph("Beam parameter statistics from error_analysis_dat", styles["Heading2"]))
    story.append(Paragraph("Average and sample standard deviation are computed across all runs. The last two columns give average - 3*std and average + 3*std.", styles["Normal"]))
    story.append(Spacer(1, 6))
    stats_table = [["Parameter", "Explanation", "Average", "Std", "Avg - 3*Std", "Avg + 3*Std"]] + stats_rows
    _add_table(story, stats_table, [88, 205, 85, 85, 95, 95], repeat_rows=1, font_size=7)

    doc.build(story)
    print(f"Error summary PDF saved to {output_pdf}")
    return output_pdf


# 7. Orbit Figure Plotting
def orbit_figure(base_name, runs, run_id):
    """
    Generates orbit error plots from multiple simulation runs.
    """
    print("Generating orbit error plots...")
    output_txt_file = f"orbit_error_{run_id}.txt"
    
    all_data = pd.DataFrame()

    for i in range(1, runs + 1):
        folder = f"{base_name}_{i}"
        tbl_files = glob.glob(os.path.join(folder, 'TIMESTEPEMITTANCE.TBL'))
        if not tbl_files:
            print(f"Warning: No TIMESTEPEMITTANCE.TBL file found in {folder}")
            continue
        
        tbl_path = tbl_files[0]

        with open(tbl_path, 'r') as f:
            lines = f.readlines()
        
        line_ini = -1
        for idx, line in enumerate(lines):
            if "DATA" in line:
                line_ini = idx
                break

        if line_ini == -1:
            print(f"Warning: 'DATA' keyword not found in {tbl_path}")
            continue
            
        df = pd.read_csv(tbl_path, skiprows=line_ini + 2, sep=r'\s+', header=None, engine='python')
        
        # Read header from the line after "DATA" and remove leading semicolon
        header = lines[line_ini + 1].strip().lstrip(';').split()

        if len(df.columns) != len(header):
            print(f"Error: Column count mismatch in file {tbl_path}.")
            print(f"Data columns: {len(df.columns)}, Header columns: {len(header)}")
            continue

        df.columns = header

        if i == 1:
            # Only copy Z(cm), <X>(mm), <Y>(mm) for the first file
            if all(col in df.columns for col in ['Z(cm)', '<X>(mm)', '<Y>(mm)']):
                all_data = df[['Z(cm)', '<X>(mm)', '<Y>(mm)']].copy()
                # Rename the columns for the first file
                all_data.columns = [f"Z(cm)_1", f"<X>(mm)_1", f"<Y>(mm)_1"]
            else:
                print(f"Warning: Key columns 'Z(cm)', '<X>(mm)', or '<Y>(mm)' not found in {tbl_path}. Skipping.")
                continue
        else:
            # For subsequent files, only take <X>(mm) and <Y>(mm)
            if '<X>(mm)' in df.columns and '<Y>(mm)' in df.columns:
                df_subset = df[['<X>(mm)', '<Y>(mm)']].copy()
                df_subset.columns = [f"<X>(mm)_{i}", f"<Y>(mm)_{i}"]
                all_data = pd.concat([all_data, df_subset], axis=1)

    # Save the combined data to a file
    all_data.to_csv(output_txt_file, sep='\t', index=False)
    print(f"Combined orbit data saved to {output_txt_file}")

    # Plotting
    z_col = f'Z(cm)_1'
    if z_col not in all_data.columns:
        print("Error: 'Z(cm)_1' column not found for plotting.")
        return

    # Plot X orbits
    plt.figure(figsize=(10, 6))
    for i in range(1, runs + 1):
        x_col = f'<X>(mm)_{i}'
        if x_col in all_data.columns:
            plt.plot(all_data[z_col], all_data[x_col])
    plt.xlabel('Z (cm)')
    plt.ylabel('<X> (mm)')
    plt.title(f'Orbit Error Analysis (X) - Run ID: {run_id}')
    plt.grid(True)
    plt.savefig(f'orbit_X_{run_id}.png', dpi=300)
    plt.show()

    # Plot Y orbits
    plt.figure(figsize=(10, 6))
    for i in range(1, runs + 1):
        y_col = f'<Y>(mm)_{i}'
        if y_col in all_data.columns:
            plt.plot(all_data[z_col], all_data[y_col])
    plt.xlabel('Z (cm)')
    plt.ylabel('<Y> (mm)')
    plt.title(f'Orbit Error Analysis (Y) - Run ID: {run_id}')
    plt.grid(True)
    plt.savefig(f'orbit_Y_{run_id}.png', dpi=300)
    plt.show()


# 8. Main entry
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python error_ana_pal.py <input_filename.inp> <error_config.yaml>")
        sys.exit(1)

    input_filename = sys.argv[1]
    yaml_filename = sys.argv[2]
    with open(yaml_filename) as f:
        params = yaml.safe_load(f)
    runs = params.get('runs', 1)
    base = os.path.splitext(input_filename)[0]
    
    
    # ################### HIGHLIGHTED CHANGE ###################
    # Clean up old directories before starting a new run
    print("Cleaning up old simulation directories...")
    old_folders = glob.glob(f"{base}_*")
    for folder in old_folders:
        if os.path.isdir(folder):
            print(f"Removing: {folder}")
            shutil.rmtree(folder)
    # ################# END HIGHLIGHTED CHANGE #################

    # unique ID
    run_id = uuid.uuid4().hex[:8]
    shutil.copyfile(yaml_filename, f"error_{run_id}.yaml")

    # create folders and prepare cases
    tbl_paths = []
    args = []
    for i in range(1, runs+1):
        folder = f"{base}_{i}"
        os.makedirs(folder, exist_ok=True)
        # copy all .inp and .T7 files
        for ext in ['*.inp','*.T7','SAVECO*']:
            for f in glob.glob(ext): shutil.copy(f, folder)
        # apply perturbation
        pert = apply_perturbations(input_filename, params, folder)
        args.append((folder, pert))

    # run in parallel
    print(f"Starting {runs} PARMELA runs in parallel...")
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(run_case, args))
    
    print("All PARMELA runs completed.")
    
    # collect
    error_dat_path = aggregate_results(results, params, run_id)

    # PDF summary
    generate_error_summary_pdf(params, error_dat_path, run_id)
    
    # Plotting
    orbit_figure(base, runs, run_id)

