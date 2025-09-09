#To evaluate the preinjector errors.
#This is pal version, that can run the multiple parmela at once. Limited by the # of CPU
#author: Erdong Wang
#Version 1.4 Sept.8th 2025


import os
import sys
import shutil
import yaml
import glob
import uuid
import subprocess
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import matplotlib.pyplot as plt


# 1. Find element indices and main frequency
def find_ele_ind(lines):
    elements = {"quad": [], "solenoid": [], "cell": [], "trwave": [],"bend": [],"steerer": []}
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
        elif line.startswith("BEND"):
            elements["bend"].append(i)
        elif line.startswith("steerer"):
            elements["steerer"].append(i)
        elif line.startswith("cell"):
            if i + 1 < len(lines) and lines[i + 1].strip().lower().startswith("trwave"):
                elements["trwave"].append(i)
            else:
                elements["cell"].append(i)
    return elements, mainfreq

# 2. Truncated normal helper
def truncated_normal(mean, sigma, bound, size):
    from scipy.stats import truncnorm
    import numpy as np
    if sigma == 0:
        return np.full(size, mean)
    a, b = (-bound - mean) / sigma, (bound - mean) / sigma
    return truncnorm.rvs(a, b, loc=mean, scale=sigma, size=size)

# 3. Generate distributions
def randseed(elements, params):
    dist = {}
    # quadrupole ,solenoid and steerer use power supply amp
    dist["quad"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["quad"]))
    dist["solenoid"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["solenoid"]))
    dist["steerer"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["solenoid"]))

    # cell
    phase_cell = truncated_normal(params["cell_rf_phase_mean"], params["cell_rf_phase_sig"], params["cell_rf_phase_bound"], len(elements["cell"]))
    amp_cell = truncated_normal(params["cell_rf_amp_mean"], params["cell_rf_amp_sig"], params["cell_rf_amp_bound"], len(elements["cell"]))
    dist["cell"] = list(zip(phase_cell, amp_cell))

    # trwave
    phase_trwave = truncated_normal(params["trwave_rf_phase_mean"], params["trwave_rf_phase_sig"], params["trwave_rf_phase_bound"], len(elements["trwave"]))
    amp_trwave = truncated_normal(params["trwave_rf_amp_mean"], params["trwave_rf_amp_sig"], params["trwave_rf_amp_bound"], len(elements["trwave"]))
    dist["trwave"] = list(zip(phase_trwave, amp_trwave))
    
    
    #bend
    dist["bend"]= truncated_normal(params["bend_amp_mean"], params["bend_amp_sig"], params["bend_amp_bound"], len(elements["bend"]))

    return dist

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
    # modify cell
    for idx, (p, a) in zip(elements["cell"], dist["cell"]):
        parts = lines[idx].split()
        if len(parts) >= 6:
            if p != 0:
                parts[4] = str(float(parts[4]) + p)
            if a != 0:
                parts[5] = str(float(parts[5]) * (1 + a))
            lines[idx] = ' '.join(parts) + '\n'
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
    shutil.move(main_file, f"error_analysis_dat_{run_id}.txt")


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
            
        df = pd.read_csv(tbl_path, skiprows=line_ini + 2, sep='\s+', header=None, engine='python')
        
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
    aggregate_results(results, params, run_id)
    
    # Plotting
    orbit_figure(base, runs, run_id)

