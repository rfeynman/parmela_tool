import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import re
import math

# ==========================================
# 1. Configuration
# ==========================================
# Set the directory path where your file is located.
file_path = "/Users/wange/Documents/Research/eRHIC injector/eRHIC baseline/Beamline/Lattice/sband2025/july"
lattice_file_name = "OUTPAR_sep.TXT"
beam_file_name = "TIMESTEPEMITTANCE.tbl"

full_lattice_path = os.path.join(file_path, lattice_file_name)
full_beam_path = os.path.join(file_path, beam_file_name)

# Output filenames
output_excel = os.path.join(file_path, f"lattice_{lattice_file_name.replace('.TXT', '')}.xlsx")
output_plot = os.path.join(file_path, f"layout_{lattice_file_name.replace('.TXT', '')}.png")

# --- Beam Plotting Options ---
beam = 1  # 1: Plot beam parameters, 0: No beam plot
# Select columns to plot from TIMESTEPEMITTANCE.tbl headers (exact strings)
# Available options from file header:
# 'T(deg)', 'Z(cm)', 'Xun', 'Yun', 'Zun', 'Xn', 'Yn', 'Zn', 
# 'Xrms(mm)', 'Yrms(mm)', 'Zrms(mm)', 'kE(MeV)', 'Del-kE(MeV)', 
# '<X>(mm)', '<Xpn>(mrad)', '<Y>(mm)', '<Ypn>(mrad)', '<Z>(cm)', 
# '<Zpn>(mrad)', 'Ezref(MV/m)'
# Example: ['Xrms(mm)', 'Yrms(mm)'] or ['kE(MeV)']
beam_y_keys = ['Xrms(mm)', 'Yrms(mm)'] 

# ==========================================
# 2. Data Processing Functions
# ==========================================

def process_parmela_file(filepath):
    if not os.path.exists(filepath):
        print(f"Error: File not found at {filepath}")
        return None, None

    with open(filepath, 'r') as f:
        lines = f.readlines()

    req_data = {
        'Solenoid': [], 'cell': [], 'quad': [],
        'steerer': [], 'bend': [], 'trwave': []
    }
    
    lattice_rows = []
    bend_angles = [] 
    current_trwave_group = None
    start_def = False
    
    # --- 1. Parse Definitions ---
    for line in lines:
        raw_line = line.strip()
        clean_line = re.sub(r'[\']', '', raw_line).strip()
        
        if not clean_line: continue
        
        # Detect sections
        if 'title' in clean_line.lower():
            start_def = True
            continue
        if 'error' in clean_line.lower() and start_def:
            if current_trwave_group:
                 req_data['trwave'].append(current_trwave_group)
                 current_trwave_group = None
            break 
        
        if not start_def: continue
        
        # Skip comments explicitly
        if clean_line.startswith('!') or clean_line.startswith(';'):
            continue

        parts = re.split(r'[\s!;]+', clean_line)
        parts = [p for p in parts if p]
        
        if not parts: continue

        element_name = parts[0].lower()
        
        # TRWAVE Logic
        if element_name == 'trwave':
            try:
                length = float(parts[1])
                phase = parts[4] if len(parts) > 4 else ''
                amp = parts[5] if len(parts) > 5 else ''
                
                if current_trwave_group is None:
                    current_trwave_group = {
                        'Element': parts[0], 'Length [cm]': length,
                        'Phase [deg]': phase, 'Amplitude': amp
                    }
                else:
                    current_trwave_group['Length [cm]'] += length
            except: pass
            continue
        else:
            if current_trwave_group:
                req_data['trwave'].append(current_trwave_group)
                current_trwave_group = None

        # Standard Elements
        try:
            if element_name == 'solenoid':
                req_data['Solenoid'].append({
                    'Element': parts[0], 'Length [cm]': parts[1], 
                    'Aperture [cm]': parts[2], 'Amp1': parts[4]
                })
            elif element_name == 'cell':
                req_data['cell'].append({
                    'Element': parts[0], 'Length [cm]': parts[1], 
                    'Phase': parts[4], 'Amp1': parts[5]
                })
            elif element_name == 'quad':
                amp_idx = 5 if len(parts) > 5 else 4
                req_data['quad'].append({
                    'Element': parts[0], 'Length [cm]': parts[1], 
                    'Aperture [cm]': parts[2] if len(parts)>2 else '', 
                    'Amp1': parts[amp_idx] if len(parts)>amp_idx else ''
                })
            elif element_name == 'steerer':
                req_data['steerer'].append({
                    'Element': parts[0], 'Length [cm]': parts[1], 
                    'Aperture [cm]': parts[2], 
                    'Amp1': parts[4] if len(parts)>4 else '', 
                    'Amp2': parts[5] if len(parts)>5 else ''
                })
            elif element_name == 'bend':
                angle = parts[5] if len(parts)>5 else '0'
                req_data['bend'].append({
                    'Element': parts[0], 'Length [cm]': parts[1], 
                    'Gap [cm]': parts[2], 'Angle[deg]': angle
                })
                bend_angles.append(float(angle))
        except (IndexError, ValueError):
            continue

    # --- 2. Parse Lattice Table ---
    start_table = False
    last_ele = None
    bend_index_counter = 0
    
    for line in lines:
        clean_line = re.sub(r'[\']', '', line).strip()
        parts = list(filter(None, re.split(r'\s+', clean_line)))
        
        if not parts: continue
        
        if 'n' in parts and 'z1' in parts and 'element' in parts:
            start_table = True
            continue
        if 'zlimit' in clean_line.lower():
            if last_ele: lattice_rows.append(last_ele)
            break
        if not start_table: continue
        if len(parts) < 4 or not parts[0].isdigit(): continue
        
        try:
            n = int(parts[0])
            z1 = float(parts[1])
            name = parts[2]
            z2 = float(parts[3])
            phase = float(parts[5]) if len(parts) > 5 else 0.0
            amp = float(parts[6]) if len(parts) > 6 else 0.0
            angle = 0.0
            
            if 'bend' in name.lower():
                if bend_index_counter < len(bend_angles):
                    angle = bend_angles[bend_index_counter]
                    bend_index_counter += 1
            
            current_data = {
                'n': n, 'z1': z1, 'element': name, 'z2': z2, 
                'dz': z2 - z1, 'phase': phase, 'amp': amp, 'angle': angle
            }
            
            if last_ele is None:
                last_ele = current_data
            else:
                if last_ele['element'] == name:
                    last_ele['z2'] = z2
                    last_ele['dz'] = last_ele['z2'] - last_ele['z1']
                    last_ele['n'] = n 
                else:
                    lattice_rows.append(last_ele)
                    last_ele = current_data
        except ValueError:
            continue

    if last_ele and (not lattice_rows or lattice_rows[-1] != last_ele):
        lattice_rows.append(last_ele)

    return req_data, lattice_rows

def process_beam_file(filepath):
    """
    Parses the TIMESTEPEMITTANCE.tbl file.
    """
    if not os.path.exists(filepath):
        print(f"Error: Beam file not found at {filepath}")
        return None

    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    headers = []
    data_start_idx = -1
    
    # 1. Find Headers (TITLES section)
    in_titles = False
    for i, line in enumerate(lines):
        clean = line.strip()
        if clean == "TITLES":
            in_titles = True
            continue
        if clean == "ENDTITLES":
            in_titles = False
            continue
        if in_titles:
            headers.append(clean)
        
        # 2. Find Start of Data
        if clean == "DATA":
            data_start_idx = i + 1
            break
            
    if data_start_idx == -1 or not headers:
        print("Error: Could not find TITLES or DATA block in beam file.")
        return None

    # 3. Read Data
    data_lines = lines[data_start_idx:]
    # Filter out lines that don't look like data (e.g. comments or empty)
    valid_data = []
    for line in data_lines:
        parts = line.strip().split()
        # Basic check: if number of columns matches headers
        if len(parts) == len(headers):
             valid_data.append(parts)
        elif len(parts) > 0 and parts[0].replace('.','',1).isdigit():
             # Try to match if length mismatch is due to spacing, though strict matching is safer
             # For now assume strict matching
             if len(parts) == len(headers):
                 valid_data.append(parts)

    try:
        df = pd.DataFrame(valid_data, columns=headers, dtype=float)
        return df
    except Exception as e:
        print(f"Error creating beam DataFrame: {e}")
        return None


# ==========================================
# 3. Geometric Figure Generation
# ==========================================

def generate_lattice_plot(df_lattice, save_path, df_beam=None, beam_keys=None):
    
    # Determine Layout based on beam flag
    if df_beam is not None and beam_keys:
        nrows = 3
        figsize = (36, 18) # Increased height slightly to accommodate equal subplots comfortably
        # All three figures have the same relative height (1:1:1)
        gridspec_kw = {'height_ratios': [1, 1, 1]} 
    else:
        nrows = 2
        figsize = (36, 12)
        gridspec_kw = {'height_ratios': [1, 1]}

    fig, axes = plt.subplots(nrows, 1, figsize=figsize, gridspec_kw=gridspec_kw, sharex=False)
    plt.subplots_adjust(hspace=0.3)
    
    # Assign Axes
    ax1 = axes[0] # Schematic
    if nrows == 3:
        ax_beam = axes[1]
        ax2 = axes[2] # Top View
    else:
        ax2 = axes[1] # Top View

    # Common Style Settings
    lw_box = 2.5
    lw_drift = 2.0
    font_size = 16
    
    # ===========================
    # Plot 1: Linear Schematic (Up Plot)
    # ===========================
    legend_patches1 = {}
    max_s = 0
    
    for idx, row in df_lattice.iterrows():
        s_start = row['z1'] / 100.0
        s_end = row['z2'] / 100.0
        length_m = s_end - s_start
        name = row['element'].lower()
        max_s = max(max_s, s_end)
        
        if 'drift' in name:
            line, = ax1.plot([s_start, s_end], [0, 0], color='black', linewidth=lw_drift)
            if 'Drift' not in legend_patches1: legend_patches1['Drift'] = line

        elif 'cathode' in name:
            circle = patches.Circle((s_start, 0), radius=0.1, edgecolor='tab:orange', facecolor='none', linewidth=lw_box)
            ax1.add_patch(circle)
            legend_patches1['Cathode'] = circle
            
        elif 'bend' in name:
            rect = patches.Rectangle((s_start, 0), length_m, 0.15, facecolor='none', edgecolor='tab:red', linewidth=lw_box)
            ax1.add_patch(rect)
            legend_patches1['Bend'] = rect

        elif 'solenoid' in name:
            rect = patches.Rectangle((s_start, -0.1), length_m, 0.2, edgecolor='tab:cyan', facecolor='none', linewidth=lw_box)
            ax1.add_patch(rect)
            legend_patches1['Solenoid'] = rect
            
        elif 'quad' in name:
            # Changed 'tab:blue' to 'blue' as requested
            rect = patches.Rectangle((s_start, -0.1), length_m, 0.2, edgecolor='blue', facecolor='none', linewidth=lw_box)
            ax1.add_patch(rect)
            legend_patches1['Quad'] = rect
            
        elif 'cell' in name or 'trwave' in name:
            h = 0.1
            y_start = -h / 2.0
            rect = patches.Rectangle((s_start, y_start), length_m, h, facecolor='none', edgecolor='tab:green', linewidth=lw_box)
            ax1.add_patch(rect)
            ax1.plot([s_start, s_end], [y_start, y_start + h], color='tab:green', linestyle='--', linewidth=1.5)
            ax1.plot([s_start, s_end], [y_start + h, y_start], color='tab:green', linestyle='--', linewidth=1.5)
            legend_patches1['Cell/Trwave'] = rect
            
        elif 'steerer' in name:
            st_start = s_start - 0.15
            st_len = (s_end + 0.02) - st_start
            rect = patches.Rectangle((st_start, -0.075), st_len, 0.15, facecolor='none', edgecolor='tab:purple', linewidth=lw_box)
            ax1.add_patch(rect)
            legend_patches1['Steerer'] = rect

    ax1.set_xlabel("S [m]", fontsize=font_size) # Added X label back for Schematic
    ax1.set_ylabel("Y [m]", fontsize=font_size)
    ax1.set_title("Schematic Layout (S vs Y)", fontsize=font_size+2)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.grid(True, linestyle='--', color='darkgray', alpha=0.7)
    ax1.set_xlim(-1.0, max_s + 1.0)
    ax1.legend(legend_patches1.values(), legend_patches1.keys(), loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=14)

    # ===========================
    # Plot 2: Beam Parameters (Middle Plot - Optional)
    # ===========================
    if nrows == 3:
        # X-axis for beam is Z(cm) -> convert to meters
        # Assuming column 'Z(cm)' exists based on PARMELA standard
        if 'Z(cm)' in df_beam.columns:
            beam_x = df_beam['Z(cm)'] / 100.0
        elif '<Z>(cm)' in df_beam.columns:
             beam_x = df_beam['<Z>(cm)'] / 100.0
        else:
            # Fallback to index or guess
            beam_x = df_beam.iloc[:, 1] / 100.0 # Assuming Z is 2nd col
            
        # Changed 'tab:blue' to 'blue' in color cycle
        colors = ['tab:red', 'blue', 'tab:orange', 'tab:green']
        
        # Plot first parameter on primary Y axis
        if len(beam_keys) > 0:
            key1 = beam_keys[0]
            if key1 in df_beam.columns:
                ax_beam.plot(beam_x, df_beam[key1], color=colors[0], linewidth=2, label=key1)
                ax_beam.set_ylabel(key1, color=colors[0], fontsize=font_size)
                ax_beam.tick_params(axis='y', labelcolor=colors[0], labelsize=12)
        
        # Plot second parameter on secondary Y axis if requested
        if len(beam_keys) > 1:
            key2 = beam_keys[1]
            if key2 in df_beam.columns:
                ax_beam2 = ax_beam.twinx()
                ax_beam2.plot(beam_x, df_beam[key2], color=colors[1], linewidth=2, label=key2)
                ax_beam2.set_ylabel(key2, color=colors[1], fontsize=font_size)
                ax_beam2.tick_params(axis='y', labelcolor=colors[1], labelsize=12)
        
        ax_beam.set_xlabel("S [m]", fontsize=font_size)
        ax_beam.set_xlim(-1.0, max_s + 1.0) # Align with top plot
        ax_beam.grid(True, linestyle='--', color='darkgray', alpha=0.7)
        ax_beam.set_title("Beam Parameters", fontsize=font_size+2)
        # Increase tick label size for beam plot x-axis
        ax_beam.tick_params(axis='x', labelsize=12)


    # ===========================
    # Plot 3: Top View Layout (Bottom Plot)
    # ===========================
    legend_patches2 = {}
    
    total_bend_deg = df_lattice[df_lattice['element'].str.lower().str.contains('bend')]['angle'].sum()
    cur_angle = -np.radians(total_bend_deg)
    
    cur_pos = np.array([0.0, 0.0]) # [Global Z, Global X]
    
    for idx, row in df_lattice.iterrows():
        length_m = row['dz'] / 100.0
        name = row['element'].lower()
        p1 = cur_pos.copy()
        
        if 'bend' in name:
            bend_deg = row['angle']
            bend_rad = np.radians(bend_deg)
            
            if abs(bend_rad) > 1e-9:
                R = length_m / bend_rad
                cx = p1[0] - R * np.sin(cur_angle)
                cy = p1[1] + R * np.cos(cur_angle)
                
                new_angle = cur_angle + bend_rad
                p2_x = cx + R * np.sin(new_angle)
                p2_y = cy - R * np.cos(new_angle)
                p2 = np.array([p2_x, p2_y])
                
                thetas = np.linspace(cur_angle, new_angle, 20)
                arc_x = cx + R * np.sin(thetas)
                arc_y = cy - R * np.cos(thetas)
                
                ax2.plot(arc_x, arc_y, color='tab:red', linewidth=lw_box)
                
                width = 0.15
                arc_x_in = cx + (R - width/2) * np.sin(thetas)
                arc_y_in = cy - (R - width/2) * np.cos(thetas)
                arc_x_out = cx + (R + width/2) * np.sin(thetas)
                arc_y_out = cy - (R + width/2) * np.cos(thetas)
                poly_x = np.concatenate([arc_x_in, arc_x_out[::-1]])
                poly_y = np.concatenate([arc_y_in, arc_y_out[::-1]])
                poly_xy = np.column_stack((poly_x, poly_y))
                
                poly = patches.Polygon(poly_xy, closed=True, facecolor='none', edgecolor='tab:red', linewidth=lw_box)
                ax2.add_patch(poly)
                
                cur_pos = p2
                cur_angle = new_angle
                if 'Bend' not in legend_patches2:
                     legend_patches2['Bend'] = patches.Patch(facecolor='none', edgecolor='tab:red', label='Bend')
            else:
                step = np.array([np.cos(cur_angle), np.sin(cur_angle)]) * length_m
                cur_pos = p1 + step
                ax2.plot([p1[0], cur_pos[0]], [p1[1], cur_pos[1]], color='tab:red', linewidth=lw_box)

        else:
            step = np.array([np.cos(cur_angle), np.sin(cur_angle)]) * length_m
            p2 = p1 + step
            
            if 'drift' in name:
                ax2.plot([p1[0], p2[0]], [p1[1], p2[1]], color='black', linewidth=lw_drift)
            elif 'cathode' in name:
                circle = patches.Circle(p1, radius=0.1, edgecolor='tab:orange', facecolor='none', linewidth=lw_box)
                ax2.add_patch(circle)
                if 'Cathode' not in legend_patches2:
                    legend_patches2['Cathode'] = patches.Patch(facecolor='none', edgecolor='tab:orange', label='Cathode')
            else:
                edge_color = 'gray'
                width = 0.2
                
                if 'quad' in name: edge_color = 'blue'; width=0.2
                elif 'solenoid' in name: edge_color = 'tab:cyan'; width=0.2
                elif 'cell' in name or 'trwave' in name: edge_color = 'tab:green'; width=0.1
                elif 'steerer' in name: edge_color = 'tab:purple'; width=0.15
                
                draw_length = length_m
                draw_p1 = p1
                
                if 'steerer' in name:
                    direction = np.array([np.cos(cur_angle), np.sin(cur_angle)])
                    draw_p1 = p1 - direction * 0.15
                    draw_length = length_m + 0.17
                
                perp_vec = np.array([-np.sin(cur_angle), np.cos(cur_angle)]) * (-width/2.0)
                p_anchor = draw_p1 + perp_vec
                
                rect = patches.Rectangle(
                    p_anchor, draw_length, width, 
                    angle=np.degrees(cur_angle), 
                    facecolor='none', edgecolor=edge_color, linewidth=lw_box
                )
                ax2.add_patch(rect)
                
                if 'cell' in name or 'trwave' in name:
                    v_len = np.array([np.cos(cur_angle), np.sin(cur_angle)]) * draw_length
                    v_wid = np.array([-np.sin(cur_angle), np.cos(cur_angle)]) * width
                    c1 = p_anchor 
                    c2 = p_anchor + v_len 
                    c3 = p_anchor + v_len + v_wid 
                    c4 = p_anchor + v_wid 
                    
                    ax2.plot([c1[0], c3[0]], [c1[1], c3[1]], color=edge_color, linestyle='--', linewidth=1.5)
                    ax2.plot([c4[0], c2[0]], [c4[1], c2[1]], color=edge_color, linestyle='--', linewidth=1.5)
                
                label = name.capitalize()
                if 'cell' in name or 'trwave' in name: label = "Cell/Trwave"
                if label not in legend_patches2:
                    legend_patches2[label] = patches.Patch(facecolor='none', edgecolor=edge_color, label=label)

            cur_pos = p2

    ax2.set_xlabel("Global Z [m]", fontsize=font_size)
    ax2.set_ylabel("Global X [m]", fontsize=font_size)
    ax2.set_title("Top-View layout", fontsize=font_size+2)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.grid(True, linestyle='--', color='darkgray', alpha=0.7)
    ax2.legend(legend_patches2.values(), legend_patches2.keys(), loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=14)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved: {save_path}")
    plt.show()

# ==========================================
# 4. Main Execution
# ==========================================

if __name__ == "__main__":
    print(f"Processing Lattice: {full_lattice_path}")
    req_data, lattice_data = process_parmela_file(full_lattice_path)
    
    df_beam = None
    if beam == 1:
        print(f"Processing Beam: {full_beam_path}")
        df_beam = process_beam_file(full_beam_path)

    if req_data and lattice_data:
        df_lattice = pd.DataFrame(lattice_data)
        df_lattice['n'] = range(1, len(df_lattice) + 1)

        print(f"Creating Excel: {output_excel}")
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            df_lattice.to_excel(writer, sheet_name='lattice', index=False)
            
            row_ptr = 0
            for key in ['Solenoid', 'cell', 'quad', 'steerer', 'bend', 'trwave']:
                data = req_data.get(key, [])
                if data:
                    pd.DataFrame(data).to_excel(writer, sheet_name='requirement', startrow=row_ptr, index=False)
                    row_ptr += len(data) + 3

        print(f"Creating Plot: {output_plot}")
        generate_lattice_plot(df_lattice, output_plot, df_beam, beam_y_keys)
        print("Done.")
    else:
        print("Failed to process lattice data.")