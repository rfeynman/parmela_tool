#To evaluate the preinjector errors.
#This is pal version, that can run the multiple parmela at once. Limited by the # of CPU
#author: Erdong Wang
#Version 1.3 Jun.24th 2025
# useage:
#   python error_ana.py rr6.inp error.yaml 

import os
import sys
import shutil
import yaml
import glob
import uuid
import subprocess
from concurrent.futures import ProcessPoolExecutor

# 1. Find element indices and main frequency
def find_ele_ind(lines):
    elements = {"quad": [], "solenoid": [], "cell": [], "trwave": []}
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
    # quadrupole and solenoid use power supply amp
    dist["quad"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["quad"]))
    dist["solenoid"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["solenoid"]))

    # cell
    phase_cell = truncated_normal(params["cell_rf_phase_mean"], params["cell_rf_phase_sig"], params["cell_rf_phase_bound"], len(elements["cell"]))
    amp_cell = truncated_normal(params["cell_rf_amp_mean"], params["cell_rf_amp_sig"], params["cell_rf_amp_bound"], len(elements["cell"]))
    dist["cell"] = list(zip(phase_cell, amp_cell))

    # trwave
    phase_trwave = truncated_normal(params["trwave_rf_phase_mean"], params["trwave_rf_phase_sig"], params["trwave_rf_phase_bound"], len(elements["trwave"]))
    amp_trwave = truncated_normal(params["trwave_rf_amp_mean"], params["trwave_rf_amp_sig"], params["trwave_rf_amp_bound"], len(elements["trwave"]))
    dist["trwave"] = list(zip(phase_trwave, amp_trwave))

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
        subprocess.run(["parmela", os.path.basename(pert_file)], check=True)
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
    return os.path.join(folder, "TIMESTEPEMITTANCE.TBL")

# 6. Aggregate results
def aggregate_results(tbl_paths, params, run_id):
    # prepare header in main analysis file
    main_file = "error_analysis_dat.txt"
    header = "T(deg) Z(cm) Xun(mm-mrad) Yun(mm-mrad) Zun(mm-mrad) Xn(mm-mrad) Yn(mm-mrad) Zn(mm-mrad) Xrms(mm) Yrms(mm) Zrmz(mm) <kE>(MeV) Del-Erms <X>(mm) <Xpn>(mrad) <Y>(mm) <Ypn>(mrad) <Z>(cm) <Zpn>(rad) EZref(MV/m)"
    with open(main_file, 'w') as out:
        out.write(header + '\n')
    # append each -3 line
    for tbl in tbl_paths:
        with open(tbl, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 3:
                line = lines[-3].strip()
                with open(main_file, 'a') as out:
                    out.write(line + '\n')
    # rename
    shutil.move(main_file, f"error_analysis_dat_{run_id}.txt")

# 7. Main entry
if __name__ == "__main__":
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
        for ext in ['*.inp','*.T7']:
            for f in glob.glob(ext): shutil.copy(f, folder)
        # apply perturbation
        pert = apply_perturbations(input_filename, params, folder)
        args.append((folder, pert))

    # run in parallel
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(run_case, args))
    # collect
    aggregate_results(results, params, run_id)
