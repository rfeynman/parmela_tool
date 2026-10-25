#To evaluate the preinjector errors.
#This is slow version, that can run the single parmela.
#author: Erdong Wang
#Version 0.8 Jun.3rd 2025
# useage:
#   python error_ana.py rr6.inp error.yaml

import os
import sys
import shutil
import yaml
import numpy as np
import subprocess
from scipy.stats import truncnorm

# 1. Read and copy the input file
def read_and_copy_input(filename):
    base, ext = os.path.splitext(filename)
    new_filename = f"{base}_erranaly{ext}"
    shutil.copyfile(filename, new_filename)
    return new_filename

# 2. Find element line indices
def find_ele_ind(lines):
    elements = {"quad": [], "solenoid": [], "cell": [], "trwave": []}
    mainfreq = None
    i = 0
    while i < len(lines):
        line = lines[i].strip().lower()
        if line.startswith("run"):
            parts = lines[i].split()
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
        i += 1
    return elements, mainfreq

# Helper: truncated normal
def truncated_normal(mean, sigma, bound, size):
    if sigma == 0:
        return np.full(size, mean)
    a, b = (-bound - mean) / sigma, (bound - mean) / sigma
    return truncnorm.rvs(a, b, loc=mean, scale=sigma, size=size)

# 3. Random seed generator
def randseed(elements, params):
    dist = {}
    dist["quad"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["quad"]))
    dist["solenoid"] = truncated_normal(params["ps_amp_mean"], params["ps_amp_sig"], params["ps_amp_bound"], len(elements["solenoid"]))

    phase_cell = truncated_normal(params["cell_rf_phase_mean"], params["cell_rf_phase_sig"], params["cell_rf_phase_bound"], len(elements["cell"]))
    amp_cell = truncated_normal(params["cell_rf_amp_mean"], params["cell_rf_amp_sig"], params["cell_rf_amp_bound"], len(elements["cell"]))
    dist["cell"] = list(zip(phase_cell, amp_cell))

    phase_trwave = truncated_normal(params["trwave_rf_phase_mean"], params["trwave_rf_phase_sig"], params["trwave_rf_phase_bound"], len(elements["trwave"]))
    amp_trwave = truncated_normal(params["trwave_rf_amp_mean"], params["trwave_rf_amp_sig"], params["trwave_rf_amp_bound"], len(elements["trwave"]))
    dist["trwave"] = list(zip(phase_trwave, amp_trwave)) 

    print("Generated random distributions:")
    for k, v in dist.items():
        print(f"{k}: {v}")

    return dist

# 4. Modify trwave lines and 84 lines following each

def trave_cngele(lines, elements, dist, mainfreq):
    for i, (p, a) in zip(elements["trwave"], dist["trwave"]):
        for j in range(i, i + 86):
            if j >= len(lines):
                break
            parts = lines[j].split()
            if len(parts) >= 10:
                if p != 0:
                    parts[4] = str(float(parts[4]) + p)
                if a != 0:
                    parts[5] = str(float(parts[5]) * (1 + a))
                lines[j] = ' '.join(parts) + '\n'

# 5. Modify cell only

def cell_cngele(lines, elements, dist, mainfreq):
    for i, (p, a) in zip(elements["cell"], dist["cell"]):
        parts = lines[i].split()
        if len(parts) >= 10:
            if p != 0:
                print("cell original phase:", parts[4])
                parts[4] = str(float(parts[4]) + p * mainfreq / float(parts[9]))
                print("cell modified phase:", parts[4])
            if a != 0:
                print("cell original amp:", parts[5])
                parts[5] = str(float(parts[5]) * (1 + a))
                print("cell modified amp:", parts[5])
            lines[i] = ' '.join(parts) + '\n'

# 6. Modify solenoid

def solenoid_cngele(lines, elements, dist):
    for i, a in zip(elements["solenoid"], dist["solenoid"]):
        if a != 0:
            parts = lines[i].split()
            print("solenoid original amp:", parts[4])
            parts[4] = str(float(parts[4]) * (1 + a))
            print("solenoid modified amp:", parts[4])
            lines[i] = ' '.join(parts) + '\n'

# 7. Modify quad

def quad_cngele(lines, elements, dist):
    for i, a in zip(elements["quad"], dist["quad"]):
        if a != 0:
            parts = lines[i].split()
            print("quad original amp:", parts[4])
            parts[4] = str(float(parts[4]) * (1 + a))
            print("quad modified amp:", parts[4])
            lines[i] = ' '.join(parts) + '\n'

# 8. Run Parmela
def run_parmela(filename):
    try:
        subprocess.run(["parmela", filename], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running parmela: {e}")

# 9. Analyze results
def analysis(ps_amp_sig, ps_amp_mean,
             cell_rf_phase_sig, cell_rf_phase_mean, cell_rf_phase_bound,
             cell_rf_amp_sig, cell_rf_amp_mean, cell_rf_amp_bound,
             trwave_rf_phase_sig, trwave_rf_phase_mean, trwave_rf_phase_bound,
             trwave_rf_amp_sig, trwave_rf_amp_mean, trwave_rf_amp_bound):
    output_file = f"TIMESTEPEMITTANCE_cellrfphase{cell_rf_phase_sig}_cellrfamp{cell_rf_amp_sig}_trwrfphase{trwave_rf_phase_sig}_trwrfamp{trwave_rf_amp_sig}_ps{ps_amp_sig}.TBL"
    if os.path.exists("TIMESTEPEMITTANCE.TBL"):
        shutil.copyfile("TIMESTEPEMITTANCE.TBL", output_file)
        with open(output_file, 'r') as f:
            lines = f.readlines()
            result_line = lines[-3].strip()
        if not os.path.exists("error_analysis_dat.txt"):
            with open("error_analysis_dat.txt", 'w') as out:
                header = "T(deg)         Z(cm)        Xun(mm-mrad)   Yun(mm-mrad)   Zun(mm-mrad)  Xn(mm-mrad)    Yn(mm-mrad)    Zn(mm-mrad)     Xrms(mm)       Yrms(mm)      Zrmz(mm)       <kE>(MeV)      Del-Erms       <X>(mm)         <Xpn>(mrad)   <Y>(mm)        <Ypn>(mrad)  <Z>(cm)       <Zpn>(rad)   EZref(MV/m)   ps_amp_sig ps_amp_mean cell_rf_phase_sig cell_rf_phase_mean cell_rf_phase_bound cell_rf_amp_sig cell_rf_amp_mean cell_rf_amp_bound trwave_rf_phase_sig trwave_rf_phase_mean trwave_rf_phase_bound trwave_rf_amp_sig trwave_rf_amp_mean trwave_rf_amp_bound \n"
                out.write(header)
        with open("error_analysis_dat.txt", 'a') as out:
            out.write(f"{result_line} {ps_amp_sig} {ps_amp_mean} {cell_rf_phase_sig} {cell_rf_phase_mean} {cell_rf_phase_bound} {cell_rf_amp_sig} {cell_rf_amp_mean} {cell_rf_amp_bound} {trwave_rf_phase_sig} {trwave_rf_phase_mean} {trwave_rf_phase_bound} {trwave_rf_amp_sig} {trwave_rf_amp_mean} {trwave_rf_amp_bound} \n")

def main():
    input_filename = sys.argv[1]
    yaml_filename = sys.argv[2]
    with open(yaml_filename, 'r') as f:
        params = yaml.safe_load(f)

    base_yaml = os.path.splitext(os.path.basename(yaml_filename))[0]
    import uuid

    run_id = uuid.uuid4().hex[:8]
    shutil.copyfile(yaml_filename, f"error_{run_id}.yaml")
    final_analysis_file = f"error_analysis_dat_{run_id}.txt"
    for i in range(params.get("runs", 1)):
        copied_file = read_and_copy_input(input_filename)

        with open(copied_file, 'r') as f:
            lines = f.readlines()

        elements, mainfreq = find_ele_ind(lines)
        dists = randseed(elements, params)

        trave_cngele(lines, elements, dists, mainfreq)
        cell_cngele(lines, elements, dists, mainfreq)
        solenoid_cngele(lines, elements, dists)
        quad_cngele(lines, elements, dists)

        with open(copied_file, 'w') as f:
            f.writelines(lines)

        run_parmela(copied_file)
        analysis(
                 params["ps_amp_sig"], params["ps_amp_mean"],
                 params["cell_rf_phase_sig"], params["cell_rf_phase_mean"], params["cell_rf_phase_bound"],
                 params["cell_rf_amp_sig"], params["cell_rf_amp_mean"], params["cell_rf_amp_bound"],
                 params["trwave_rf_phase_sig"], params["trwave_rf_phase_mean"], params["trwave_rf_phase_bound"],
                 params["trwave_rf_amp_sig"], params["trwave_rf_amp_mean"], params["trwave_rf_amp_bound"]
        )


    shutil.move("error_analysis_dat.txt", final_analysis_file)

if __name__ == "__main__":
    main()
