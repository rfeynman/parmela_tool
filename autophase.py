#!/usr/bin/env python3
import sys
import subprocess
import numpy as np
import os

# Usage:
#   python optimize_parmela.py <input_file.inp> <mode> [params...]
# Modes:
#   g: gradient descent with adaptive step size
#     params: <init_step> [sect] [max_iters] [lr]
#   p: parabola fitting
#     params: <delphase> [sect] [iterations]
# e.g.:python3 autophase.py rr6.inp p 10 13 10 ; inital delta phase is 10 degree, section is 14 (start from 0), maximum iterations is 10
# Convergence and optimization settings
tol = 2e-7  # dE convergence tolerance
momentum = 0.9  # momentum for gradient descent
step_up = 1.2   # factor to increase step when improving
step_down = 0.5 # factor to decrease step when not improving

def print_usage_and_exit():
    print("Usage: python optimize_parmela.py <input_file.inp> <g|p> [...]")
    sys.exit(1)

# Validate args
if len(sys.argv) < 3:
    print_usage_and_exit()

# Parse basic args
base_inp = sys.argv[1]
mode = sys.argv[2].lower()
root, ext = os.path.splitext(base_inp)
if ext.lower() != ".inp":
    print_usage_and_exit()
default_temp = f"{root}_temp{ext}"
default_tbl = "TIMESTEPEMITTANCE.TBL"

# Prepare temp file (fresh copy)
if os.path.exists(default_temp):
    os.remove(default_temp)
subprocess.run(["cp", base_inp, default_temp], check=True)

# Helpers
def find_indices():
    idxs, counts = [], []
    lines = open(default_temp).readlines()
    for i, L in enumerate(lines):
        if L.lower().startswith("cell") and i+1 < len(lines) and lines[i+1].lower().startswith("trwave"):
            c = 0
            for j in range(i+1, len(lines)):
                if lines[j].lower().startswith("trwave"): c += 1
                else: break
            idxs.append(i)
            counts.append(c)
    return idxs, counts


def cngele(newphase, sect=0):
    idxs, counts = find_indices()
    i, n = idxs[sect], counts[sect]
    lines = open(default_temp).read().splitlines()
    parts = lines[i].split(); parts[4] = str(newphase + 90); lines[i] = " ".join(parts)
    for k in range(n):
        parts = lines[i+1+k].split(); parts[4] = str(newphase); lines[i+1+k] = " ".join(parts)
    open(default_temp, 'w').write("\n".join(lines) + "\n")


def parse_delE():
    lines = open(default_tbl).read().splitlines()
    cols = lines[-3].split()
    return float(cols[12])  # 13th column (0-based)


def run_parmela():
    subprocess.run(["parmela", default_temp], check=True)

# Gradient Descent Optimizer w/ adaptive step
def optimize_g(init_step, sect, max_it, lr):
    diff_step = init_step
    idxs, _ = find_indices()
    phase = float(open(default_temp).read().splitlines()[idxs[sect]].split()[4]) - 90.0
    vel = 0.0
    prevE = None

    for it in range(1, max_it + 1):
        # finite-difference gradient
        cngele(phase + diff_step, sect); run_parmela(); Ep = parse_delE()
        cngele(phase - diff_step, sect); run_parmela(); Em = parse_delE()
        grad = (Ep - Em) / (2 * diff_step)
        vel = momentum * vel - lr * grad
        phase += vel
        cngele(phase, sect); run_parmela(); E = parse_delE()
        print(f"Gradient it {it}: phase={phase}, ΔE={E}, grad={grad}, vel={vel}, step={diff_step}")
        # adjust step size
        if prevE is not None:
            diff_step *= step_up if E < prevE else step_down
            if abs(prevE - E) < tol:
                print(f"Converged at iter {it}, ΔE change < {tol}")
                break
        prevE = E
    return phase, E

# Parabola Fitting Optimizer w/ convergence check
def optimize_p(dp, sect, its):
    idxs, _ = find_indices()
    orig = float(open(default_temp).read().splitlines()[idxs[sect]].split()[4]) - 90.0
    # initial points
    pts = [orig, orig + dp, orig - dp]
    res = []
    prevEv = None
    for p in pts:
        cngele(p, sect); run_parmela(); Ev = parse_delE()
        res.append((p, Ev))
        print(f"Initial point: phase={p}, ΔE={Ev}")
        if prevEv is not None and abs(prevEv - Ev) < tol:
            print(f"Converged at initial evaluation, ΔE change < {tol}")
            return p, Ev
        prevEv = Ev
    # iterative fitting
    for i in range(1, its + 1):
        xs = np.array([x for x, _ in res]); ys = np.array([y for _, y in res])
        a, b, _ = np.polyfit(xs, ys, 2)
        vert = -b / (2 * a)
        cngele(vert, sect); run_parmela(); Ev = parse_delE()
        res.append((vert, Ev))
        print(f"Parabola it {i}: phase={vert}, ΔE={Ev}")
        if abs(prevEv - Ev) < tol:
            print(f"Converged at parabola it {i}, ΔE change < {tol}")
            return vert, Ev
        prevEv = Ev
        res = sorted(res, key=lambda t: t[1])[:3]
    return min(res, key=lambda t: t[1])

# Dispatch and final run
if mode == 'g':
    init_step = float(sys.argv[3]); sect = int(sys.argv[4]); mit = int(sys.argv[5]); lr = float(sys.argv[6])
    best_phase, best_de = optimize_g(init_step, sect, mit, lr)
elif mode == 'p':
    dp = float(sys.argv[3]); sect = int(sys.argv[4]); its = int(sys.argv[5])
    best_phase, best_de = optimize_p(dp, sect, its)
else:
    print_usage_and_exit()

# Apply optimal phase one more time and run
cngele(best_phase, sect)
run_parmela()

print(f"Optimal result: phase={best_phase}, ΔE={best_de}")
sys.exit(0)
