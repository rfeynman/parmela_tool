#!/usr/bin/env python3
import sys
import subprocess
import numpy as np
import os

# Usage:
#   python autocorrection.py <input_file.inp> <delta_val> <sect> <iterations>
#
# Example:
#   python3 autocorrection.py rr6.inp 30 0 10
#
# This script optimizes steerer settings to minimize beam orbit deviation.
# Two models are used. Parabola fitting and linear. If the corrector BL is larger than certain number, then choose linear fitting. Linear fitting is more effective.
# --- Configuration ---
# Convergence tolerance: optimization stops if orbit change is less than this value.
TOLERANCE = 2e-4

def print_usage_and_exit():
    """Prints the script usage information and exits."""
    print("Usage: python autocorrection.py <input_file.inp> <delta_val> <sect> <iterations>")
    sys.exit(1)

# --- Core Functions ---

def find_indices():
    """
    Finds the line numbers of all 'steerer' entries in the input file.
    
    Returns:
        list: A list of integer line numbers (0-based).
    """
    indices = []
    try:
        with open(default_temp, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.strip().lower().startswith("steerer"):
                indices.append(i)
    except FileNotFoundError:
        print(f"Error: Temporary file {default_temp} not found.")
        sys.exit(1)
    return indices

def cngele(xvalue, yvalue, sect=0):
    """
    Changes a steerer's X and Y values in the temporary input file.

    Args:
        xvalue (float): The new X value for the steerer.
        yvalue (float): The new Y value for the steerer.
        sect (int): The index of the steerer to modify.
    """
    indices = find_indices()
    if not indices or sect >= len(indices):
        print(f"Error: Steerer section {sect} not found in the file.")
        return

    line_index = indices[sect]
    with open(default_temp, 'r') as f:
        lines = f.read().splitlines()
    
    parts = lines[line_index].split()
    # Update the X (index 4) and Y (index 5) values
    parts[4] = str(xvalue)
    parts[5] = str(yvalue)
    lines[line_index] = " ".join(parts)
    
    with open(default_temp, 'w') as f:
        f.write("\n".join(lines) + "\n")

def parse_orbits():
    """
    Parses the output file to get the final X and Y orbits.

    Returns:
        tuple: A tuple containing the x_orbit (float) and y_orbit (float).
               Returns (None, None) if parsing fails.
    """
    try:
        with open(default_tbl, 'r') as f:
            lines = f.read().splitlines()
        # The orbit data is in the second to last line of the table
        cols = lines[-2].split()
        x_orbit = float(cols[13])
        y_orbit = float(cols[15])
        return x_orbit, y_orbit
    except (FileNotFoundError, IndexError, ValueError) as e:
        print(f"Error parsing orbit file {default_tbl}: {e}")
        return None, None

def run_parmela():
    """Runs the Parmela simulation."""
    try:
        # Using DEVNULL to hide Parmela's stdout for a cleaner output
        subprocess.run(["parmela", default_temp], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Error running Parmela: {e}")
        print("Please ensure 'parmela' is in your system's PATH.")
        sys.exit(1)

# --- Optimization Logic ---

def _fit_parabola(initial_val, fixed_val, sect, its, dp, optimize_x):
    """
    Helper function to run parabola optimization on one variable (X or Y).
    It finds the value for the optimized variable that minimizes the corresponding orbit.

    Args:
        initial_val (float): The starting value for the variable to be optimized.
        fixed_val (float): The fixed value for the other variable.
        sect (int): The steerer index.
        its (int): The number of fitting iterations.
        dp (float): The delta step to create the initial points for the parabola.
        optimize_x (bool): True to optimize X orbit, False to optimize Y orbit.

    Returns:
        tuple: The best value found and its corresponding orbit.
    """
    opt_var_name = "x" if optimize_x else "y"
    orbit_name = "x-orbit" if optimize_x else "y-orbit"
    
    # Create three initial points to define the first parabola
    pts = [initial_val, initial_val + dp, initial_val - dp]
    results = []
    
    # Evaluate the initial points
    for p in pts:
        if optimize_x:
            cngele(p, fixed_val, sect)
        else:
            cngele(fixed_val, p, sect)
        run_parmela()
        x_orbit, y_orbit = parse_orbits()
        current_orbit = x_orbit if optimize_x else y_orbit
        if current_orbit is None: continue # Skip if parsing failed
        
        results.append((p, current_orbit))
        print(f"  Initial point: {opt_var_name}={p:.6f}, {orbit_name}={current_orbit:.6f}")

    # Iteratively find the minimum of the parabola
    for i in range(1, its + 1):
        # Sort points by the absolute value of their orbit and keep the best three
        results.sort(key=lambda t: abs(t[1]))
        results = results[:3]

        # Check for convergence
        if len(results) > 1 and abs(results[0][1] - results[1][1]) < TOLERANCE:
            print(f"  Converged: Orbit change is less than {TOLERANCE}.")
            break
            
        if len(results) < 3:
            print("  Warning: Not enough valid points to fit a parabola. Stopping this fit.")
            break

        # Fit a 2nd-degree polynomial (a*x^2 + b*x + c) to the points
        # We fit against the absolute orbit value to find the minimum magnitude
        x_vals = np.array([r[0] for r in results])
        y_vals = np.array([abs(r[1]) for r in results])
        
        a, b, _ = np.polyfit(x_vals, y_vals, 2)

        # The vertex of the parabola y = ax^2 + bx + c is at x = -b / (2a)
        # This vertex is our prediction for the optimal value.
        if a == 0:
            # If 'a' is not positive, the parabola opens downwards, so there's no minimum.
            # We'll stop and use the best point we've found so far.
            print("  Warning: Parabolic fit resulted in a line. Using best point found so far.")
            break
        
        vertex = -b / (2 * a)

        if np.abs(vertex) > 20:
            print(f"  Warning: Parabolic vertex is large ({vertex:.2f}). Switching to linear fit.")
            y_vals_signed = np.array([r[1] for r in results])
            m, c_linear = np.polyfit(x_vals, y_vals_signed, 1) # Linear fit: y = mx + c
            if m != 0:
                vertex = -c_linear / m # Find x where y=0
                print(f"  Linear fit suggests next point: {vertex:.6f}")
            else:
                print("  Warning: Linear fit slope is zero. Cannot find zero-crossing. Stopping.")
                break

        # Run simulation with the new predicted optimal value
        if optimize_x:
            cngele(vertex, fixed_val, sect)
        else:
            cngele(fixed_val, vertex, sect)
        run_parmela()
        x_orbit, y_orbit = parse_orbits()
        current_orbit = x_orbit if optimize_x else y_orbit
        if current_orbit is None: continue

        results.append((vertex, current_orbit))
        print(f"  Parabola it {i}: {opt_var_name}={vertex:.6f}, {orbit_name}={current_orbit:.6f}")

    # After all iterations, return the best value found
    best_val, best_orbit = min(results, key=lambda t: abs(t[1]))
    return best_val, best_orbit

def optimize_p(dp, sect, its):
    """
    Main optimization function. First optimizes the X value, then the Y value.

    Args:
        dp (float): The delta step for the parabola fitting.
        sect (int): The steerer index to optimize.
        its (int): The number of iterations for the fitting process.

    Returns:
        tuple: Contains best_x, best_y, final_x_orbit, final_y_orbit.
    """
    indices = find_indices()
    if not indices or sect >= len(indices):
        print(f"Error: Cannot start optimization, steerer {sect} not found.")
        return None, None, None, None
        
    line_index = indices[sect]
    
    # Get original X and Y values from the file
    with open(default_temp, 'r') as f:
        line = f.read().splitlines()[line_index]
    
    try:
        parts = line.split()
        orig_x = float(parts[4])
        orig_y = float(parts[5])
    except (ValueError, IndexError):
        print("\nError: Could not read initial X and Y values for the steerer.")
        print("Please check the format of the following line in your input file:")
        print(f" -> \"{line.strip()}\"")
        print("This line should have at least 6 columns, with numeric values in columns 5 and 6.")
        sys.exit(1)

    print(f"Starting optimization for steerer {sect} from X={orig_x}, Y={orig_y}")

    # --- Step 1: Optimize X value, keeping Y fixed at its original value ---
    print("\n--- Optimizing X value ---")
    best_x, best_x_orbit = _fit_parabola(orig_x, orig_y, sect, its, dp, optimize_x=True)
    print(f"--> Best X value found: {best_x:.6f} (x-orbit: {best_x_orbit:.6f})")

    # --- Step 2: Optimize Y value, keeping X fixed at its new best value ---
    print("\n--- Optimizing Y value ---")
    best_y, best_y_orbit = _fit_parabola(orig_y, best_x, sect, its, dp, optimize_x=False)
    print(f"--> Best Y value found: {best_y:.6f} (y-orbit: {best_y_orbit:.6f})")

    # --- Final Check: Run with both best values to get the final combined orbits ---
    print("\n--- Final Check ---")
    cngele(best_x, best_y, sect)
    run_parmela()
    final_x_orbit, final_y_orbit = parse_orbits()
    print(f"Final orbits with optimal settings: x-orbit={final_x_orbit:.6f}, y-orbit={final_y_orbit:.6f}")

    return best_x, best_y, final_x_orbit, final_y_orbit

# --- Main Execution ---
if __name__ == "__main__":
    # Validate command-line arguments
    if len(sys.argv) != 5:
        print_usage_and_exit()

    # Parse arguments
    base_inp = sys.argv[1]
    try:
        delta_val = float(sys.argv[2])
        section = int(sys.argv[3])
        iterations = int(sys.argv[4])
    except ValueError:
        print("Error: <delta_val>, <sect>, and <iterations> must be numbers.")
        print_usage_and_exit()

    root, ext = os.path.splitext(base_inp)
    if ext.lower() != ".inp":
        print(f"Error: Input file '{base_inp}' must have a .inp extension.")
        print_usage_and_exit()
        
    # Define file paths
    default_temp = f"{root}_temp{ext}"
    default_tbl = "TIMESTEPEMITTANCE.TBL"

    # Create a fresh temporary copy of the input file for modification
    try:
        subprocess.run(["cp", base_inp, default_temp], check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"Error: Could not copy '{base_inp}' to '{default_temp}'.")
        sys.exit(1)

    # Run the optimization
    bestxvalue, bestyvalue, bestxorbit, bestyorbit = optimize_p(delta_val, section, iterations)

    if bestxvalue is None:
        print("\nOptimization failed.")
        sys.exit(1)

    # Apply optimal values one last time to ensure the final state is set
    print("\nApplying final optimal values...")
    cngele(bestxvalue, bestyvalue, section)
    run_parmela()

    # Print the final results
    print("\n" + "="*35)
    print("      OPTIMIZATION COMPLETE")
    print("="*35)
    print(f"  Optimal X value: {bestxvalue:.6f}")
    print(f"  Optimal Y value: {bestyvalue:.6f}")
    print(f"  Final X orbit:   {bestxorbit:.6f}")
    print(f"  Final Y orbit:   {bestyorbit:.6f}")
    print("="*35)
    
    sys.exit(0)
