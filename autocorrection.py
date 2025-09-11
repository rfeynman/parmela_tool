#!/usr/bin/env python3
import sys
import subprocess
import numpy as np
import os
import shutil

# Usage:
#   python autocorrection.py <input_file.inp> <delta_val> <start_sect> <iterations>
#
# Example:
#   python3 autocorrection.py rr6_with_cors.inp 0.1 0 10
#
# This script sequentially optimizes all steerer settings in a file to minimize beam orbit.

# --- Configuration ---
# Convergence tolerance: optimization stops if orbit change is less than this value.
TOLERANCE = 2e-4

# --- Global State ---
# These will be initialized once in the main execution block.
steerer_indices = []
cor_indices = []

def print_usage_and_exit():
    """Prints the script usage information and exits."""
    print("Usage: python autocorrection.py <input_file.inp> <delta_val> <start_sect> <iterations>")
    sys.exit(1)

# --- Core Functions ---

def find_indices(filename):
    """
    Finds the line numbers of all 'steerer' and '!cor' entries in a file.
    
    Args:
        filename (str): The path to the input file.
        
    Returns:
        tuple: (list of steerer indices, list of !cor marker indices).
    """
    s_indices = []
    c_indices = []
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            stripped_line = line.strip().lower()
            if stripped_line.startswith("steerer"):
                s_indices.append(i)
            elif stripped_line.startswith("!cor"):
                c_indices.append(i)
    except FileNotFoundError:
        print(f"Error: Input file {filename} not found.")
        sys.exit(1)
    return s_indices, c_indices

def run_parmela():
    """Runs the Parmela simulation."""
    try:
        # Using DEVNULL to hide Parmela's stdout for a cleaner output
        subprocess.run(["parmela", default_temp], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Error running Parmela: {e}")
        print("Please ensure 'parmela' is in your system's PATH.")
        sys.exit(1)

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

def modify_steerer(xvalue, yvalue, sect, filename, truncate=False):
    """
    Changes a steerer's values and optionally truncates the file for simulation.

    Args:
        xvalue (float): The new X value for the steerer.
        yvalue (float): The new Y value for the steerer.
        sect (int): The index of the steerer to modify.
        filename (str): The file to modify.
        truncate (bool): If True, adds an 'end' command and comments out the previous section.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # If truncating, modify the file for a temporary simulation run.
    if truncate:
        # Remove any 'end' command from a previous run
        lines = [line for line in lines if line.strip().lower() != 'end']

        # If this is not the first section, comment out the previous section
        # to isolate the effect of the current steerer.
        if sect > 0:
            start_comment_idx = cor_indices[sect - 1]
            end_comment_idx = cor_indices[sect]
            
            for i in range(start_comment_idx, end_comment_idx):
                # BUG FIX: Do not comment out the previous steerer line. Its correction must remain active.
                if i == steerer_indices[sect - 1]:
                    continue
                # Only add '!' if the line is not empty and not already a comment.
                if lines[i].strip() and not lines[i].lstrip().startswith('!'):
                    lines[i] = '!' + lines[i]

    # Update the specified steerer's X and Y values
    if sect < len(steerer_indices):
        line_index = steerer_indices[sect]
        parts = lines[line_index].split()
        parts[4] = str(xvalue)
        parts[5] = str(yvalue)
        lines[line_index] = " ".join(parts) + "\n"
    else:
        print(f"Error: Steerer section {sect} is out of bounds.")
        return

    # If truncating, add a new 'end' command to stop the simulation after this section
    if truncate:
        if sect < len(cor_indices):
            cor_line_index = cor_indices[sect]
            end_pos = min(cor_line_index + 4, len(lines))
            lines.insert(end_pos, "end\n")
        else:
            lines.append("end\n")
            print(f"Warning: No '!cor' marker for section {sect}. Appending 'end' to file.")

    with open(filename, 'w') as f:
        f.writelines(lines)

# --- Optimization Logic ---

def _fit_linear(initial_val, fixed_val, sect, its, dp, optimize_x):
    """
    Helper function to find the optimal value for one variable (X or Y) using a linear fit.
    """
    opt_var_name = "x" if optimize_x else "y"
    orbit_name = "x-orbit" if optimize_x else "y-orbit"
    
    # Create two initial points to define the first line.
    # Use a relative step (multiplicative) if initial_val is not zero to avoid large jumps.
    # Fall back to an absolute step (additive) if initial_val is zero.
    if initial_val != 0:
        second_pt = initial_val * (1 + dp)
    else:
        second_pt = dp
    
    pts = [initial_val, second_pt]
    results = []
    
    # Evaluate the initial points
    for p in pts:
        if optimize_x:
            modify_steerer(p, fixed_val, sect, default_temp, truncate=True)
        else:
            modify_steerer(fixed_val, p, sect, default_temp, truncate=True)
        run_parmela()
        x_orbit, y_orbit = parse_orbits()
        current_orbit = x_orbit if optimize_x else y_orbit
        if current_orbit is None: continue
        
        results.append((p, current_orbit))
        print(f"  Initial point: {opt_var_name}={p:.6f}, {orbit_name}={current_orbit:.6f}")

    # Iteratively find the zero-crossing using a linear fit
    for i in range(1, its + 1):
        results.sort(key=lambda t: abs(t[1]))
        results = results[:2]

        if abs(results[0][1]) < TOLERANCE:
            print(f"  Converged: Orbit magnitude is less than {TOLERANCE}.")
            break
            
        if len(results) < 2:
            print("  Warning: Not enough valid points to fit a line.")
            break

        x_vals = np.array([r[0] for r in results])
        y_vals_signed = np.array([r[1] for r in results])
        
        m, c_linear = np.polyfit(x_vals, y_vals_signed, 1)

        if m == 0:
            print("  Warning: Linear fit slope is zero. Cannot find zero-crossing.")
            break
        
        next_val = -c_linear / m
        
        if optimize_x:
            modify_steerer(next_val, fixed_val, sect, default_temp, truncate=True)
        else:
            modify_steerer(fixed_val, next_val, sect, default_temp, truncate=True)
        run_parmela()
        x_orbit, y_orbit = parse_orbits()
        current_orbit = x_orbit if optimize_x else y_orbit
        if current_orbit is None: continue

        results.append((next_val, current_orbit))
        print(f"  Linear fit it {i}: {opt_var_name}={next_val:.6f}, {orbit_name}={current_orbit:.6f}")

    best_val, best_orbit = min(results, key=lambda t: abs(t[1]))
    return best_val, best_orbit

def optimize_p(dp, sect, its):
    """
    Main optimization function for a single section.
    """
    # Get original X and Y values from the current state of the temp file
    with open(default_temp, 'r') as f:
        line = f.read().splitlines()[steerer_indices[sect]]
    
    try:
        parts = line.split()
        orig_x = float(parts[4])
        orig_y = float(parts[5])
    except (ValueError, IndexError):
        print("\nError: Could not read initial X and Y values for the steerer.")
        sys.exit(1)

    print(f"Starting optimization for steerer {sect} from X={orig_x}, Y={orig_y}")

    # --- Step 1: Optimize X value ---
    print("\n--- Optimizing X value ---")
    best_x, best_x_orbit = _fit_linear(orig_x, orig_y, sect, its, dp, optimize_x=True)
    print(f"--> Best X value found: {best_x:.6f} (x-orbit: {best_x_orbit:.6f})")

    # --- Step 2: Optimize Y value ---
    print("\n--- Optimizing Y value ---")
    best_y, best_y_orbit = _fit_linear(orig_y, best_x, sect, its, dp, optimize_x=False)
    print(f"--> Best Y value found: {best_y:.6f} (y-orbit: {best_y_orbit:.6f})")

    # --- Final Check ---
    modify_steerer(best_x, best_y, sect, default_temp, truncate=True)
    run_parmela()
    final_x_orbit, final_y_orbit = parse_orbits()
    print(f"Final orbits for section {sect}: x-orbit={final_x_orbit:.6f}, y-orbit={final_y_orbit:.6f}")

    return best_x, best_y, final_x_orbit, final_y_orbit

# --- Main Execution ---
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print_usage_and_exit()

    base_inp = sys.argv[1]
    try:
        delta_val = float(sys.argv[2])
        start_section = int(sys.argv[3])
        iterations = int(sys.argv[4])
    except ValueError:
        print("Error: <delta_val>, <start_sect>, and <iterations> must be numbers.")
        print_usage_and_exit()

    root, ext = os.path.splitext(base_inp)
    if ext.lower() != ".inp":
        print(f"Error: Input file '{base_inp}' must have a .inp extension.")
        sys.exit(1)
        
    default_temp = f"{root}_temp{ext}"
    default_tbl = "TIMESTEPEMITTANCE.TBL"

    # Create a fresh temporary copy of the input file to work with.
    try:
        shutil.copy(base_inp, default_temp)
    except (FileNotFoundError, IOError) as e:
        print(f"Error: Could not copy '{base_inp}' to '{default_temp}': {e}")
        sys.exit(1)

    # Find all steerers and markers once from the clean file.
    steerer_indices, cor_indices = find_indices(default_temp)
    #print(f"steerer num='{steerer_indices}'; cor num='{cor_indices}'")
    if not steerer_indices:
        print("Error: No 'steerer' elements found in the input file.")
        sys.exit(1)
    
    if not cor_indices:
        print("Error: No '!cor' markers found in the input file.")
        print("Please add markers like '!cor 0', '!cor 1', etc., to your file after each section.")
        sys.exit(1)

    num_sections = len(cor_indices)

    # --- Main Automation Loop ---
    for sect in range(start_section, num_sections):
        print("\n" + "="*50)
        print(f"      STARTING OPTIMIZATION FOR SECTION {sect}")
        print("="*50)

        # Run optimization for the current section. This will leave default_temp
        # in a temporary state (with comments and an 'end' command).
        bestx, besty, _, _ = optimize_p(delta_val, sect, iterations)
        
        if bestx is None:
            print(f"Optimization failed for section {sect}. Stopping.")
            break
        
        # Finalize the state of default_temp for this section.
        # 1. Read the current temporary file state.
        with open(default_temp, 'r') as f:
            lines = f.readlines()
            
        # 2. Remove the temporary 'end' command added by the optimizer.
        lines = [line for line in lines if line.strip().lower() != 'end']

        # 3. Permanently set the new best steerer values in the lines.
        #    The comments added by modify_steerer() will remain.
        line_index = steerer_indices[sect]
        parts = lines[line_index].split()
        parts[4] = str(bestx)
        parts[5] = str(besty)
        lines[line_index] = " ".join(parts) + "\n"
        
        # 4. Write the finalized, corrected lines back to the temp file.
        with open(default_temp, 'w') as f:
            f.writelines(lines)

        print(f"\nApplied final correction for section {sect}: X={bestx:.6f}, Y={besty:.6f}")

    print("\n" + "="*50)
    print("      ALL SECTIONS OPTIMIZED")
    print("="*50)
    
    # For the final run, re-enable all simulation sections by removing the '!'
    # comments that were added during the optimization steps.
    print("Re-enabling all sections for final simulation...")
    with open(default_temp, 'r') as f:
        lines = f.readlines()

    uncommented_lines = []
    # These are the keywords for lines that get commented out by the optimizer.
    keywords_to_uncomment = ('!scheff', '!restart', '!save','!start')
    for line in lines:
        stripped_line = line.lstrip()
        # Check if the line is one of the commented-out types.
        if stripped_line.lower().startswith(keywords_to_uncomment):
            # Find the position of the '!' and remove it, preserving original indentation.
            first_char_index = len(line) - len(stripped_line)
            uncommented_lines.append(line[:first_char_index] + line[first_char_index+1:])
        else:
            uncommented_lines.append(line)
            
    with open(default_temp, 'w') as f:
        f.writelines(uncommented_lines)

    # Run a final simulation with all corrections applied to get the final orbit.
    print("Running final simulation with all corrections...")
    run_parmela()
    final_x_orbit, final_y_orbit = parse_orbits()

    print("\n--- FINAL RESULTS ---")
    print(f"  Final X orbit:   {final_x_orbit:.6f}")
    print(f"  Final Y orbit:   {final_y_orbit:.6f}")
    print(f"  Final corrected input file is in: {default_temp}")
    print("="*50)
    
    sys.exit(0)


