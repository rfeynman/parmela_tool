import sys
import pandas as pd
import matplotlib.pyplot as plt

def plot_from_file(file_path):
    """
    Reads orbit data from a text file and generates plots.
    """
    try:
        # Extract run_id from the filename
        run_id = file_path.split('orbit_error_')[-1].split('.txt')[0]
    except IndexError:
        run_id = "unknown"

    try:
        data = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        return

    z_col = 'Z(cm)_1'
    if z_col not in data.columns:
        print(f"Error: Base Z column '{z_col}' not found in {file_path}.")
        return

    # Find all X and Y columns dynamically
    x_cols = [col for col in data.columns if col.startswith('<X>(mm)_')]
    y_cols = [col for col in data.columns if col.startswith('<Y>(mm)_')]

    # Plot X orbits
    plt.figure(figsize=(10, 6))
    for x_col in x_cols:
        plt.plot(data[z_col], data[x_col])
    plt.xlabel('Z (cm)')
    plt.ylabel('<X> (mm)')
    plt.title(f'Orbit Error Analysis (X) - Run ID: {run_id}')
    plt.grid(True)
    plt.savefig(f'orbit_X_{run_id}.png', dpi=300)
    print(f"Saved high-resolution plot to orbit_X_{run_id}.png")
    plt.show()

    # Plot Y orbits
    plt.figure(figsize=(10, 6))
    for y_col in y_cols:
        plt.plot(data[z_col], data[y_col])
    plt.xlabel('Z (cm)')
    plt.ylabel('<Y> (mm)')
    plt.title(f'Orbit Error Analysis (Y) - Run ID: {run_id}')
    plt.grid(True)
    plt.savefig(f'orbit_Y_{run_id}.png', dpi=300)
    print(f"Saved high-resolution plot to orbit_Y_{run_id}.png")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python orbitplot.py <path_to_orbit_error_file.txt>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    plot_from_file(input_file)
