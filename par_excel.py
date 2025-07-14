import pandas as pd
import os

# Load the Excel file
file_path = "/Users/wange/Documents/Research/eRHIC injector/eRHIC baseline/Beamline/Lattice/sband2025/py_pro.xlsx"  # Replace with the path to your file
df = pd.read_excel(file_path, sheet_name="Sheet1")

# Initialize a list to store processed rows
processed_rows = []
current_element = None
start_index = None

# Iterate through the dataframe to group and process rows
for i, row in df.iterrows():
    if row["element"] != current_element:  # If a new element group starts
        # If there's an active group, finalize it
        if current_element is not None and start_index is not None:
            df.loc[start_index, "z2"] = df.loc[i - 1, "z2"]
            processed_rows.append(df.iloc[start_index])

        # Start a new group
        current_element = row["element"]
        start_index = i

# Process the final group
if current_element is not None and start_index is not None:
    df.loc[start_index, "z2"] = df.loc[len(df) - 1, "z2"]
    processed_rows.append(df.iloc[start_index])

# Convert the list of processed rows into a new dataframe
result_df = pd.DataFrame(processed_rows)
directory=os.path.dirname(file_path)
output_file=os.path.join(directory, "processed_file.xlsx")

# Save the processed dataframe to a new Excel file
result_df.to_excel(output_file, index=False)

print("Processing complete. File saved as 'processed_file.xlsx'")
