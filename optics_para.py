import pandas as pd
import os

# Load the Excel file
filename = "rr7.inp"
file_path = "/Users/wange/Documents/Research/eRHIC injector/eRHIC baseline/Beamline/Lattice/sband2025/"  # Replace with the path to your file
file_path_name = os.path.join(file_path, filename)
ele_name_list = ["Solenoid","cell","quad","trwave"]

def get_ele_value(ele_name_list, file_path_name):
    # Initialize an empty list to store the lines
    lines = []
    # Read the file and filter lines
    with open(file_path_name, 'r') as file:
        for line in file:
            for name in ele_name_list:
                if line.startswith(name):
                    parts = line.split()
                    if len(parts) >= 4:
                        if line.startswith("Solenoid"):
                            filtered_line = f"{parts[0]} {parts[1]} {parts[2]} {parts[4]}"
                            lines.append(filtered_line)
                        elif line.startswith("cell"):
                            filtered_line = f"{parts[0]} {parts[1]}  {parts[4]} {parts[5]}"
                            lines.append(filtered_line)
                        elif line.startswith("quad"):
                            filtered_line = f"{parts[0]} {parts[1]} {parts[2]} {parts[4]}"
                            lines.append(filtered_line)
                        elif line.startswith("trwave"):
                            filtered_line = f"{parts[0]} {parts[1]} {parts[4]} {parts[5]}"
                            lines.append(filtered_line)

    # Create a DataFrame from the filtered lines
    df = pd.DataFrame(lines, columns=["Line"])

    # Display the DataFrame
    print(df)
    #output_file_path = os.path.join(file_path, "df.csv")
    #df.to_csv(output_file_path, index=False, header=False)
    return df

def ele_process(df_org):
    # Split the lines into separate columns
    df_split = df_org['Line'].str.split(expand=True)
    df_split.columns = ['Element', 'Length', 'Apt_or_phase','Amp']

    # Convert the Length column to numeric
    df_split['Length'] = pd.to_numeric(df_split['Length'])

    # Process only 'trwave' lines to sum the Length for identical lines
    df_trwave = df_split[df_split['Element'] == 'trwave']
    df_other = df_split[df_split['Element'] != 'trwave']

    df_trwave_grouped = df_trwave.groupby(['Element', 'Apt_or_phase', 'Amp']).agg({'Length': 'sum'}).reset_index()

    # Combine the columns back into a single string
    #df_trwave_grouped['Line'] = df_trwave_grouped['Element'] + ' ' + df_trwave_grouped['Length'].astype(str) + ' '+df_trwave_grouped['Apt_or_phase']+' ' + df_trwave_grouped['Amp']
    #df_other['Line'] = df_other['Element'] + ' ' + df_other['Length'].astype(str) + ' '+df_other['Apt_or_phase']+' ' + df_other['Amp']

    # Concatenate the processed trwave lines with the other lines
    df_result = pd.concat([df_other, df_trwave_grouped])

    # Display the DataFrame
    print(df_result)

    return df_result
     
def main():
    df_org=get_ele_value(ele_name_list, file_path_name) 
    req=ele_process(df_org)
    output_file_path = os.path.join(file_path, "requirement.xlsx")
    req.to_excel(output_file_path, index=False)

if __name__ == "__main__":
    main()
    