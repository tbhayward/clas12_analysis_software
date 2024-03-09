import os
import re

# Define base directories and output file path
base_dir = "/volatile/clas12/thayward/UU_validation/packages"
output_file_path = "output/analysis_statistics.txt"

# Ensure the output directory exists
os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

# Function to insert default values if a key is missing
def insert_default(data, keys, default):
    for key in keys:
        if key not in data:
            data[key] = [default] * 5

# Function to parse asymmetry files and update the data dictionary
def parse_asymmetries(file_path, q2y_bin, data):
    with open(file_path) as f:
        content = f.read()
        
        # Regular expression to match the data format
        pattern = re.compile(r"(\w+) = \{\{([^\}]+)\}\};")
        
        for match in pattern.finditer(content):
            variable_name, values = match.groups()
            processed_values = [list(map(float, val.split(', '))) for val in values.split('}, {')]
            data[variable_name] = processed_values

# Function to parse LaTeX tables and update the data dictionary
def parse_latex(file_path, data):
    with open(file_path) as f:
        content = f.read()
        
        # Regular expression to match LaTeX tables for kinematics
        table_pattern = re.compile(r"\\begin\{table\}.*?\\end\{table\}", re.DOTALL)
        value_pattern = re.compile(r"(\d+(?:\.\d+)?~&~)+")
        
        for table in table_pattern.findall(content):
            # Attempt to extract the Q2y and z bin from the caption
            caption_match = re.search(r"\\caption\{.*?Q2y(\d+)z(\d+).*?\}", table)
            if not caption_match:
                continue  # Skip if we can't find a caption match
            
            q2y_bin, z_bin = caption_match.groups()
            
            # Pre-generate keys for missing values insertion
            keys = [f"Q2y{q2y_bin}z{z_bin}_{var}" for var in ["Q2", "W", "xB", "y", "z", "zeta", "PT", "xF", "t"]]
            insert_default(data, keys, -100)
            
            for row_match in value_pattern.finditer(table):
                values = row_match.group().split("~&~")
                for i, val in enumerate(values[:-1]):  # Exclude the last empty split
                    key = keys[i]
                    data[key].append(float(val))

# Initialize the data dictionary
data = {}

# Process each directory (Q2y bin)
for q2y_bin in range(1, 18):
    dir_path = os.path.join(base_dir, str(q2y_bin))
    
    # Default insertion for missing asymmetry data
    for z_bin in range(1, 6):
        keys = [
            f"Q2y{q2y_bin}z{z_bin}chi2FitsALUoffset",
            f"Q2y{q2y_bin}z{z_bin}chi2FitsALUsinphi",
            f"Q2y{q2y_bin}z{z_bin}MLMFitsALUsinphi",
            f"Q2y{q2y_bin}z{z_bin}MLMFitsAULsinphi",
            f"Q2y{q2y_bin}z{z_bin}MLMFitsAULsin2phi",
            f"Q2y{q2y_bin}z{z_bin}MLMFitsALL",
            f"Q2y{q2y_bin}z{z_bin}MLMFitsALLcosphi",
            f"Q2y{q2y_bin}z{z_bin}MLMFitsAUUcosphi",
            f"Q2y{q2y_bin}z{z_bin}MLMFitsAUUcos2phi",
            # Add more keys as per your requirement
        ]
        insert_default(data, keys, [-1, 0, 0])
    
    # Process files in the directory
    for file_name in os.listdir(dir_path):
        file_path = os.path.join(dir_path, file_name)
        
        if file_name.startswith("asymmetries"):
            parse_asymmetries(file_path, q2y_bin, data)
        elif file_name.startswith("kinematics"):
            parse_latex(file_path, data)

# Write the processed data to the output file
with open(output_file_path, "w") as f_out:
    for key, value_lists in data.items():
        # Ensure each sublist is formatted as a Mathematica list; handle edge cases where lists might be empty or incomplete
        formatted_values = ", ".join(
            "{" + ", ".join(map(lambda x: str(x) if x != -1 else "{-1, 0, 0}", values)) + "}"
            for values in value_lists
        )
        # Write the formatted string to the output file, ensuring Mathematica list syntax is used
        f_out.write(f"{key} = {{{formatted_values}}};\n")


