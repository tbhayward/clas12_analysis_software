import re

# Function to extract kinematic variables and z-bin from LaTeX content
def extract_kinematic_variables_and_zbin(latex_content):
    # Define a regular expression pattern to extract z-bin from table caption
    caption_pattern = re.compile(r"\\caption{The mean kinematic variables in each of the bins for the extracted \$Q2y(\d)z(\d)\$ asymmetries. Values given in GeV or GeV\$\\^2\$ where appropriate.}")
    
    # Find the first match of z-bin in the caption
    caption_match = caption_pattern.search(latex_content)
    if caption_match:
        package_number, z_bin = caption_match.groups()
    else:
        print("No caption match found. Check the LaTeX file format.")
        return None, None
    
    # Define a regular expression pattern to match rows in the LaTeX tables
    pattern = re.compile(r"\\hline\s+(\d+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([-\d.]+)~&~([-\d.]+)~&~([-\d.]+)\\\\\s+\\hline")
    
    # Find all matches of the pattern in the LaTeX content
    matches = pattern.findall(latex_content)

    # Initialize kinematic variables with default values of -100
    kinematic_variables = {
        'Q2': [-100] * 5,
        'W': [-100] * 5,
        'xB': [-100] * 5,
        'y': [-100] * 5,
        'z': [-100] * 5,
        'zeta': [-100] * 5,
        'PT': [-100] * 5,
        'xF': [-100] * 5,
        't': [-100] * 5,
    }
    
    # Fill in the actual values from matches
    for match in matches:
        index = int(match[0]) - 1
        for i, key in enumerate(kinematic_variables.keys(), 1):
            kinematic_variables[key][index] = float(match[i])

    return kinematic_variables, package_number, z_bin

# Function to write variables to a file
def write_variables_to_file(kinematic_variables, package_number, z_bin, output_file):
    with open(output_file, 'a') as f:
        for key, values in kinematic_variables.items():
            variable_name = f"Q2y{package_number}z{z_bin}_{key}"
            values_str = ', '.join(map(str, values))
            f.write(f"{variable_name} = {{{values_str}}};\n")

# Main function to process the LaTeX file
def process_latex_file(latex_file_path, output_file):
    with open(latex_file_path, 'r') as file:
        latex_content = file.read()

    kinematic_variables, package_number, z_bin = extract_kinematic_variables_and_zbin(latex_content)
    if kinematic_variables:
        write_variables_to_file(kinematic_variables, package_number, z_bin, output_file)
        print(f"Processed and wrote variables for package {package_number}, z-bin {z_bin}.")

# Example usage
latex_file_path = "/volatile/clas12/thayward/UU_validation/packages/1/clas12_analysis_software/analysis_scripts/asymmetry_extraction/output/results/kinematics_rga_fa18_inb_epX_skimmed_timeStamp_03_06_104135.txt"  # Update this path
output_file = "output2/analysis_statistics.txt"
process_latex_file(latex_file_path, output_file)
