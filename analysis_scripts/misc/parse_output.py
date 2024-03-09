# Import regular expressions module
import re

# Function to extract kinematic variables from LaTeX content
def extract_kinematic_variables(latex_content):
    # Define a regular expression pattern to match rows in the LaTeX tables
    pattern = re.compile(r"\\hline\s+(\d+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([-\d.]+)~&~([-\d.]+)~&~([-\d.]+)\\\\\s+\\hline")
    
    # Find all matches of the pattern in the LaTeX content
    matches = pattern.findall(latex_content)
    
    # Initialize a dictionary to store kinematic variables with default values of -100
    kinematic_variables = {
        'Q2': [-100, -100, -100, -100, -100],
        'W': [-100, -100, -100, -100, -100],
        'xB': [-100, -100, -100, -100, -100],
        'y': [-100, -100, -100, -100, -100],
        'z': [-100, -100, -100, -100, -100],
        'zeta': [-100, -100, -100, -100, -100],
        'PT': [-100, -100, -100, -100, -100],
        'xF': [-100, -100, -100, -100, -100],
        't': [-100, -100, -100, -100, -100],
    }
    
    # Fill in the actual values from matches
    for match in matches:
        # match[0] contains the bin number, so it is used to index the lists
        index = int(match[0]) - 1
        for i, key in enumerate(kinematic_variables.keys(), 1):
            kinematic_variables[key][index] = match[i]
    
    return kinematic_variables

# Function to format and write the extracted variables to a file
def write_variables_to_file(kinematic_variables, package_number, z_bin, output_file):
    with open(output_file, 'a') as f:  # 'a' mode for appending to the file
        for key, values in kinematic_variables.items():
            variable_name = f"Q2y{package_number}z{z_bin}_{key}"
            values_str = ', '.join(map(str, values))
            f.write(f"{variable_name} = {{{values_str}}};\n")

# Main function to process a single LaTeX file
def process_latex_file(latex_file_path, package_number, output_file):
    # Extract the z_bin from the filename
    z_bin = latex_file_path.split('z')[-1][0]  # Assuming format is always like "Q2yXzY"
    
    # Read LaTeX content from the file
    with open(latex_file_path, 'r') as file:
        latex_content = file.read()
    
    # Extract kinematic variables from the LaTeX content
    kinematic_variables = extract_kinematic_variables(latex_content)
    
    # Write the extracted and formatted variables to the output file
    write_variables_to_file(kinematic_variables, package_number, z_bin, output_file)

# Example usage
latex_file_path = "/volatile/clas12/thayward/UU_validation/packages/1/clas12_analysis_software/analysis_scripts/asymmetry_extraction/output/results/kinematics_rga_fa18_inb_epX_skimmed_timeStamp_03_06_104135.txt"  # Update this to your actual file path
package_number = 1  # Update based on your directory structure or extract dynamically
output_file = "output2/analysis_statistics.txt"
process_latex_file(latex_file_path, package_number, output_file)
