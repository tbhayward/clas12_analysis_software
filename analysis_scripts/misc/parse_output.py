import re

# Improved Function to extract kinematic variables and z-bin from LaTeX content
def extract_kinematic_variables_and_zbin(latex_content):
    # Adjusted regex pattern to match the actual caption format in the uploaded file
    caption_pattern = re.compile(r"\\caption{The mean kinematic variables in each of the bins for the extracted \$Q2y(\d)z(\d)\$ asymmetries. Values given in GeV or GeV\$\\^2\$ where appropriate.}")
    
    # Initialization
    kinematic_variables, package_number, z_bin = {}, None, None
    
    # Finding the caption and extracting z-bin and package number
    caption_match = caption_pattern.findall(latex_content)
    if caption_match:
        # Extract the first (and presumably only) match
        package_number, z_bin = caption_match[0]
    else:
        print("No caption match found.")
        return None, None, None  # Safely return if no caption match is found

    # Define a regex pattern for table rows
    row_pattern = re.compile(r"\\hline\s+\d+~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([-\d.]+)~&~([-\d.]+)~&~([-\d.]+)\\\\\s+\\hline")
    matches = row_pattern.findall(latex_content)

    # Prepare variables with default values
    keys = ['Q2', 'W', 'xB', 'y', 'z', 'zeta', 'PT', 'xF', 't']
    for key in keys:
        kinematic_variables[key] = [-100, -100, -100, -100, -100]

    # Fill in the actual values
    for i, match in enumerate(matches):
        for j, value in enumerate(match):
            kinematic_variables[keys[j]][i] = float(value)

    return kinematic_variables, package_number, z_bin

# Assuming the rest of your script remains unchanged, just update the process_latex_file function to handle None values:
def process_latex_file(latex_file_path, output_file):
    with open(latex_file_path, 'r') as file:
        latex_content = file.read()

    kinematic_variables, package_number, z_bin = extract_kinematic_variables_and_zbin(latex_content)
    if kinematic_variables and package_number and z_bin:
        write_variables_to_file(kinematic_variables, package_number, z_bin, output_file)
        print(f"Successfully processed and wrote variables for package {package_number}, z-bin {z_bin}.")
    else:
        print(f"Failed to process file: {latex_file_path}")

# Update these paths according to your setup
latex_file_path = "/volatile/clas12/thayward/UU_validation/packages/1/clas12_analysis_software/analysis_scripts/asymmetry_extraction/output/results/kinematics_rga_fa18_inb_epX_skimmed_timeStamp_03_06_104135.txt"
output_file = "output2/analysis_statistics.txt"
process_latex_file(latex_file_path, output_file)
