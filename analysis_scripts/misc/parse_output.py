import re

# Improved Function to extract kinematic variables and z-bin from LaTeX content
def extract_kinematic_variables_and_zbin(latex_content):
    # Updated regex pattern to correctly match the caption based on the provided LaTeX format
    caption_pattern = re.compile(r"\\caption{The mean kinematic variables in each of the bins for the extracted \$Q2y(\d)z(\d)\$ asymmetries. Values given in GeV or GeV\\$\\^2\\$ where appropriate.}")
    
    # Initialize variables
    package_number, z_bin = None, None
    
    # Attempt to find all matches for captions in the LaTeX content
    captions = caption_pattern.findall(latex_content)
    if captions:
        # Assuming the first match is representative for this script's purpose
        package_number, z_bin = captions[0]
    else:
        print("No caption match found. Check the LaTeX file format.")
        # Return early with None values to signal failure
        return None, None, None
    
    # Define regex pattern to match table rows
    pattern = re.compile(r"\\hline\s+(\d+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([\d.]+)~&~([-\d.]+)~&~([-\d.]+)~&~([-\d.]+)\\\\\s+\\hline")
    
    # Find all matches of the pattern in the LaTeX content
    matches = pattern.findall(latex_content)

    # Initialize kinematic variables with default values of -100
    kinematic_variables = { 'Q2': [-100] * 5, 'W': [-100] * 5, 'xB': [-100] * 5, 'y': [-100] * 5, 'z': [-100] * 5, 'zeta': [-100] * 5, 'PT': [-100] * 5, 'xF': [-100] * 5, 't': [-100] * 5 }
    
    # Fill in actual values from matches
    for match in matches:
        index = int(match[0]) - 1
        for i, key in enumerate(kinematic_variables.keys(), 1):
            kinematic_variables[key][index] = float(match[i])

    return kinematic_variables, package_number, z_bin

# Assuming the rest of your script remains unchanged, just update the process_latex_file function to handle None values:
def process_latex_file(latex_file_path, output_file):
    with open(latex_file_path, 'r') as file:
        latex_content = file.read()

    kinematic_variables, package_number, z_bin = extract_kinematic_variables_and_zbin(latex_content)
    if kinematic_variables and package_number and z_bin:
        write_variables_to_file(kinematic_variables, package_number, z_bin, output_file)
        print(f"Processed and wrote variables for package {package_number}, z-bin {z_bin}.")
    else:
        print(f"Failed to process file: {latex_file_path}")

# Update these paths according to your setup
latex_file_path = "/volatile/clas12/thayward/UU_validation/packages/1/clas12_analysis_software/analysis_scripts/asymmetry_extraction/output/results/kinematics_rga_fa18_inb_epX_skimmed_timeStamp_03_06_104135.txt"
output_file = "output2/analysis_statistics.txt"
process_latex_file(latex_file_path, output_file)
