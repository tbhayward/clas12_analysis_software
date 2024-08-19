import csv
import sys

def calculate_total_charge(filename):
    charges = {
        'NH3': 0,
        'C': 0,
        'CH2': 0,
        'Helium Bath': 0,
        'Empty Target': 0
    }

    # Define the runs to skip
    skipped_runs = {}

    current_section = None
    run_data = []  # List to store run, target_pol, target_pol_sigma for NH3 only

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            # Skip comments
            if row[0].startswith('#'):
                if 'NH3' in row[0]:
                    current_section = 'NH3'
                elif 'Carbon' in row[0]:
                    current_section = 'C'
                elif 'CH2' in row[0]:
                    current_section = 'CH2'
                elif 'Helium Bath' in row[0]:
                    current_section = 'Helium Bath'
                elif 'Empty Target' in row[0]:
                    current_section = 'Empty Target'
                else:
                    current_section = None
                continue
            
            run_number = int(row[0])
            charge = float(row[2]+row[3])

            # Skip the specified runs
            if run_number in skipped_runs:
                continue

            # Accumulate charge for the correct section
            if current_section == 'NH3':
                charges['NH3'] += charge

                # Extract target polarization and its uncertainty (assumed to be in 5th and 6th columns)
                try:
                    target_pol = float(row[4])
                    target_pol_sigma = float(row[5])
                    run_data.append((run_number, target_pol, target_pol_sigma))
                except IndexError:
                    # Handle rows with missing columns
                    print(f"Skipping run {run_number} due to missing data in the row.")
                    continue

            elif current_section == 'Empty Target' and run_number == 16194:
                charges['Empty Target'] += charge
            elif current_section == 'Empty Target' and run_number == 16186:
                continue
            elif current_section:
                charges[current_section] += charge

    # Calculate the total charge
    total_charge = sum(charges.values())
    
    # Calculate fractions for each target
    fractions = {key: charges[key] / total_charge if total_charge > 0 else 0 for key in charges}

    return charges, total_charge, fractions, run_data

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python calculate_total_charge.py <input_file>")
        sys.exit(1)

    filename = sys.argv[1]
    charges, total_charge, fractions, run_data = calculate_total_charge(filename)
    
    print(f"Total accumulated charge: {total_charge} nC\n")
    print(f"Total accumulated charge for NH3: {charges['NH3']} nC ({fractions['NH3']:.3%} of total)")
    print(f"Total accumulated charge for Carbon: {charges['C']} nC ({fractions['C']:.3%} of total)")
    print(f"Total accumulated charge for CH2: {charges['CH2']} nC ({fractions['CH2']:.3%} of total)")
    print(f"Total accumulated charge for Helium Bath: {charges['Helium Bath']} nC ({fractions['Helium Bath']:.3%} of total)")
    print(f"Total accumulated charge for Empty Target (run 16194): {charges['Empty Target']} nC ({fractions['Empty Target']:.3%} of total)\n")

    # Generate and print the list of runs with target polarization and uncertainty for NH3 only
    run_data_list = ", ".join([f"{{{run}, {pol}, {sigma}}}" for run, pol, sigma in run_data])
    print(f"Run data for C++/Mathematica: {{{run_data_list}}}")