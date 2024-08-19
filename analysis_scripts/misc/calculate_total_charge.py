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
    skipped_runs = {16317, 16742}

    current_section = None
    run_data = []  # List to store run, target_pol, target_pol_sigma

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0].startswith('#'):
                if 'NH3' in row[0]:
                    current_section = 'NH3'
                elif 'Carbon Runs' in row[0]:
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
            charge = float(row[1])

            # Skip the specified runs
            if run_number in skipped_runs:
                continue

            if current_section == 'Empty Target' and run_number == 16194:
                charges[current_section] += charge
            elif current_section == 'Empty Target' and run_number == 16186:
                continue
            elif current_section:
                charges[current_section] += charge

            # Extract target polarization and its uncertainty (assumed to be in 5th and 6th columns)
            if row[0].startswith('#'):
                target_pol = float(row[4])
                target_pol_sigma = float(row[5])
                run_data.append((run_number, target_pol, target_pol_sigma))

    # Calculate the total charge
    total_charge = sum(charges.values())
    
    # Calculate fractions for each target
    fractions = {
        'NH3': charges['NH3'] / total_charge,
        'C': charges['C'] / total_charge,
        'CH2': charges['CH2'] / total_charge,
        'Helium Bath': charges['Helium Bath'] / total_charge,
        'Empty Target': charges['Empty Target'] / total_charge
    }

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

    # Generate and print the list of runs with target polarization and uncertainty
    run_data_list = ", ".join([f"{{{run}, {pol}, {sigma}}}" for run, pol, sigma in run_data])
    print(f"Run data for C++/Mathematica: {{{run_data_list}}}")