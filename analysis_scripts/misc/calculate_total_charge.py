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

    charges_aligned = {
        'NH3': 0,
        'C': 0,
        'CH2': 0,
        'Helium Bath': 0,
        'Empty Target': 0
    }

    charges_unaligned = {
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
                elif 'CH2' in row[0]]:
                    current_section = 'CH2'
                elif 'Helium Bath' in row[0]]:
                    current_section = 'Helium Bath'
                elif 'Empty Target' in row[0]]:
                    current_section = 'Empty Target'
                else:
                    current_section = None
                continue
            
            run_number = int(row[0])
            charge_pos = float(row[2])
            charge_neg = float(row[3])
            charge = charge_pos + charge_neg

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

                    # Accumulate charges based on alignment
                    if target_pol > 0:
                        charges_aligned['NH3'] += charge_pos
                        charges_unaligned['NH3'] += charge_neg
                    elif target_pol < 0:
                        charges_aligned['NH3'] += charge_neg
                        charges_unaligned['NH3'] += charge_pos

                except IndexError:
                    # Handle rows with missing columns
                    print(f"Skipping run {run_number} due to missing data in the row.")
                    continue

            elif current_section == 'Empty Target' and run_number == 16194:
                charges['Empty Target'] += charge
                charges_aligned['Empty Target'] += charge
                charges_unaligned['Empty Target'] += charge
            elif current_section == 'Empty Target' and run_number == 16186:
                continue
            elif current_section:
                charges[current_section] += charge
                charges_aligned[current_section] += charge
                charges_unaligned[current_section] += charge

    # Calculate the total charge
    total_charge = sum(charges.values())
    total_charge_aligned = sum(charges_aligned.values())
    total_charge_unaligned = sum(charges_unaligned.values())
    
    # Calculate fractions for each target
    fractions = {key: charges[key] / total_charge if total_charge > 0 else 0 for key in charges}
    fractions_aligned = {key: charges_aligned[key] / total_charge_aligned if total_charge_aligned > 0 else 0 for key in charges_aligned}
    fractions_unaligned = {key: charges_unaligned[key] / total_charge_unaligned if total_charge_unaligned > 0 else 0 for key in charges_unaligned}

    return charges, total_charge, fractions, run_data, charges_aligned, total_charge_aligned, fractions_aligned, charges_unaligned, total_charge_unaligned, fractions_unaligned

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python calculate_total_charge.py <input_file>")
        sys.exit(1)

    filename = sys.argv[1]
    charges, total_charge, fractions, run_data, charges_aligned, total_charge_aligned, fractions_aligned, charges_unaligned, total_charge_unaligned, fractions_unaligned = calculate_total_charge(filename)
    
    print(f"Total accumulated charge: {total_charge:.3f} nC\n")
    print(f"Total accumulated charge for NH3: {charges['NH3']:.3f} nC ({fractions['NH3']:.3%} of total)")
    print(f"Total accumulated charge for Carbon: {charges['C']:.3f} nC ({fractions['C']:.3%} of total)")
    print(f"Total accumulated charge for CH2: {charges['CH2']:.3f} nC ({fractions['CH2']:.3%} of total)")
    print(f"Total accumulated charge for Helium Bath: {charges['Helium Bath']:.3f} nC ({fractions['Helium Bath']:.3%} of total)")
    print(f"Total accumulated charge for Empty Target (run 16194): {charges['Empty Target']:.3f} nC ({fractions['Empty Target']:.3%} of total)\n")

    print(f"Total accumulated charge for NH3 with aligned polarizations: {charges_aligned['NH3']:.3f} nC ({fractions_aligned['NH3']:.3%} of total aligned)")
    print(f"Total accumulated charge for Carbon with aligned polarizations: {charges_aligned['C']:.3f} nC ({fractions_aligned['C']:.3%} of total aligned)")
    print(f"Total accumulated charge for CH2 with aligned polarizations: {charges_aligned['CH2']:.3f} nC ({fractions_aligned['CH2']:.3%} of total aligned)")
    print(f"Total accumulated charge for Helium Bath with aligned polarizations: {charges_aligned['Helium Bath']:.3f} nC ({fractions_aligned['Helium Bath']:.3%} of total aligned)")
    print(f"Total accumulated charge for Empty Target with aligned polarizations (run 16194): {charges_aligned['Empty Target']:.3f} nC ({fractions_aligned['Empty Target']:.3%} of total aligned)\n")

    print(f"Total accumulated charge for NH3 with unaligned polarizations: {charges_unaligned['NH3']:.3f} nC ({fractions_unaligned['NH3']:.3%} of total unaligned)")
    print(f"Total accumulated charge for Carbon with unaligned polarizations: {charges_unaligned['C']:.3f} nC ({fractions_unaligned['C']:.3%} of total unaligned)")
    print(f"Total accumulated charge for CH2 with unaligned polarizations: {charges_unaligned['CH2']:.3f} nC ({fractions_unaligned['CH2']:.3%} of total unaligned)")
    print(f"Total accumulated charge for Helium Bath with unaligned polarizations: {charges_unaligned['Helium Bath']:.3f} nC ({fractions_unaligned['Helium Bath']:.3%} of total unaligned)")
    print(f"Total accumulated charge for Empty Target with unaligned polarizations (run 16194): {charges_unaligned['Empty Target']:.3f} nC ({fractions_unaligned['Empty Target']:.3%} of total unaligned)\n")

    # Generate and print the list of runs with target polarization and uncertainty for NH3 only
    run_data_list = ", ".join([f"{{{run}, {pol}, {sigma}}}" for run, pol, sigma in run_data])
    # print(f"Run data for C++/Mathematica: {{{run_data_list}}}")