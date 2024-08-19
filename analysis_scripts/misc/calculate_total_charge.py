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
    # skipped_runs = {16317, 16742}
    skipped_runs = {3284729834729}

    nh3_abridged_charge = 0
    c_abridged_charge = 0

    current_section = None

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

            # Track the charge for the "NH3 abridged" and "C abridged" runs
            if run_number in nh3_abridged_charge:
                nh3_abridged_charge += charge
            if run_number in c_abridged_charge:
                c_abridged_charge += charge

    # Total charge using the abridged NH3 and C values
    total_charge_abridged = nh3_abridged_charge + c_abridged_charge + charges['CH2'] + charges['Helium Bath'] + charges['Empty Target']

    # Original total charge
    total_charge = sum(charges.values())
    
    # Calculate fractions using the total charge with abridged NH3 and C
    fractions_abridged = {
        'NH3': nh3_abridged_charge / total_charge_abridged,
        'C': c_abridged_charge / total_charge_abridged,
        'CH2': charges['CH2'] / total_charge_abridged,
        'Helium Bath': charges['Helium Bath'] / total_charge_abridged,
        'Empty Target': charges['Empty Target'] / total_charge_abridged
    }

    return charges, total_charge, fractions_abridged, nh3_abridged_charge, c_abridged_charge

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python calculate_total_charge.py <input_file>")
        sys.exit(1)

    filename = sys.argv[1]
    charges, total_charge, fractions_abridged, nh3_abridged_charge, c_abridged_charge = calculate_total_charge(filename)
    
    print(f"Total accumulated charge: {total_charge} nC\n")
    print(f"Total accumulated charge for NH3: {charges['NH3']} nC ({(charges['NH3'] / total_charge):.3%} of total)")
    print(f"Total accumulated charge for Carbon: {charges['C']} nC ({(charges['C'] / total_charge):.3%} of total)")
    print(f"Total accumulated charge for CH2: {charges['CH2']} nC ({(charges['CH2'] / total_charge):.3%} of total)")
    print(f"Total accumulated charge for Helium Bath: {charges['Helium Bath']} nC ({(charges['Helium Bath'] / total_charge):.3%} of total)")
    print(f"Total accumulated charge for Empty Target (run 16194): {charges['Empty Target']} nC ({(charges['Empty Target'] / total_charge):.3%} of total)\n\n")
    
    print(f"Total accumulated charge for NH3 abridged runs: {nh3_abridged_charge} nC ({fractions_abridged['NH3']:.3%} of total abridged)")
    print(f"Total accumulated charge for C abridged runs: {c_abridged_charge} nC ({fractions_abridged['C']:.3%} of total abridged)")
    print(f"Total accumulated charge for CH2: {charges['CH2']} nC ({fractions_abridged['CH2']:.3%} of total abridged)")
    print(f"Total accumulated charge for Helium Bath: {charges['Helium Bath']} nC ({fractions_abridged['Helium Bath']:.3%} of total abridged)")
    print(f"Total accumulated charge for Empty Target (run 16194): {charges['Empty Target']} nC ({fractions_abridged['Empty Target']:.3%} of total abridged)")