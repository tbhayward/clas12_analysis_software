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

    # Placeholder run lists for "NH3 abridged" and "C abridged"
    nh3_abridged_runs = {16096, 16102, 16103, 16106, 16107, 16112, 16113, 16116, 16117, 16122, 16292, 16293, 16297, 16697, 16702}
    c_abridged_runs = {16096, 16102, 16103, 16106, 16107, 16112, 16113, 16116, 16117, 16122, 16292, 16293, 16297, 16697, 16702}

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

            if current_section == 'Empty Target' and run_number == 16194:
                charges[current_section] += charge
            elif current_section == 'Empty Target' and run_number == 16186:
                continue
            elif current_section:
                charges[current_section] += charge

            # Track the charge for the "NH3 abridged" and "C abridged" runs
            if run_number in nh3_abridged_runs:
                nh3_abridged_charge += charge
            if run_number in c_abridged_runs:
                c_abridged_charge += charge

    total_charge = sum(charges.values())
    fractions = {key: (value / total_charge) for key, value in charges.items()}

    return charges, total_charge, fractions, nh3_abridged_charge, c_abridged_charge

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python calculate_total_charge.py <input_file>")
        sys.exit(1)

    filename = sys.argv[1]
    charges, total_charge, fractions, nh3_abridged_charge, c_abridged_charge = calculate_total_charge(filename)
    
    print(f"Total accumulated charge: {total_charge} nC\n")
    print(f"Total accumulated charge for NH3: {charges['NH3']} nC ({fractions['NH3']:.3%} of total)\n")
    print(f"Total accumulated charge for Carbon: {charges['C']} nC ({fractions['C']:.3%} of total)\n")
    print(f"Total accumulated charge for CH2: {charges['CH2']} nC ({fractions['CH2']:.3%} of total)\n")
    print(f"Total accumulated charge for Helium Bath: {charges['Helium Bath']} nC ({fractions['Helium Bath']:.3%} of total)\n")
    print(f"Total accumulated charge for Empty Target (run 16194): {charges['Empty Target']} nC ({fractions['Empty Target']:.3%} of total)\n\n\n")
    
    print(f"Total accumulated charge for NH3 abridged runs: {nh3_abridged_charge} nC\n")
    print(f"Total accumulated charge for C abridged runs: {c_abridged_charge} nC")