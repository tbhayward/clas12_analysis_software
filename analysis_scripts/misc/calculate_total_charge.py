import csv
import sys

def calculate_total_charge(filename, nh3_filter=None):
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
            charge_pos = float(row[2])
            charge_neg = float(row[3])
            charge = charge_pos + charge_neg

            # Skip the specified runs
            if run_number in skipped_runs:
                continue

            # Accumulate charge for the correct section
            if current_section == 'NH3':
                # If a filter is set, only accumulate charges for runs in the NH3 period
                if nh3_filter and not (nh3_filter[0] <= run_number <= nh3_filter[1]):
                    continue
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

def print_charge_results(charges, total_charge, fractions, period_name="Total"):
    print(f"{period_name} accumulated charge: {total_charge:.3f} nC\n")
    print(f"{period_name} accumulated charge for NH3: {charges['NH3']:.3f} nC ({fractions['NH3']:.3%} of total)")
    print(f"{period_name} accumulated charge for Carbon: {charges['C']:.3f} nC ({fractions['C']:.3%} of total)")
    print(f"{period_name} accumulated charge for CH2: {charges['CH2']:.3f} nC ({fractions['CH2']:.3%} of total)")
    print(f"{period_name} accumulated charge for Helium Bath: {charges['Helium Bath']:.3f} nC ({fractions['Helium Bath']:.3%} of total)")
    print(f"{period_name} accumulated charge for Empty Target (run 16194): {charges['Empty Target']:.3f} nC ({fractions['Empty Target']:.3%} of total)\n")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python calculate_total_charge.py <input_file>")
        sys.exit(1)

    filename = sys.argv[1]
    
    # Define NH3 periods
    nh3_periods = {
        'Period 1': (16137, 16148),
        'Period 2': (16156, 16178),
        'Period 3': (16211, 16228),
        'Period 4': (16231, 16260),
        'Period 5': (16318, 16333),
        'Period 6': (16335, 16357),
        'Period 7': (16709, 16720),
        'Period 8': (16721, 16766),
        'Period 9': (16767, 16772)
    }

    # Calculate for the entire data set
    charges, total_charge, fractions, run_data = calculate_total_charge(filename)
    print_charge_results(charges, total_charge, fractions)

    # Calculate for each NH3 period individually
    for period, run_range in nh3_periods.items():
        charges, total_charge, fractions, _ = calculate_total_charge(filename, nh3_filter=run_range)
        print_charge_results(charges, total_charge, fractions, period_name=period)