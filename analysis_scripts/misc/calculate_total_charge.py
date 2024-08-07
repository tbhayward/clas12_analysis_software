import csv
import sys

def calculate_total_charge(filename):
    charges = {
        'NH3': 0,
        'C': 0,
        'CH2': 0,
        'ET': 0,
        'Helium Bath': 0
    }

    current_section = None
    helium_bath_section = False

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
                elif 'Empty Target' in row[0]:
                    current_section = 'Helium Bath'
                    helium_bath_section = True
                else:
                    current_section = None
                    helium_bath_section = False
                continue
            
            run_number = int(row[0])
            charge = float(row[1])

            if run_number == 16194:
                charges['ET'] += charge
            elif current_section == 'Helium Bath' and run_number not in [16194, 16186]:
                charges['Helium Bath'] += charge
            elif current_section:
                charges[current_section] += charge

    return charges

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python calculate_total_charge.py <input_file>")
        sys.exit(1)

    filename = sys.argv[1]
    charges = calculate_total_charge(filename)
    
    print(f"Total accumulated charge for NH3: {charges['NH3']} nC")
    print(f"Total accumulated charge for Carbon: {charges['C']} nC")
    print(f"Total accumulated charge for CH2: {charges['CH2']} nC")
    print(f"Total accumulated charge for ET (run 16194): {charges['ET']} nC")
    print(f"Total accumulated charge for Helium Bath: {charges['Helium Bath']} nC")