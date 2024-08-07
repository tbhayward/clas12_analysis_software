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

            if current_section:
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
    print(f"Total accumulated charge for Helium Bath: {charges['Helium Bath']} nC")
    print(f"Total accumulated charge for Empty Target: {charges['Empty Target']} nC")