# load_binning_scheme.py

from collections import namedtuple

# Define Binning at module level.
Binning = namedtuple("Binning", ["xBmin", "xBmax", "Q2min", "Q2max", "tmin", "tmax"])

def load_binning_scheme(csv_file_path):
    """
    Reads the integrated_bin_v2.csv file and returns a list of bin boundaries.
    
    The CSV file has two header lines, and then each data line looks something like:
    
      1    (0   0   0)    0.062    0.090    0.076    1.000    1.200    1.089    0.110    0.150    0.129    0.006    0.327
    
    Here the columns (after splitting on whitespace) are:
      index 0: bin number
      index 1: part of bin name (e.g., "(0")
      index 2: another part (e.g., "0")
      index 3: closing part of the bin name (e.g., "0)")
      index 4: xBmin
      index 5: xBmax
      index 6: xBavg (not used)
      index 7: Q2min
      index 8: Q2max
      index 9: Q2avg (not used)
      index 10: |tmin|
      index 11: |tmax|
      index 12: |tavg| (not used)
      etc.
      
    We extract columns 4, 5, 7, 8, 10, and 11 (with the t values made positive).
    """
    binning_list = []
    
    with open(csv_file_path, "r") as f:
        lines = f.readlines()
    
    # Skip the first two header lines
    data_lines = lines[2:]
    
    for i, line in enumerate(data_lines):
        tokens = line.split()  # splits on any whitespace
        if len(tokens) < 12:
            continue  # not enough tokens to parse a valid row
        try:
            # Based on the structure above:
            xBmin = float(tokens[4])
            xBmax = float(tokens[5])
            Q2min = float(tokens[7])
            Q2max = float(tokens[8])
            tmin  = abs(float(tokens[10]))
            tmax  = abs(float(tokens[11]))
        except Exception as e:
            print(f"Error parsing line {i+3}: {line.strip()} -> {e}")
            continue
        
        binning_list.append(Binning(xBmin, xBmax, Q2min, Q2max, tmin, tmax))
    
    return binning_list