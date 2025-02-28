import csv
from collections import namedtuple

def load_binning_scheme(csv_file_path):
    """
    Reads the integrated_bin_v2.csv file and returns a list of bin boundaries.
    The CSV is assumed to have the following columns (starting on row 3):
      bin, bin, xBmin, xBmax, xBavg, Q2min, Q2max, Q2avg, |tmin|, |tmax|, |tavg|, ...
    We extract xBmin, xBmax, Q2min, Q2max, and the absolute values of tmin and tmax.
    """
    Binning = namedtuple("Binning", ["xBmin", "xBmax", "Q2min", "Q2max", "tmin", "tmax"])
    binning_list = []
    with open(csv_file_path, "r") as f:
        reader = csv.reader(f)
        # Skip the first two rows (header information)
        header1 = next(reader)
        header2 = next(reader)
        print("Header row 1:", header1)
        print("Header row 2:", header2)
        for i, row in enumerate(reader):
            print(f"Row {i}:", row)
            try:
                xBmin = float(row[2])
                xBmax = float(row[3])
                Q2min = float(row[5])
                Q2max = float(row[6])
                # The CSV provides |tmin| and |tmax|; we ensure they are positive.
                tmin = abs(float(row[8]))
                tmax = abs(float(row[9]))
            except (ValueError, IndexError) as e:
                print(f"Skipping row {i} due to error: {e}")
                continue
            binning_list.append(Binning(xBmin, xBmax, Q2min, Q2max, tmin, tmax))
    return binning_list