import ROOT
import sys

ROOT.gROOT.SetBatch(True)

def get_current_from_runnum(runnum):
    # Sp18 Out: 30 nA (3211-3293), 45 nA (3867-3987)
    if 3211 <= runnum <= 3293:
        return 30
    #endif
    if 3867 <= runnum <= 3987:
        return 45
    #endif

    # Sp18 Inb: overrides, then ranges
    if runnum == 3418:
        return 70
    #endif
    if runnum == 3421 or runnum == 3422:
        return 35
    #endif
    if runnum == 3429:
        return 50
    #endif

    # 35 nA from 3306-3411
    if 3306 <= runnum <= 3411:
        return 35
    #endif

    # 50 nA from 3431-4325
    if 3431 <= runnum <= 4325:
        return 50
    #endif

    # Sp19 Inb: all 50 nA except 6616
    if runnum == 6616:
        return 5
    #endif
    if runnum >= 0:
        return 50
    #endif

    raise ValueError("Invalid run number encountered: " + str(runnum))
#endif

def main():
    if len(sys.argv) != 2:
        print("Usage: python print_event_currents.py input.root")
        sys.exit(1)
    #endif

    input_file = sys.argv[1]

    root_file = ROOT.TFile.Open(input_file)
    if not root_file or root_file.IsZombie():
        print("ERROR: could not open file " + input_file)
        sys.exit(1)
    #endif

    tree = root_file.Get("PhysicsEvents")
    if not tree:
        print("ERROR: could not find tree PhysicsEvents in " + input_file)
        sys.exit(1)
    #endif

    currents = []

    for entry in tree:
        runnum = int(entry.runnum)
        current = get_current_from_runnum(runnum)
        currents.append(str(current))
    #endfor

    print(" ".join(currents))

    root_file.Close()
#endif

if __name__ == "__main__":
    main()
#endif