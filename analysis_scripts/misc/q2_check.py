#!/usr/bin/env python3

import sys
import ROOT

def count_q2(path, label):
    f = ROOT.TFile.Open(path, "READ")
    if f is None or f.IsZombie():
        raise RuntimeError("Could not open {}".format(path))
    #endif

    t = f.Get("PhysicsEvents")
    if t is None:
        raise RuntimeError("Could not find PhysicsEvents in {}".format(path))
    #endif

    if t.GetBranch("Q2") is None:
        raise RuntimeError("{} is missing Q2".format(label))
    #endif

    n_total = 0
    n_lt2 = 0
    n_ge2 = 0

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        q2 = float(getattr(t, "Q2"))

        n_total += 1

        if q2 < 2.0:
            n_lt2 += 1
        else:
            n_ge2 += 1
        #endif
    #endfor

    print("")
    print(label)
    print("  entries      = {}".format(n_total))
    print("  Q2 < 2       = {}".format(n_lt2))
    print("  Q2 >= 2      = {}".format(n_ge2))
    print("  frac Q2 < 2  = {:.6f}".format(float(n_lt2) / float(n_total) if n_total > 0 else 0.0))
    print("  frac Q2 >= 2 = {:.6f}".format(float(n_ge2) / float(n_total) if n_total > 0 else 0.0))

    f.Close()
#enddef

def main():
    if len(sys.argv) != 3:
        print("Usage: python q2_check.py data.root rec_mc.root")
        sys.exit(1)
    #endif

    count_q2(sys.argv[1], "data")
    count_q2(sys.argv[2], "reconstructed MC")
#enddef

if __name__ == "__main__":
    main()
#endif