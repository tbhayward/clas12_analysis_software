#!/usr/bin/env python3

import ROOT
import math
import sys

path = sys.argv[1]
f = ROOT.TFile.Open(path, "READ")
t = f.Get("PhysicsEvents")

n = min(20, t.GetEntries())

for i in range(n):
    t.GetEntry(i)

    weight = float(getattr(t, "weight"))

    vals = []
    for name in ["x", "Q2", "W", "e_p"]:
        if t.GetBranch(name):
            vals.append("{}={:.6g}".format(name, float(getattr(t, name))))

    if t.GetBranch("e_p"):
        nu_106 = 10.604 - float(getattr(t, "e_p"))
        vals.append("nu_10p604={:.6g}".format(nu_106))

    print("entry {:4d}  weight={:.6g}  {}".format(i, weight, "  ".join(vals)))

f.Close()