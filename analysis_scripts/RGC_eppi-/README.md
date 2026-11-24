# Exclusive e' p pi- analysis

Analysis scaffold alongside RGC_enpi+. This directory is for analysis of the completed data campaign; it does not replace processing scripts or claim that extraction code is implemented yet.

- analysis_note/: compilable Overleaf source, split into requested analysis sections.
- imports/datasets.json: the 19 merged ROOT files and branch conventions.
- imports/clas12_run_info_ND3_source.csv: processing-branch metadata snapshot. Do not overwrite shared development metadata from this file without reviewing it.
- Other directories are prepared for channel selection, PID, dilution, kinematics, corrections, radiation, systematics and extraction work.

Use the detected proton as p1 and FD pi- as p2. The pion production angle is phi2, not its lab azimuth p2_phi or a dihadron angle. RGB has Pt=0. Only RGC ND3 carries target polarization. Normalize the likelihood separately within each period/configuration and share the intended ratios across data samples.

The note includes provisional ND3 polarization calculations with an explicit sign-convention audit. Resolve input discrepancies and final accepted-run charge coverage before adopting those scale errors in a fit. The requested unspecified analysis/results sections remain placeholders. No MC analysis, new asymmetry extraction, or numerical MB prediction is claimed.
