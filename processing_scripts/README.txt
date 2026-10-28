STAGE 1: EVENT-ALIGNED MC PILOT (NOT PRODUCTION-READY)

This package contains replacement Groovy processors for the DVCS and pi0 MC
pilot, an efficient C++ FIFO-to-ROOT converter and a single-file pilot runner.

The processors retain all pairedNames()/pairedMethods() physics fields from the
latest supplied processing_mc_rec_*.groovy. Each processed HIPO event writes one
MCEvents row with gen_ and rec_ branches, validity flags and weight. No
MC::RecMatch-based join is used. Particle-level diagnostic sidecars are removed.

The .events.tsv pathname is a NAMED PIPE, not an on-disk TSV. The pilot runner
creates it, launches the converter and runs Groovy. The small selected-event
legacy .txt file is retained to check physics parity. It is not converted.

INSTALL: From your CLAS12 checkout root, create processing_scripts/rebuild_pilot
and copy all these files into that directory. Do not replace the old processors.

COMPILE: pilot_paired.py compiles convert_paired_stream.cpp with root-config.

PILOT: python3 processing_scripts/rebuild_pilot/pilot_paired.py --mode dvcs_mc --hipo /absolute/path/to/one_dvcsgen.hipo --period fa18_inb --energy 10.604 --events 10000

PILOT: python3 processing_scripts/rebuild_pilot/pilot_paired.py --mode pi0_mc --hipo /absolute/path/to/one_aaogen.hipo --period fa18_inb --energy 10.604 --events 10000

KNOWN BLOCKERS: Groovy compilation against JLab runtime, fitter correctness,
pi0 composite-photon reconstruction, output branch parity and MC truth validity
MUST be verified on ifarm. This package deliberately contains NO 20-worker
production launcher; running an unvalidated full campaign is prohibited.

IMPORTANT: The original MC processors pre-read the first HIPO event and
therefore skipped it. The pilot removes this pre-read. Validate counts.

The selected-event text remains small and is intentionally kept for comparison.
No existing campaign files are modified or deleted.
