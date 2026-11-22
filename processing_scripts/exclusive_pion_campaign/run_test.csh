#!/bin/csh -f
# Prepare once first, then run this command separately on each farm host.
if ( $#argv != 1 ) then
    echo 'Usage: csh -f processing_scripts/exclusive_pion_campaign/run_test.csh /volatile/clas12/thayward/eppi-/test_NAME'
    exit 1
endif
set campaign = "$1"
set scripts = "$cwd/processing_scripts/exclusive_pion_campaign"
# Invoke the prepared, frozen worker implementation, rather than mutable sources.
if ( ! -e "$campaign/manifest.json" ) then
    echo 'Prepare the shared campaign once before starting workers on either host.'
    exit 1
endif
python3 -c 'import json,sys; m=json.load(open(sys.argv[1])); sys.exit(0 if m["mode"]==sys.argv[2] else "Campaign mode mismatch")' "$campaign/manifest.json" test
if ( $status != 0 ) exit 1
csh -f "$campaign/build/source/environment.csh" python3 "$campaign/build/source/campaign.py" worker --campaign "$campaign" --workers 10
exit $status
