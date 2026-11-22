#!/bin/csh -f
# Run from the analysis repository root, inside tmux on either ifarm host.
if ( $#argv != 1 ) then
    echo 'Usage: csh -f processing_scripts/exclusive_pion_campaign/run_test.csh /volatile/clas12/thayward/eppi-/test_NAME'
    exit 1
endif
set campaign = "$1"
set repo = "$cwd"
set scripts = "$repo/processing_scripts/exclusive_pion_campaign"
csh -f "$scripts/environment.csh" python3 "$scripts/campaign.py" prepare --repo "$repo" --campaign "$campaign" --mode test
if ( $status != 0 ) exit 1
csh -f "$scripts/environment.csh" python3 "$scripts/campaign.py" launch --campaign "$campaign"
exit $status
