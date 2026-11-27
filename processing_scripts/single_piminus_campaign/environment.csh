#!/bin/csh -f
# Deliberately avoid processing.csh and run-groovy: both mutate shared files.
if ( -e /etc/profile.d/modules.csh ) source /etc/profile.d/modules.csh
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module unload pymods
module load python/3.12.4
if ( $status != 0 ) exit 1
module load root/6.36.04
if ( $status != 0 ) exit 1
module load groovy/4.0.3
if ( $status != 0 ) exit 1
module load qadb/3.4.1
if ( $status != 0 ) exit 1
exec $argv:q
