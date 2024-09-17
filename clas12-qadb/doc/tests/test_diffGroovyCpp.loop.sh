#!/bin/bash
# loop over runs in specified dataset, running test_diffGroovyCpp.sh on each run

set -e

if [ -z "$QADB" ]; then
  echo "ERROR: you must source environ.sh first"; exit 1
fi

if [ $# -lt 1 ]; then
  echo "USAGE: $0 [pass/dataset]"; exit 2
fi

datasetSubdir=$1

mkdir -p ${QADB}/tmp

echo """

begin TEST:
- test DumpQADB for dataset ${datasetSubdir}
- check for differences between C++ and Groovy APIs
- a diff will be dumped for any differences found, and the script will stop

"""

grep -E '^RUN: ' ${QADB}/qadb/${datasetSubdir}/qaTree.json.table |\
  awk '{print $2}' |\
  while read run; do
    ${QADB}/tests/test_diffGroovyCpp.sh DumpQADB $run
  done
