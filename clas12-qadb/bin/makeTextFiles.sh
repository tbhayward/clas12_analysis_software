#!/bin/bash
# print out text file tables in `text/`
if [ -z "$QADB" ]; then
  echo "ERROR: you must source environ.sh first"; exit
fi
pushd $QADB

# run groovy script $1, filter out header lines, redirect to output file $2
function exe() {
  printf "\n\nEXECUTING $* ...\n\n"
  $*
  printf "\nPRODUCED $(ls -t text/*|head -n1)\n\n"
}

# format text file
function col() {
  column -t $1 > $1.tmp
  mv $1{.tmp,}
}

# build
pushd util
make clean
make
popd

# execution
exe run-groovy util/printGoldenRuns.groovy
exe run-groovy util/printGoldenFiles.groovy
exe util/printSummary.exe

# re-format
col text/summary.txt

popd
