#!/bin/bash
# convert json files into human-readable tables
# run this from top-level directory of `clas12-qadb`
if [ -z "$QADB" ]; then
  echo "ERROR: you must source environ.sh first"; exit
fi
pushd $QADB
for file in $(find -P qadb -name "qaTree.json"); do
  run-groovy util/parseQaTree.groovy $file
done
popd
