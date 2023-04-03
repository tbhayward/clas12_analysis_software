#!/bin/bash
# convert json files into human-readable tables
# run this from top-level directory of `clasqaDB`
if [ -z "$QADB" ]; then
  echo "ERROR: you must source env.sh first"; exit
fi
pushd $QADB
for file in qadb/qa.*/qaTree.json; do
  run-groovy util/parseQaTree.groovy $file
done
popd
