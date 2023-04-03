#!/bin/bash
# prints list of golden files
if [ -z "$QADB" ]; then
  echo "ERROR: you must source env.sh first"; exit
fi
run-groovy ${QADB}/util/printGoldenFiles.groovy
