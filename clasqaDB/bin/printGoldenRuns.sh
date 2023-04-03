#!/bin/bash
# prints list of golden runs
if [ -z "$QADB" ]; then
  echo "ERROR: you must source env.sh first"; exit
fi
run-groovy ${QADB}/util/printGoldenRuns.groovy
