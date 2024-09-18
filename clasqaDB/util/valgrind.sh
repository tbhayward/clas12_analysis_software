#!/bin/bash
# run valgrind to test QADB.h

if [ $# -ne 1 ]; then echo "specify hipo file"; exit; fi
if [ -z "$QADB" ]; then
  echo "ERROR: you must source env.sh first"; exit
fi
mkdir -p ${QADB}/tmp

pushd ${QADB}/srcC/examples
make clean && make
valgrind \
  --log-file="${QADB}/tmp/mem.log" \
  --leak-check=full \
  --num-callers=50 \
  --verbose \
  --suppressions=$ROOTSYS/etc/valgrind-root.supp \
cutAsymmetry.exe $*

popd
