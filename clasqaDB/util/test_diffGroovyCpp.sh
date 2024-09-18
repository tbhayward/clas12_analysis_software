#!/bin/bash
# run test program to cross check groovy and c++ readers

if [ -z "$QADB" ]; then
  echo "ERROR: you must source env.sh first"; exit
fi

if [ $# -lt 2 ]; then
  echo "USAGE: $0 [test num] [argument]"
  echo "- test 1 [hipo file]"
  echo "- test 2 [run number]"
  echo "- provide a 3rd argument to suppress vimdiff (batch test)"
  exit
fi

testnum=$1
arg=$2
echo "ARGUMENT = $arg"

mkdir -p ${QADB}/tmp

suffix=$arg
if [ $testnum -eq 1 ]; then
  suffix=$(echo $arg | sed 's/^.*\///g')
fi

# groovy test
echo "EXECUTE GROOVY TEST $testnum"
pushd ${QADB}/src/examples
run-groovy test${testnum}.groovy $arg > ${QADB}/tmp/groovy.${suffix}.out
popd

# c++ test
echo "EXECUTE C++ TEST $testnum"
pushd ${QADB}/srcC/examples
test${testnum}.exe $arg > ${QADB}/tmp/cpp.${suffix}.out
popd

# vimdiff
if [ $# -lt 3 ]; then
  vimdiff ${QADB}/tmp/{cpp,groovy}.${suffix}.out
else
  diff ${QADB}/tmp/{cpp,groovy}.${suffix}.out > ${QADB}/tmp/diff.${suffix}.out
  echo produced ${QADB}/tmp/diff.${suffix}.out
fi
