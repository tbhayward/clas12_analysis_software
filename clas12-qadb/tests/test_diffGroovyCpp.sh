#!/bin/bash
# run test program to cross check groovy and c++ readers

set -e

if [ -z "$QADB" ]; then
  echo "ERROR: you must source environ.sh first"
  exit 1
fi

if [ $# -lt 2 ]; then
  echo """
USAGE: $0 [test name] [run number]

- [test name] can be any of:
$(ls src/tests | sed 's/^test//' | sed 's/\.groovy$//')

- [run number] is a run number (only used for some tests)
  """
  exit 2
fi

testname=$1
run=$2

mkdir -p ${QADB}/tmp

# groovy test
echo "EXECUTE GROOVY TEST $testname RUN $run"
pushd ${QADB}/src/tests
groovy -cp "$JYPATH" test${testname}.groovy $run > ${QADB}/tmp/groovy.${run}.out
popd

# c++ test
echo "EXECUTE C++ TEST $testname RUN $run"
pushd ${QADB}/srcC/tests
./test${testname}.exe $run > ${QADB}/tmp/cpp.${run}.out
popd

# compare
echo "DIFFERENCE:"
diff ${QADB}/tmp/{cpp,groovy}.${run}.out
exit $?
