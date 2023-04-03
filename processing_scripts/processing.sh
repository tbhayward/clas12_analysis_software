#!/bin/sh

cd clasqaDB
. env.sh
cd ..
coatjava/bin/run-groovy "$1" "$2" "$3" "$4"
