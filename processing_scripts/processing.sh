#!/bin/sh

# Check if the current shell is csh
if echo $SHELL | grep -q "csh"; then
    # Run csh commands
    cd clasqaDB
    source env.csh
    cd ..
    coatjava/bin/run-groovy "$1" "$2" "$3" "$4"
elif echo $SHELL | grep -q "bash"; then
    # Run sh commands
    cd clasqaDB
    . env.sh
    cd ..
    coatjava/bin/run-groovy "$1" "$2" "$3" "$4"
else
    echo "Currently only c-shell and bash supported."
fi
