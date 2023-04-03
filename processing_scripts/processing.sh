#!/bin/sh

# Check if the current shell is csh
if expr "$SHELL" : ".*csh" > /dev/null; then
    # Run csh commands
    cd clasqaDB
    source env.csh
    cd ..
    coatjava/bin/run-groovy "$1" "$2" "$3" "$4"
    echo "$1 $2 $3 $4"
else
    # Run sh commands
    # Assuming you have an env.sh file for sh shell
    cd clasqaDB
    . env.sh
    cd ..
    coatjava/bin/run-groovy "$1" "$2" "$3" "$4"
fi
