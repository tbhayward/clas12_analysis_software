#!/bin/sh

# Run the appropriate commands based on the user's shell
if echo $SHELL | grep -q "csh"; then
    # Run csh commands in a subshell
    ( 
        exec csh -c '
        cd clasqaDB
        source env.csh
        cd ..
        coatjava/bin/run-groovy '"\$1 \$2 \$3 \$4"'
        '
    )
else
    # Run sh commands
    (
        cd clasqaDB
        . env.sh
        cd ..
        coatjava/bin/run-groovy "$1" "$2" "$3" "$4"
    )
fi
