#!/bin/csh

# This script creates a new workflow using the swif2 command
# provided in the '/site/bin/' directory.

# Set the alias for swif
alias swif /site/bin/swif2

# Check if the user provided a workflow name
if ($#argv != 1) then
    echo "Usage: $0 workflow_name"
    exit 1
endif

# Get the workflow name from the command line argument
set workflow_name = $argv[1]

# Create the workflow using the swif command
swif create $workflow_name

# Inform the user that the workflow was created
echo "Workflow '$workflow_name' created successfully."
