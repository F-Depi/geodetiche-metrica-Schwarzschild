#!/bin/bash

# Check if at least two arguments are passed
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input1> <input2> [optional_inputs...]"
    exit 1
fi

# Assign the first two inputs
l=$1
E=$2

# Collect any additional inputs (if any)
shift 2  # Shift the positional parameters to the left by 2
optional_inputs="$@"

# Run the C program with both the first two inputs and any optional inputs
./main.x "$l" "$E" $optional_inputs

# Check if the C program ran successfully
if [ "$?" -ne 0 ]; then
    echo "C program failed to run."
    exit 1
fi

# Run the Python script with only the first two inputs
python3 ch2_plots.py "$l" "$E"

# Check if the Python script ran successfully
if [ "$?" -ne 0 ]; then
    echo "Python script failed to run."
    exit 1
fi

echo "Both scripts ran successfully."

