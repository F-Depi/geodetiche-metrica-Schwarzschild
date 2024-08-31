#!/bin/bash

# Check if at least two arguments are passed
if [ "$#" -lt 1 ]; then
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
python3 plot_orbit.py "$l" "$E"

# Check if the Python script ran successfully
if [ "$?" -ne 0 ]; then
    echo "Python script failed to run."
    exit 1
fi

# Format l and E to the right decimal places for the filename
l_formatted=$(printf "%.3f" "$l")
E_formatted=$(printf "%.5f" "$E")
file="l${l_formatted}_E${E_formatted}.csv"

while true; do
    # Ask the user if they want to save the data
    read -p "Do you want to save the data? (y/n): " user_choice

    if [[ "$user_choice" == "n" || "$user_choice" == "no" ]]; then

        # Check if the file exists before attempting to delete it
        if [ -f "data/$file" ]; then
            rm "data/$file"
            if [ "$?" -ne 0 ]; then
                echo "Failed to delete the file."
            else
                echo "File deleted successfully."
            fi
        else
            echo "File does not exist: $file"
        fi
    elif [[ "$user_choice" == "y" || "$user_choice" == "yes" ]]; then
        read -p "name the file: " folder
        mkdir "data/keep/$folder"
        mv "data/$file" "data/keep/$folder/$file"
        echo "./main.x $l $E $optional_inputs" > "data/keep/$folder/param.txt"
        echo "Data saved in data/keep/$folder/"

    else
        echo "Invalid input. Please enter 'y' or 'n'."
        continue
    fi

    break

done

