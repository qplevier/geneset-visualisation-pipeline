#!/bin/bash

folder_path=".snakemake/log"  # replace with your folder path

# Get the modification date of the most recent file
latest_date=$(ls -ltr --time-style=long-iso "$folder_path" | tail -n 1 | awk '{print $6}')

# List all files
for file in "$folder_path"/*
do
    # Get the modification date of the file
    file_date=$(ls -ltr --time-style=long-iso "$file" | awk '{print $6}')

    # If the file was not modified on the same day as the most recent file, print its name
    if [ "$file_date" != "$latest_date" ]
    then
        rm "$file"
    fi
done
