#!/bin/bash

# Define the base directory containing the yearly folders
base_dir="./"

# Define the output directory for monthly averages
output_dir="$base_dir/clim_30y_mean"
mkdir -p "$output_dir"

# List of years to process
years=$(seq 1993 2022)

# Loop over each month
for month in $(seq -w 1 12); do
    echo "Processing month $month"
    
    # Create a list to hold file paths for the current month across all years
    file_paths=()
    
    # Loop over each year
    for year in $years; do
        file_pattern="mercatorfreebiorys2v4_global_mean_${year}${month}.nc"
        file_path="$base_dir/$year/$file_pattern"
        if [ -f "$file_path" ]; then
            file_paths+=("$file_path")
        else
            echo "File not found: $file_path"
        fi
    done

    if [ ${#file_paths[@]} -gt 0 ]; then
        # Join the file paths into a single string separated by spaces
        input_files="${file_paths[@]}"
        
        # Define the output file path
        output_file="$output_dir/clim_30y_mean_${month}.nc"
        
        # Construct the NCO command
        nco_command="ncra $input_files -O $output_file"
        
        # Execute the NCO command
        echo "Executing: $nco_command"
        eval $nco_command
    else
        echo "No files found for month $month"
    fi

    echo "Finished processing month $month"
done

