#!/bin/bash

# Directory containing daily input files
input_dir="."

# Output directory for monthly averages
output_dir="./monthly_mean"
mkdir -p $output_dir

# Month start and end days (1-based index)
declare -A month_start=( [1]=0 [2]=31 [3]=59 [4]=90 [5]=120 [6]=151 )
declare -A month_end=( [1]=30 [2]=58 [3]=89 [4]=119 [5]=150 [6]=179 )

#declare -A month_start=( [1]=0 [2]=31 [3]=59 [4]=90 [5]=120 [6]=151 [7]=181 [8]=212 [9]=243 [10]=273 [11]=304 [12]=334 )
#declare -A month_end=( [1]=30 [2]=58 [3]=89 [4]=119 [5]=150 [6]=180 [7]=211 [8]=242 [9]=272 [10]=303 [11]=333 [12]=364 )

# Calculate and save monthly averages
for month in {1..6}; do
    output_file="${output_dir}/monthly_mean_${month}.nc"
    # Initialize empty temporary file list
    temp_files=()
    
    # Collect relevant files for the month
    for day in $(seq ${month_start[$month]} ${month_end[$month]}); do
        day_file=$(printf "${input_dir}/avg_NWP12_%04d.nc" $((day+1)))
        if [ -f $day_file ]; then
            temp_files+=($day_file)
        fi
    done
    
    # Merge and average the collected files
    if [ ${#temp_files[@]} -gt 0 ]; then
        ncrcat "${temp_files[@]}" temp_${month}.nc
        if [ $? -ne 0 ]; then
            echo "Error in ncrcat for month $month"
            exit 1
        fi
        ncra -O temp_${month}.nc $output_file
        if [ $? -ne 0 ]; then
            echo "Error in ncra for month $month"
            exit 1
        fi
        rm temp_${month}.nc  # Remove temporary file
        echo "Monthly average for month ${month} saved: $output_file"
    else
        echo "No files found for month ${month}"
    fi
done

