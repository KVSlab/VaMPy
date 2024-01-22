#!/bin/bash

# Array of zero-padded case numbers

declare -a cases=( "0003"  "0004" "0006" "0007" "0008" )

declare -a cases=("0003" "0004" "0006" "0007" "0008" "0009" "0019" "0020" "0021" "0023" "0027" "0028" "0031" "0032" "0034" "0035" "0074" "0076" "0077" "0078" "0081" "1029" "1030" "1031" "1032" "1033" "1035" "1037" "1039" "2022")

declare -a cases=( "0009" )

#declare -a cases=("0005" "0009" "0019" "0020" "0021" "0023" "0025" "0027" "0028" "0030" )
#declare -a cases=("0031" "0032" "0033" "0034" "0035" "0074" "0076" "0077" "0078" "0080" "0081" )
declare -a cases=( "1029" "1030" "1031" "1032" "1033" "1035" "1037" "1039" "2022")
declare -a cases=( "0029" )
declare -a cases=( "0026" )
declare -a cases=( "1038" )

declare -a cases=( "0021" "0024" "0025" "0026" "0027" "0028" "0029" "0030" "0031" "0032" "0033")
declare -a cases=( "0034" "0035" "0074" "0076"  "0078" "0080" "0081"  )

declare -a cases=( "0007" "0077" )

declare -a cases=( "0006" )
# Loop through each case
for case in "${cases[@]}"; do
    # Generate a new sbatch script for each case by replacing the placeholder
    sed "s/CASE_NUMBER/$case/g" MERGE.sh > "case_$case.sh"

    # Submit the job
    sbatch "case_$case.sh"
done

