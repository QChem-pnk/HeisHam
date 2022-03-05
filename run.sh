#!/bin/bash

# Check number of arguments
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "At least one argument is needed."
  echo "Usage:"
  echo "${0} input [output]"
  exit
fi

if [[ ! -e inp/$1 ]]; then # Check if file exist with extension
	if [[ ! -e inp/${1}.inp ]]; then # Check if file exist without extension
		if [[ ! -e ${1} ]]; then # Check if file exist in another directory
			echo "Can't find the input file"
			exit
		else
			input="${1}"
			output="out/$(basename ${1} .inp).out"
		fi
	else
		input="inp/${1}.inp"
		output="out/${1}.out"
	fi
else
	input="inp/${1}"
	output="out/$(basename ${1} .inp).out"
fi

if [[ $# -eq 2 ]]; then # Check if an output file is given
    output="${2}"
fi

outdir="$(dirname "${output}")"

# Create output folder if it does not exist
if [[ ! -e $outdir ]]; then
    mkdir $outdir
fi

# Run program
./bin/magnetic.exe ${input} -o ${output}
