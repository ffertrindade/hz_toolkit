#!/bin/bash
# Run NGSadmix for several k

## Checking arguments
if [[ $# -lt 3 ]]; then
    echo "Usage: runNGSadmix.sh <input.beagle.gz> <output_sufix> <num_threads>"
    echo "<input.beagle.gz> beagle input file, e.g. ~/hz_toolkit/example.beagle.gz"
    echo "<output_sufix> sufix output file, e.g. ~/ngsadmix/example"
    echo "<num_threads> number of threads, e.g. 3"
    exit 1
fi

## Files and parameters
input=$1 # beagle input file, e.g. ~/hz_toolkit/example.beagle.gz
output=$2 # sufix output file, e.g. ~/ngsadmix/example
threads=$3 # number of threads, e.g. 3
minMaf="0.05"
maxiter="100000"
misTol="0.5"

## Running
echo "##### k2... #####"
NGSadmix -likes $input -K 2 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k2

echo "##### k3... #####"
NGSadmix -likes $input -K 3 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k3

echo "##### k4... #####"
NGSadmix -likes $input -K 4 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k4

echo "##### k5... #####"
NGSadmix -likes $input -K 5 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k5

echo "##### k5... #####"
NGSadmix -likes $input -K 6 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k6

echo "##### k5... #####"
NGSadmix -likes $input -K 7 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k7

echo "##### k5... #####"
NGSadmix -likes $input -K 8 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k8

echo "##### k5... #####"
NGSadmix -likes $input -K 9 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k9

echo "##### k5... #####"
NGSadmix -likes $input -K 10 -P $threads -minMaf $minMaf -maxiter $maxiter -misTol $misTol -seed 1 -o $output\_k10