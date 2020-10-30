#!/usr/bin/bash
# Run NGSadmix for several k

# Files
input=~/angsd/calling/run003_gl_batch001.filtered.beagle.gz
output=~/ngsadmix/batch001_run003_filtered

# Running
echo "##### k2... #####"
NGSadmix -likes $input -K 2 -P 10 -minMaf 0.05 -seed 1 -o $output\_k2

echo "##### k3... #####"
NGSadmix -likes $input -K 3 -P 10 -minMaf 0.05 -seed 1 -o $output\_k3

echo "##### k4... #####"
NGSadmix -likes $input -K 4 -P 10 -minMaf 0.05 -seed 1 -o $output\_k4

echo "##### k5... #####"
NGSadmix -likes $input -K 5 -P 10 -minMaf 0.05 -seed 1 -o $output\_k5
