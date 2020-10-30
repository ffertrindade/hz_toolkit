# hz_toolkit
Scripts for hybrid zone genomic analyses

angsd4aHMMinput.sh: Run angsd to generate maf file for the AIMs

createInputAncestryHMM.py: This script estimates allele counts for two parental panels and hybrid individuals and geerates a file that can be used as input for aHMM
Usage: createInputAncestryHMM.py mafPop1 nIndPop1 mafPop2 nIndPop2 bamlistHyb nIndHyb outFile

getAIMs.sh: Create AIMs list based on a Fst threshold

getSitesHzar_v2.sh: Perform cline analyses, this script runs hzar_par_v2.R

hzar_par_v2.R: Perform cline analyses

runNGSadmix.sh: Run NGSadmix for several k

plotNGSadmix.R: Plot NGSadmix results
Usage: Rscript plotNGSadmix.R qopt_file label_file k_num

example/: Folder with example files
