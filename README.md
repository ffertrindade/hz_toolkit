# hz_toolkit
Scripts for hybrid zone genomic analyses

runRealSFS.py: Estimate 1dSFS, 2dSFS, folded or unfolded and use them to estimate Fst.
Usage: runRealSFS.py folded/unfolded bamlist_pop1 bamlist_pop2 anc chr_file threads

getAIMs.sh: Create AIMs list based on a Fst threshold and filter them by a minimum distance. Set the parameters inside the script.

angsd4aHMMinput.sh: Run angsd to generate maf file for the AIMs. Set the parameters inside the script.

createInputAncestryHMM.py: This script estimates allele counts for two parental panels and hybrid individuals and geerates a file that can be used as input for aHMM
Usage: createInputAncestryHMM.py mafPop1 nIndPop1 mafPop2 nIndPop2 [mafPopN nIndPopN] bamlistHyb nIndHyb outFile recRate

getSitesHzar_v2.sh: Perform cline analyses, this script runs hzar_par_v2.R. Set the parameters inside the script. Comments in Portuguese.

hzar_par_v2.R: Perform cline analyses, to be used inside getSitesHzar_v2.sh.

runNGSadmix.sh: Run NGSadmix for several k. Set the parameters inside the script.

plotNGSadmix.R: Plot NGSadmix results
Usage: Rscript plotNGSadmix.R qopt_file label_file k_num

example/: Folder with example files

