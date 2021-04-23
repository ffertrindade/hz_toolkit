#!/usr/bin/python3.6
# References: https://github.com/ANGSD/angsd/issues/259, http://www.popgen.dk/angsd/index.php/Fst#Note_about_fst_for_folded_spectra
# My alternative ANGSD instalation: /apps/conda/fjtrindade/envs/angsd_v0921/bin/
# Estimate 1dSFS, 2dSFS, folded or unfolded and use them to estimate Fst

import sys
import os
import re

## parsing arguments
my_args = sys.argv
if len(my_args) == 0 or len(my_args) < 7:
	print("Usage: runRealSFS.py folded/unfolded bamlist_pop1 bamlist_pop2 anc chr_file threads")
	exit()
else:
	fold = my_args[1]
	bam1 = my_args[2] # ~/angsd/bamlist_gut_allPure.txt
	pop1 = re.sub('.txt', '', re.sub('bamlist_', '', bam1))
	bam2 = my_args[3] # ~/angsd/bamlist_geo_v2.txt
	pop2 = re.sub('.txt', '', re.sub('bamlist_', '',bam2))
	anc = my_args[4] # ~/consensos/LCO-3.fa #/home/labgenoma4/duda_grupo/references/pampas_chromium_l5000.fasta #felCat9-edited.fasta
	chr = my_args[5] #($(cat ../chr.txt))
	threads = my_args[6]

	if os.path.isfile(bam1) and os.path.isfile(bam2) and os.path.isfile(anc) and os.path.isfile(chr):
		pass
	else:
		print("File doesnt exist!")
		exit()

## function to print and run a command line
def run_and_print(cmd):
	print(cmd)
	os.system(cmd)

## function to run Fst
def run_fst(pop1, pop2, chr, fold, threads):
	with open(chr) as file:
		for row in file:

			row = row.strip()
			index = "realSFS fst index {0}.saf.idx {1}.saf.idx -sfs {0}.{1}.{2}.sfs -r {3} -P {4} -fstout {0}.{1}.{2}.{3}".format(pop1, pop2, fold, row, threads)
			stats = "realSFS fst stats {0}.{1}.{2}.{3}.fst.idx > global.{0}.{1}.{2}.{3}".format(pop1, pop2, fold, row)
			stats2 = "realSFS fst stats2 {0}.{1}.{2}.{3}.fst.idx -win 50000 -step 10000 > slw.{0}.{1}.{2}.{3}".format(pop1, pop2, fold, row)
			awk = "awk 'NR>1 {print $3, $5}'"
			plot = "{4} slw.{0}.{1}.{2}.{3} > plot.{0}.{1}.{2}.{3}".format(pop1, pop2, fold, row, awk)
			site = "realSFS fst stats2 {0}.{1}.{2}.{3}.fst.idx -win 1 -step 1 | {4} > st.{0}.{1}.{2}.{3}".format(pop1, pop2, fold, row, awk)

			# prepare the fst for easy window analysis
			if fold == "unfolded":
				run_and_print(index)
			else:
				index = index+" -fold 1"
				run_and_print(index)

			# get the global estimate
			run_and_print(stats)

			# slidingwindow file to plot per window
			run_and_print(stats2)
			run_and_print(plot)

			# slidingwindow file to plot per site
			run_and_print(site)

## first calculate unfolded saf for each population
cmd = "angsd -b {0} -anc {1} -out {2} -doSaf 1 -GL 1 -doMaf 1 -doMajorMinor 1 -P {3} -minMapQ 24 -minQ 20 -minInd 20 -sites ~/paleomix/genmap/Oge-1_final_mappability1_v2.bed".format(bam1, anc, pop1, threads)
run_and_print(cmd)
cmd = "angsd -b {0} -anc {1} -out {2} -doSaf 1 -GL 1 -doMaf 1 -doMajorMinor 1 -P {3} -minMapQ 24 -minQ 20 -minInd 14 -sites ~/paleomix/genmap/Oge-1_final_mappability1_v2.bed".format(bam2, anc, pop2, threads)
run_and_print(cmd)

## then calculate SFS priors and Fst
if fold == "unfolded":

	##calculate the unfolded 2dsfs prior by using the unfolded saf
	cmd = "realSFS {0}.saf.idx {1}.saf.idx -anc {2} -P {3} > {0}.{1}.unfolded.sfs".format(pop1, pop2, anc, threads)
	run_and_print(cmd)

	##calculate Fst
	run_fst(pop1, pop2, chr, fold, threads)

	##calculate the unfolded 1dsfs prior by using the unfolded saf; use "-bootstrap 100" to bootstrap
	cmd = "realSFS {0}.saf.idx -anc {1} -P {2} > {0}.unfolded.sfs".format(pop1, anc, threads)
	run_and_print(cmd)
	cmd = "realSFS {0}.saf.idx -anc {1} -P {2} > {0}.unfolded.sfs".format(pop2, anc, threads)
	run_and_print(cmd)

elif fold == "folded":

	##calculate the folded 2dsfs prior by using the unfolded saf
	cmd = "realSFS {0}.saf.idx {1}.saf.idx -anc {3} -fold 1 -P {2} > {0}.{1}.folded.sfs".format(pop1, pop2, threads, anc)
	run_and_print(cmd)

	##calculate Fst
	run_fst(pop1, pop2, chr, fold, threads)

	##calculate the folded 1dsfs prior by using the unfolded saf
	cmd = "realSFS {0}.saf.idx -anc {2} -fold 1 -P {1} > {0}.folded.sfs".format(pop1, threads, anc)
	run_and_print(cmd)
	cmd = "realSFS {0}.saf.idx -anc {2} -fold 1 -P {1} > {0}.folded.sfs".format(pop2, threads, anc)
	run_and_print(cmd)

else:
	print("You must choose <folded> or <unfolded> spectrum.")
	exit()

## running for all chromosomes
#realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs $pop1.$pop2.sfs -P $threads -fstout $pop1.$pop2
#realSFS fst stats $pop1.$pop2.fst.idx > global.$pop1.$pop2
#realSFS fst stats2 $pop1.$pop2.fst.idx -win 50000 -step 10000 > slidingwindow.$pop1.$pop2
#awk 'NR>1 {print $3, $5}' slidingwindow.$pop1.$pop2 > slidingwindow.$pop1.$pop2.plot

## file to plot autos histogram
#awk '{print $2}' slidingwindow.$pop1.$pop2.ChrA*.plot slidingwindow.$pop1.$pop2.ChrB*.plot slidingwindow.$pop1.$pop2.ChrC*.plot slidingwindow.$pop1.$pop2.ChrD*.plot slidingwindow.$pop1.$pop2.ChrE*.plot slidingwindow.$pop1.$pop2.ChrF*.plot > slidingwindow.$pop1.$pop2.Autos.txt
#pop1=batch001-002_allPure_gut_map1; pop2=batch001-002_allPure_geo_map1; awk '{print $2}' plot.$pop1.$pop2.unfolded.chrA* plot.$pop1.$pop2.unfolded.chrB* plot.$pop1.$pop2.unfolded.chrC* plot.$pop1.$pop2.unfolded.chrD* plot.$pop1.$pop2.unfolded.chrE* > slw.$pop1.$pop2.Autos.txt
