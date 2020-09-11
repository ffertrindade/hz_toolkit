#!/usr/bin/bash
# rodar angsd para gerar info baseada nos AIMs

pure_pop1='batch001_pure_guttulus'
pure_pop2='batch001_pure_geoffroyi'
hibrids='batch001_hibrids'
aims='gut_ld.geo_ld.AIMs.fst0.3000' # arquivo gerado por getAIMs.sh

# index aims file
#sort -k1 $aims.txt | sed 's/ /\t/g' > $aims.sorted.txt
sed -ie 's/ /\t/g' $aims.txt
/apps/conda/pkgs/angsd-0.921-h3ef6ad9_2/bin/angsd sites index $aims.txt

# allele frequencies for the pure pops
/apps/conda/pkgs/angsd-0.921-h3ef6ad9_2/bin/angsd -bam $pure_pop1.txt -GL 2 -doMajorMinor 1 -doMaf 1 -minMapQ 24 -minQ 20 -minInd 5 -minMaf 0.05 -out $pure_pop1.freq -P 5 -sites $aims.txt
/apps/conda/pkgs/angsd-0.921-h3ef6ad9_2/bin/angsd -bam $pure_pop2.txt -GL 2 -doMajorMinor 1 -doMaf 1 -minMapQ 24 -minQ 20 -minInd 5 -minMaf 0.05 -out $pure_pop2.freq -P 5 -sites $aims.txt

# genotype likelihoods for the hybrids
/apps/conda/pkgs/angsd-0.921-h3ef6ad9_2/bin/angsd -bam $hibrids.txt -GL 2 -doMajorMinor 1 -doMaf 1 -minMapQ 24 -minQ 20 -minMaf 0.05 -doGlf 2 -out $hibrids.gl -P 5 -sites $aims.txt
