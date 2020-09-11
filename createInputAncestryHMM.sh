#!/usr/bin/bash
# pegar arquivo de frequencia das populacoes parentais (gerado por angsd4aHMMinput.sh) e estimar allele counts
# pegar arquivo gl dos híbridos (gerado por angsd4aHMMinput.sh) e estimar allele counts
# colunas: chr pos pop1A pop1a pop2A pop2a morgans sam1A sam1a sam2A sam2a (...)

## input files
maf_pop1_file='batch001_pure_guttulus.freq.mafs.gz'
maf_pop2_file='batch001_pure_geoffroyi.freq.mafs.gz'
beagle_file='batch001_hibrids.gl.beagle.gz'

## parameters
num_ind_pop1='21'
num_ind_pop2='16'
num_ind_hibr='14'

## comparing files sites position - to prepare aHMM file only for overlaid called sites
echo "### Selecting commom sites ###"
sites_pop1=($(zcat $maf_pop1_file | head | awk 'NR>1 {print $1"_"$2}'))
sites_pop2=($(zcat $maf_pop2_file | head | awk 'NR>1 {print $1"_"$2}'))
sites_hibr=($(zcat $beagle_file | head | awk 'NR>1 {print $1}'))

sites_tmp=($(echo ${sites_pop1[@]} ${sites_pop2[@]} | tr ' ' '\n' | sort | uniq -d))
sites_final=($(echo ${sites_tmp[@]} ${sites_hibr[@]} | tr ' ' '\n' | sort | uniq -d))
unset sites_pop1 sites_pop2 sites_hibr sites_tmp
echo ${sites_final[@]}

## selecting common sites - to prepare aHMM file only for overlaid called sites
for (( i=0; i<"${#sites_final[@]}"; i++ )); do
	a=$(echo ${sites_final[i]} | sed 's/_/ /g' | awk '{print $1}')
	b=$(echo ${sites_final[i]} | sed 's/_/ /g' | awk '{print $2}')

	zcat $maf_pop1_file | awk -v a=$a -v b=$b '$1 == a && $2 == b {print $0}' >> maf_pop1.txt
	zcat $maf_pop2_file | awk -v a=$a -v b=$b '$1 == a && $2 == b {print $0}' >> maf_pop2.txt
	zcat $beagle_file | awk -v a=${sites_final[i]} '$1 == a {print $0}' >> beagle.txt
done

## comparing alleles - to verify the minor/major order in each parental pop and hybrids
echo "### Verifying minor/major order ###"
set_ale_pop2=() # 0= 1-freq; 1= freq
set_ale_hibr=()
for (( i=0; i<"${#sites_final[@]}"; i++ )); do
	a=$(echo ${sites_final[i]} | sed 's/_/ /g' | awk '{print $1}')
	b=$(echo ${sites_final[i]} | sed 's/_/ /g' | awk '{print $2}')

	alleles_pop1=$(awk -v a=$a -v b=$b '$1 == a && $2 == b {print $3"_"$4}' maf_pop1.txt)
	alleles_pop2=$(awk -v a=$a -v b=$b '$1 == a && $2 == b {print $3"_"$4}' maf_pop2.txt)
	alleles_hibr=$(awk -v a=${sites_final[i]} '$1 == a {print $2"_"$3}' beagle.txt | sed 's/0/A/' | sed 's/1/C/' | sed 's/2/G/' | sed 's/3/T/')
	echo $alleles_pop1 $alleles_pop2 $alleles_hibr

	if [ "$alleles_pop1" = "$alleles_pop2" ]; then
		set_ale_pop2+=('0')
	else
		set_ale_pop2+=('1')
	fi

	if [ "$alleles_pop1" = "$alleles_hibr" ]; then
		set_ale_hibr+=('0')
	else
		set_ale_hibr+=('1')
	fi
done

## round function (to use latter on allele count)
round() {
	printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc)
}

## reading major freqs of each parental pop
echo "### Calculating allele counts for pop1 and pop2 ###"
freq_pop1_maj=($(awk '{print $5}' maf_pop1.txt))
freq_pop2_maj=($(awk '{print $5}' maf_pop2.txt))

## allele counts for each parental pop
allelle_count_pop1_maj=()
allelle_count_pop1_min=()
allelle_count_pop2_maj=()
allelle_count_pop2_min=()
for (( i=0; i<"${#sites_final[@]}"; i++ )); do

	freq_pop1_min=$(echo "1 - ${freq_pop1_maj[i]}" | bc)
	allelle_count_pop1_maj+=($(round ${freq_pop1_maj[i]}*$num_ind_pop1*2 0))
	allelle_count_pop1_min+=($(round $freq_pop1_min*$num_ind_pop1*2 0))

	if [ "${set_ale_pop2[i]}" = "0" ]; then # 0= 1-freq; 1= freq
		freq_pop2_min=$(echo "1 - ${freq_pop2_maj[i]}" | bc)
		allelle_count_pop2_maj+=($(round ${freq_pop2_maj[i]}*$num_ind_pop2*2 0))
		allelle_count_pop2_min+=($(round $freq_pop2_min*$num_ind_pop2*2 0))
	else
		freq_pop2_maj_tmp=$(echo "1 - ${freq_pop2_maj[i]}" | bc)
		allelle_count_pop2_maj+=($(round $freq_pop2_maj_tmp*$num_ind_pop2*2 0))
		allelle_count_pop2_min+=($(round ${freq_pop2_maj[i]}*$num_ind_pop2*2 0))
	fi
done

## allele counts for hibrid individuals

allelle_count_hibr_maj
allelle_count_hibr_min
