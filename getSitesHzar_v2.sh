#!/usr/bin/bash
# pegar sitios com frequencia alelica de varias localidades e ver quais sobrepoem
# rodar apos filtro de ld

#files=($(ls -l bamlist_l* | awk -F "_" '{print $2}' | sed 's/.txt//g'))
#sample_to_prune="lat23"
maf_files=$(cat maf_files.txt)
number_local="7"
sites_lat="sites_7lat"
hzar_file=($(cat ../loc_dist_v2.csv))
final_output=run001_input_7lat

# running angsd
#for (( i=0; i<"${#files[@]}"; i++ )); do
#	echo -e "######################"
#	echo -e "Running angsd for file ${files[i]} (...)"
#	echo -e "######################"

#	angsd -bam bamlist_${files[i]}.txt -GL 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 2e-6 -minMapQ 24 -minQ 20 -minInd 1 -minMaf 0.05 -doGlf 2 -out ${files[i]} -P 5
#done

# running ngsLD
#for (( i=0; i<"${#files[@]}"; i++ )); do

	#ind=$(wc -l bamlist_${files[i]}.txt)
	#site=$(zcat ${files[i]}.mafs.gz | awk 'NR>1 {print $1}' | wc -l)
	#zcat ${files[i]}.beagle.gz | awk 'NR>1 {print $1}' | sed 's/_/\t/g' > ${files[i]}.pos
#	echo -e "######################"
	#echo -e "Running ngsLD for file ${files[i]}: $ind individuals, $site sites, ${files[i]}.pos position sites (...)"
#	echo -e "######################"
	#ngsLD --geno ${files[i]}.beagle.gz --n_ind $ind --n_sites $site --pos ${files[i]}.pos --out ${files[i]}.LD --n_threads 10
#done

# running prune_graph only for the file with fewer sites (>=2 ind)
#perl /home/labgenoma4/programs/ngsLD/scripts/prune_graph.pl --in_file $sample_to_prune.LD --max_kb_dist 5 --min_weight 0.5 --out $sample_to_prune.unlinked

# juntando arquivos para achar posicoes em comum e filtrando elas por LD
#zcat $maf_files | awk 'NR>1 {print $1":"$2}' | sort | uniq -cd | awk -v i="$number_local" '$1 == i {print $2}' > $sites_lat.txt
#mv $sites_lat.txt $sites_lat.final.txt
#grep -Ff $sample_to_prune.unlinked $sites_lat.txt > $sites_lat.final.txt

# montar o arquivo dos sitios a serem usados de cada localidade
sites_final=($(cat $sites_lat.final.txt))
for (( j=0; j<"${#sites_final[@]}"; j++ )); do

	echo -e "###################"
	echo -e "Preparing files for ${sites_final[j]} analysis (...)"
        echo -e "###################"

        ale=$(echo -e ${sites_final[j]} | sed 's/:/-/g')
	echo "locID,dist,freq,samples" > $final_output.$ale.txt
	chr=$(echo -e ${sites_final[j]} | awk -F ":" '{print $1}')
	pos=$(echo -e ${sites_final[j]} | awk -F ":" '{print $2}')

	for (( i=1; i<"${#hzar_file[@]}"; i++ )); do

		echo -e "${hzar_file[i]}"
		lat=$(echo -e ${hzar_file[i]} | awk -F ";" '{print $1}')
		dist=$(echo -e ${hzar_file[i]} | awk -F ";" '{print $2}')

		freq=($(zcat ../$lat.mafs.gz | awk -v a="$chr" -v b="$pos" '$1 == a && $2 == b {print $1, $2, $4, $5, $7}'))

		if [ -z "$freq" ]; then
			echo -e "$chr $pos N 0.00 0" >> $lat.sites.txt
			freq_snp="0.00"
			samples="0"
		else
			echo -e ${freq[@]} >> $lat.sites.txt
			freq_snp=$(awk -v a="$chr" -v b="$pos" '$1 == a && $2 == b {print $4}' $lat.sites.txt)
			samples=$(awk -v a="$chr" -v b="$pos" '$1 == a && $2 == b {print $5}' $lat.sites.txt)
		fi

		echo $lat,$dist,$freq_snp,$samples >> $final_output.$ale.txt
	done

	echo "Running cline models for $ale"
	/apps/conda/fjtrindade/envs/hzarR/bin/Rscript ~/bin/hzar_par_v2.R $final_output.$ale.txt > $final_output.$ale.hzar.log

done

