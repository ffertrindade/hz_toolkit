#!/usr/bin/bash
# gera lista de sitios AIMs de acordo com threshold de Fst
# usa de input o arquivo site.* gerado em runRealSFS.sh

path="/media/labgenoma4/DATAPART1/fjtrindade/fst" # caminho de onde esta o arquivo site.*, sem o "/"
pop1="gut_ld"
pop2="geo_ld"
threshold="0.3000"
chr=($(cat $path/chr.txt))

for (( i=0; i<"${#chr[@]}"; i++ )); do
        echo -e "Running for ${chr[i]} (...)"

	pos=($(awk '{print $1}' $path/site.$pop1.$pop2.${chr[i]}))
	fst=($(awk '{print $2}' $path/site.$pop1.$pop2.${chr[i]}))

	for (( j=0; j<"${#pos[@]}"; j++ )); do
		if (( $(echo "${fst[j]} > $threshold" | bc -l) )); then
			echo -e "${chr[i]} ${pos[j]}" >> $pop1.$pop2.AIMs.fst$threshold.txt
		fi

	done

done

