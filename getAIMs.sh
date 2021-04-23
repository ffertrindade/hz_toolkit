#!/usr/bin/bash
# Create AIMs list based on a Fst threshold
# Uses st.* file generated by runRealSFS.py as input

## input files
path="/media/labgenoma4/DATAPART1/fjtrindade/fst/batch001-002_oge-1_allPure_map/map1_v2" # st.* file path, without the last "/"
pop1="batch001-002-oge-1_allPure_gut"
pop2="batch001-002-oge-1_allPure_geo"
fold="unfolded"
threshold="0.10000"
chr=($(cat $path/chr_oge1.txt)) # list of chromosomes
dist=1000

## getting AIMs
for (( i=0; i<"${#chr[@]}"; i++ )); do
        echo -e "Running for ${chr[i]} (...)"

	pos=($(awk '{print $1}' $path/st.$pop1.$pop2.$fold.${chr[i]}))
	fst=($(awk '{print $2}' $path/st.$pop1.$pop2.$fold.${chr[i]}))

	for (( j=0; j<"${#pos[@]}"; j++ )); do
		if (( $(echo "${fst[j]} > $threshold" | bc -l) )); then
			echo -e "${chr[i]} ${pos[j]}" >> $pop1.$pop2.AIMs.fst$threshold.txt
		fi

	done

done

## filtering AIMs
chr=($(awk '{print $1}' $pop1.$pop2.AIMs.fst$threshold.txt | sort | uniq)) # list of chromosomes

for (( i=0; i<"${#chr[@]}"; i++ )); do
	echo -e "Running for ${chr[i]} (...)"

	pos=($(awk -v b="${chr[i]}" '$1 == b {print $2}' $pop1.$pop2.AIMs.fst$threshold.txt))
	echo -e "${chr[i]}\t${pos[0]}" >> $pop1.$pop2.AIMs.fst$threshold.filtered.txt

	check='1'
	j='0'
	while [ "$j" -lt "${#pos[@]}" ]; do

		pos1=$((${pos[j]}+$dist))
		a=$(($j+$check))

		if [ "$a" -lt "${#pos[@]}" ]; then

			if (( $(echo "${pos[a]} < $pos1" | bc -l) )); then
				((check++))
			else
				echo -e "${chr[i]}\t${pos[a]}" >> $pop1.$pop2.AIMs.fst$threshold.filtered.txt
				check=1
				j=$a
			fi

		else
			break

		fi

        done

done
