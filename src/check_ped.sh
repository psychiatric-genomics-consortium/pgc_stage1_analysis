#!/usr/bin/env bash
# This script checks pedigrees for datasets with family structure (e.g. trio datasets)
#
# Overview:
# 1) Basic variant and sample QC to prepare for LD pruning
# 2) LD prune
# 3) Run filter_ped.py, remove flagged samples
#
# Possible files created by filter_ped.py (from Picopili)
#  - ${bfile}.unmarkedparents.txt: List of related pairs that are possibly parent/offspring
#  - ${bfile}.wrongparents.txt: List of possibly incorrect parent/offspring pairs from fam file
#  - ${bfile}.remove.crossFID.txt: Samples to remove for cross-FID cryptic relatedness
#  - ${bfile}.remove.nonFID.txt: Samples unrelated to reported FIDs
#
#
# Author: Nikolas Baya, Feb 2020

bfile=${1?Missing param: No bfile prefix provided} # PLINK bfiles on which to run pedigree check

filter_ped_path=/psych/ripke/stage1_pgc/picopili/filter_ped.py 

# TODO: Set working directory

# Check that input files and filter_ped.py scripts exist
input_files=( ${bfile}.{bed,bim,fam} )
for file in ${input_files[@]}; do
	test ! -f ${file} && echo "Error: ${file} does not exist" && exit 1
done

test ! -f ${filter_ped_path} && echo "filter_ped.py path does not exist: '${filter_ped_path}'" > /dev/stderr && exit 1

log=${bfile}.check_ped.log # for redirecting all stdout, stderr

echo -e "...Performing basic variant & sample QC..." | tee -a ${log}

plink --bfile ${bfile} --silent --geno 0.05 --make-bed --out tmp0_${bfile} >> ${log} 2>&1

plink --bfile tmp0_${bfile} --silent --mind 0.02 --make-bed --out tmp1_${bfile} >> ${log} 2>&1

plink --bfile tmp1_${bfile} --silent --geno 0.02 --make-bed --out tmp2_${bfile} >> ${log} 2>&1

plink --bfile tmp2_${bfile} --silent --hwe 1e-6 --make-bed --out tmp3_${bfile} >> ${log} 2>&1

plink --bfile tmp3_${bfile} --silent --maf 0.05 --make-bed --out tmp4_${bfile} >> ${log} 2>&1

tmp4_files=( tmp4_${bfile}.{bed,bim,fam} )
for file in ${tmp4_files[@]}; do
	test ! -f ${file} && echo "Error: Basic variant and sample QC step failed." | tee -a ${log} && exit 1
done

echo -e "...LD pruning..." | tee -a ${log}

# LD pruning
ld_th=0.2
ld_wind=200
ld_move=$(($ld_wind/2))

i=1

plink --bfile tmp4_${bfile} --indep-pairwise $ld_wind $ld_move $ld_th  --silent --allow-no-sex --make-founders require-2-missing --out tmp4_${bfile}_prune$i >> ${log} 2>&1

nprune_old=$( wc -l < tmp4_${bfile}.bim )
nprune_new=$( wc -l < tmp4_${bfile}_prune$i.prune.in )

while [ ${nprune_old} -gt ${nprune_new} ]; do
	i=$(( i+1 ))
	echo "Pruning pass $i" >> ${log} 2>&1
	plink --bfile tmp4_${bfile} --extract tmp4_${bfile}_prune$((i-1)).prune.in --indep-pairwise $ld_wind $ld_move $ld_th --silent --allow-no-sex --make-founders require-2-missing --out tmp4_${bfile}_prune$i >> ${log} 2>&1

	nprune_old=$nprune_new
	nprune_new=$( wc -l < tmp4_${bfile}_prune$i.prune.in )
	echo "nprune_old: ${nprune_old}, nprune_new: ${nprune_new}" >> ${log} 2>&1
done

test ! -f tmp4_${bfile}_prune$i.prune.in  && echo "Error: LD pruning step failed." | tee -a ${log}  && exit 1

echo -e "...Extracting LD pruned set and running IBD..." | tee -a ${log}
plink --bfile tmp4_${bfile} --extract tmp4_${bfile}_prune$i.prune.in --silent --genome --make-bed --allow-no-sex --out tmp4_${bfile}_finalpruned >> ${log} 2>&1

# TODO: Check for "Segmentation fault" exit code 139 (likely due to memory constraint)

test ! -f tmp4_${bfile}_finalpruned.genome  && echo "Error: PLINK IBD (--genome) step failed." | tee -a ${log} && exit 1

echo -e "...Running filter_ped.py..." | tee -a ${log}

# TODO: Check python version before running

python ${filter_ped_path} --input-ibd tmp4_${bfile}_finalpruned.genome --bfile tmp4_${bfile}_finalpruned --format plink --out ${bfile} >> ${log} 2>&1

# Possible files created by filter_ped.py
#       - ${bfile}.unmarkedparents.txt: List of related pairs that are possibly parent/offspring
#       - ${bfile}.wrongparents.txt: List of possibly incorrect parent/offspring pairs from fam file
#       - ${bfile}.remove.crossFID.txt: Samples to remove for cross-FID cryptic relatedness
#       - ${bfile}.remove.nonFID.txt: Samples unrelated to reported FIDs

if [ $? -ne 0 ]; then
	test ! -f tmp4_${bfile}_prune$i.prune.in  && echo "Error: filter_ped.py failed." | tee -a ${log} && exit 1
fi

# remove samples in ${bfile}.{unmarked,wrong}parents.txt if neither of the pair occur in ${bfile}.remove.{cross,non}FID.txt
exist=true # expanded multiple conditional statements over two lines (checking if enough relevant files exist)
[ $( ls -l ${bfile}.{unmarked,wrong}parents.txt 2> /dev/null | wc -l ) -eq 0 ] && exist=false # if both files don't exist, set exist=False
${exist} && [ $( cat ${bfile}.remove.{cross,non}FID.txt 2>/dev/null | wc -l ) -eq 0 ] && exist=false # if at least one file existed from previous line but both files don't exist, set exist=False
if ${exist}; then # if both conditions satisfied
	cat <(tail -n +2  ${bfile}.unmarkedparents.txt 2> /dev/null | cut -f 1,2 -d ' ') \
		<(tail -n +2  ${bfile}.unmarkedparents.txt 2> /dev/null | cut -f 4,5 -d ' ') \
		<(tail -n +2  ${bfile}.wrongparents.txt 2> /dev/null | cut -f 1,2 -d ' ') \
		<(tail -n +2  ${bfile}.wrongparents.txt 2> /dev/null | cut -f 4,5 -d ' ') \
		| sort > tmp_${bfile}.unmarkedwrongparents.txt

	# get samples in tmp_${bfile}.unmarkedwrongparents.txt but not in ${bfile}.remove.{cross,non}FID.txt
	toremove=( $(diff --changed-group-format='%<' --unchanged-group-format=''  tmp_${bfile}.unmarkedwrongparents.txt \
	<(cat ${bfile}.remove.{cross,non}FID.txt 2> /dev/null \
	| grep -v "FID IID" \
	| cut -f 1,2 -d ' ' \
	| sort )) )

	rm tmp_${bfile}.unmarkedwrongparents.txt

	if [ ${#toremove[@]} -gt 0 ]; then
		for i in `seq 0 $((${#toremove[@]}/2-1))`; do
			# get rows that contain the FIDs/IIDs to remove
			grep "${toremove[$((2*i))]} ${toremove[$((2*i + 1))]}" \
				<(cat ${bfile}.{unmarked,wrong}parents.txt 2>/dev/null ) \
				| cut -f 1,2 -d ' ' >> tmp_${bfile}.unmarkedwrongparents.txt # extract only the first two columns (FID1, IID1)
		done
		# last check against individuals being removed (necessary step)
		# filters individuals who overlapped with existing *remove*FID.txt files,
		# but the other half of their parent/child pair does not overlap
		diff --changed-group-format='%<' --unchanged-group-format='' \
		<(sort tmp_${bfile}.unmarkedwrongparents.txt | uniq ) \
		<(cat ${bfile}.remove.{cross,non}FID.txt 2> /dev/null | grep -v "FID IID" | cut -f 1,2 -d ' ' | sort )  > ${bfile}.remove.unmarkedwrongparents.txt # remove duplicate FIDs/IIDs to remove
		rm tmp_${bfile}.unmarkedwrongparents.txt
	fi
fi

# remove intermediate temporary files (separated for readability, redirect stderr to /dev/null)
rm tmp1_${bfile}.irem 2> /dev/null
rm tmp[0-4]_${bfile}.{bed,bim,fam,hh,log,nosex} 2> /dev/null
rm tmp4_${bfile}_prune*.{hh,log,prune.{in,out},nosex} 2> /dev/null
rm tmp4_${bfile}_finalpruned.{bed,bim,fam,hh,log,genome} 2> /dev/null

# TODO: Check edge cases for unmarked/wrong parents

# combine files with IIDs to remove
cut -d ' ' -f 1,2 ${bfile}.remove.{{cross,non}FID,unmarkedwrongparents}.txt 2> ${log} | grep -v "FID IID" > ${bfile}.remove.check_ped.txt 

plink \
	--bfile ${bfile} \
	--remove ${bfile}.remove.check_ped.txt \
	--make-founders require-2-missing \
	--make-bed \
	--out tmp_${bfile}.check_ped >> ${log} 2>&1

outfiles_ct=$( ls -l tmp_${bfile}.check_ped.{bed,bim,fam} | wc -l )
if [ ${outfiles_ct} -eq 3 ]; then
	start_ct=$( cat ${bfile}.fam | wc -l)
	end_ct=$( cat tmp_${bfile}.check_ped.fam | wc -l)
	mv tmp_${bfile}.check_ped.bed ${bfile}.check_ped.bed
	mv tmp_${bfile}.check_ped.bim ${bfile}.check_ped.bim
	mv tmp_${bfile}.check_ped.fam ${bfile}.check_ped.fam
	rm tmp_${bfile}.check_ped.{hh,log} 2> /dev/null
	echo -e "$(( ${start_ct} - ${end_ct} )) samples removed."
	filestoremove=( ${bfile}.remove.{{cross,non}FID,unmarkedwrongparents}.txt )
	for file in ${filestoremove[@]}; do
		if [ -f ${file} ]; then
			echo "  $( cut -f 1,2 -d ' ' ${file} | grep -v "FID IID" | wc -l ) samples from ${file}"
		fi
	done
	echo "Pedigree check complete."
else
	echo "Error: check_ped.sh ran unsuccessfully." | tee -a ${log}
fi
