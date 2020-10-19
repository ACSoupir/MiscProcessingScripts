#!/bin/bash

# setting defaults for the non-required inputs
threads=1
reads_wanted=15000000
delete=n
seed=333

# inputs
while getopts t:r:d:s:m: flag
do
	case "${flag}" in
		t) threads=${OPTARG};;
		r) reads_wanted=${OPTARG};;
		d) delete=${OPTARG};;
		s) seed=${OPTARG};;
		m) manifest=${OPTARG};;
	esac
done

# colors
red='tput setaf 1'
green='tput setaf 2'
reset='tput sgr0'

# echoing back inputs
echo "Number of threads to use: ${red}${threads}${reset}"
echo "Number of reads wanted in output file: ${red}${reads_wanted}${reset}"
echo "Delete original files? ${red}${delete}{reset}"

echo -e "Sampling BAM files in current directory to ${red}${reads_wanted}${reset} reads with ${red}${threads}${reset} threads...\n"

echo -e "Making ${red}'Subsample'${reset} directory"
mkdir -p Subsample

FULL_BAMs=*.bam

for file in $FULL_BAMS
do
	echo -e "\n\n"
	
	echo "Counting lines in ${red}${file}${reset}"
	echo -e "\tRunning 'samtools view -c -# ${red}${threads}${reset} ${red}${file}${reset}'"
	total_reads=$(samtools view -c -@ ${threads} ${file})
	echo -e "\tFound ${red}${total_reads}${reset} in ${red}${file}${reset}"
	ratio=$(echo "scale=4;${reads_wanted}/${total-reads}" | bc)
	echo -e "\tThis is a sample ratio of ${red}0${ratio}${reset} for ${red}${reads_wanted}${reset} reads.\n"

	echo "Sampling ${red}${file}${reset} to ${red}${reads_wanted}${reset} reads with seed '${red}${seed}${reset}'..."
	echo -e "\tRunning 'samtools view -v -@ ${red}${threads}${reset} -s ${red}333${ratio}${reset} -o Subsample/${red}${reads_wanted}.${file} ${file}${reset}'"
	samtools view -v -@ ${threads} -s 333${ratio} -o Subsample/${reads_wanted}.${file} ${file}
	sampled_reads=$(samtools view -c -@ ${threads} Subsample/${reads_wanted}.${file})
	echo -e "\tFinished sampling ${red}${file}${reset} to ${red}${sampled_reads}${reset}."
	echo -e "\tNew file has ${red}${sampled_reads}${reset}\n"

	if [ "$delete" != "${delete#[Yy]}" ]
	then
		echo -e "Deleting ${red}${file}${reset}"
		rm ${file}
	else
		echo -e "You stated you ${red}DON'T WANT${reset} to delete the files"
	fi
done
