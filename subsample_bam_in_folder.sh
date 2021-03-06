#!/bin/bash

function Help() {
   # Display Help
   echo "This script is to help subsample bam files that are in the current folder based on the inputs below."
   echo
   echo "Syntax: bash subsample_bam_in_foldre.sh [-t|r|d|s]"
   echo "options:"
   echo "t     Select number of threads to use with SAMTools. (default 1)"
   echo "r     The number of reads desired in the output BAM file. (default 5000000"
   echo "d     Whether or not to delete the original BAM file. (default (n)o)"
   echo "s     Seed used to subsampling (reproducible). (default 333)"
   echo
}

threads=1
reads_wanted=5000000
delete=n
seed=333

if [ "$1" = "--help" ] ; then
	Help
	exit 0
fi

if [ "$1" = "-h" ] ; then
	Help
	exit 0
fi

# inputs
while getopts t:r:d:s:h: flag
do
	case "${flag}" in
		t) threads=${OPTARG};;
		r) reads_wanted=${OPTARG};;
		d) delete=${OPTARG};;
		s) seed=${OPTARG};;
	esac
done

(
if [ "$delete" != "${delete#[Yy]}" ]
then
	echo "deleting files"
	echo "${threads} threads"
	echo "seed ${seed}"
	echo "${reads_wanted} final reads wanted"
else
	echo "keeping files"
	echo "${threads} threads"
	echo "seed ${seed}"
	echo "${reads_wanted} final reads wanted"
fi

echo -e "Sampling BAM files in current dirctory\nto ${reads_wanted} reads with ${threads} threads...\n"

echo "Making 'Subsample' directory"
mkdir -p Subsample

FULL_BAMS=*.bam

for file in $FULL_BAMS
do
	echo -e "\n\n"

	echo "Counting lines in ${file}..."
	echo -e "\tRunning 'samtools view -c -@ ${threads} ${file}'"
	total_reads=$(samtools view -c -@ ${threads} ${file})
	echo -e "\tfound ${total_reads} in ${file}."
	ratio=$(echo "scale=5;${reads_wanted}/${total_reads}" | bc)
	echo -e "\tThis is a sample ratio of 0${ratio} for ${reads_wanted} reads.\n"
	
	echo "Sampling ${file} to ${reads_wanted} reads with seed '${seed}'..."
	echo -e "\tRunning 'samtools view -b -@ ${threads} -s ${seed}${ratio} -o Subsample/${reads_wanted}.${file} ${file}'"
	samtools view -@ ${threads} -bs ${seed}${ratio} -o Subsample/${reads_wanted}.${file} ${file}
	sampled_reads=$(samtools view -c -@ ${threads} Subsample/${reads_wanted}.${file})
	echo -e "\tFinished sampling ${file} to ${sampled_reads}."
	echo -e "\tNew file has ${sampled_reads}\n"
	
	if [ "$delete" != "${delete#[Yy]}" ]
	then
		echo "Deleting ${file}..."
		echo "rm is commented to prevented deleting of files"
		# rm ${file}
	else
		echo "You stated you DONT WANT to delete the file"
	fi
done
) | tee Subsampling_BAM_in_Folder.log