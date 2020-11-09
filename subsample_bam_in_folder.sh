##!/bin/bash

echo "Enter number of threads to use:"
read threads
echo "Enter number of reads wanted in output file:"
read reads_wanted
echo "Delete original files?"
read delete

if [ "$delete" != "${delete#[Yy]}" ]
then
	echo -e "You stated you DO WANT to delete the file"
else
	echo -e "You stated you DONT WANT to delete the file"
fi

echo -e "Sampling BAM files in current dirctory\nto ${reads_wanted} reads with ${threads} threads...\n"

echo -e "Making 'Subsample' directory"
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
	
	echo "Sampling ${file} to ${reads_wanted} reads with seed '333'..."
	echo -e "\tRunning 'samtools view -b -@ ${threads} -s 333${ratio} -o Subsample/${reads_wanted}.${file} ${file}'"
	samtools view -@ ${threads} -bs 333${ratio} -o Subsample/${reads_wanted}.${file} ${file}
	sampled_reads=$(samtools view -c -@ ${threads} Subsample/${reads_wanted}.${file})
	echo -e "\tFinished sampling ${file} to ${sampled_reads}."
	echo -e "\tNew file has ${sampled_reads}\n"
	
	if [ "$delete" != "${delete#[Yy]}" ]
	then
		echo -e "Deleting ${file}..."
		echo -e "rm is commented to prevented deleting of files"
		# rm ${file}
	else
		echo -e "You stated you DONT WANT to delete the file"
	fi
done
