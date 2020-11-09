#!/bin/bash
source ~/.bash_profile

threads=1
reads_wanted=15000000
delete=n
seed=333

# inputs
while getopts t:r:d:s:m: flag
do
	case "${flag}" in
		t) threads=${OPTARG};;
		r) reads_wanted =${OPTARG};;
		d) delete=${OPTARG};;
		s) seed=${OPTARG};;
		m) manifest=${OPTARG};;
	esac
done

# colors 
red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

# echoing back inputs
echo "Number of threads to use: ${red}${threads}${reset}"
echo "Number of reads wanted in output file: ${red}${reads_wanted}${reset}"
echo "Delete original files? ${red}${delete}${reset}"
echo "Manifest file to use: ${red}${manifest}${reset}"

echo -e "Sampling BAM files in current dirctory to ${red}${reads_wanted}${reset} reads with ${red}${threads}${reset} threads...\n"

echo -e "Making ${red}'Subsample'${reset} directory"
mkdir -p Subsample
echo -e "Making temporary directory"
mkdir -p temp
echo -e "Making ${red}clippedReads${reset} foler"
mkdir -p clippedReads

cut -f3 ${manifest} > ids.txt
tail -n +2 ids.txt > ids.tmp
mv ids.tmp ids.txt

ids=$(cut -f3 ${manifest})
remove=object_id
ids_clean=${ids[@]/$remove}

for object_id in $ids_clean
do
	echo -e "Downloading ICGC object ID ${red}${object_id}${reset}..."
	score-client --profile collab download --object-id ${object_id} --output-dir temp
	echo -e "Download of ${red}${object_id}${reset} complete"
	
	echo -e "Beginning BAM Processing and Sampling.."
	FULL_BAMS=$(ls temp/*.bam | xargs -n 1 basename)
	for bam in $FULL_BAMS
	do
		echo -e "\tStarting softclipped read extraction..."
		extractSoftclippedRetain -q 10 -l 20 temp/${bam} > clippedReads/${bam}.fastq.gz &
		echo -e "\tCreating BAM index file for ${red}${bam}${reset}..."
		samtools index -b -@ ${threads} temp/${bam}
		echo -e "\t${red}${bam}${reset} index created\n"

		echo -e "\tCounting records.."
		total_reads=$(samtools view -c -@ ${threads} temp/${bam})
		echo -e "\tfound ${red}${total_reads}${reset} in ${red}${bam}${reset}."
		ratio=$(echo "scale=5;${reads_wanted}/${total_reads}" | bc)
		echo -e "\tThis is a sample ratio of ${red}0${ratio}${reset} for ${red}${reads_wanted}${reset} reads.\n"
		
		echo -e "\tSampling ${red}${bam}${reset} to ${red}${reads_wanted}${reset} reads with seed '${red}333${reset}'..."
		samtools view -@ ${threads} -bs 333${ratio} -o Subsample/${reads_wanted}.${bam} temp/${bam}
		sampled_reads=$(samtools view -c -@ ${threads} Subsample/${reads_wanted}.${bam})
		echo -e "\tFinished sampling ${red}${bam}${reset} to ${red}${reads_wanted}${reset}."
		echo -e "\t${red}${reads_wanted}.${bam}${reset} file has ${red}${sampled_reads}${reset}\n"
		
		if [ "$delete" != "${delete#[Yy]}" ]
		then
			echo -e "Deleting ${red}${bam}${reset}..."
			rm temp/${bam}
		else
			echo -e "${red}You state you DON'T WANT to delete the full BAM files${reset}"
		fi
	done	
done

rm ids.txt

: '
#only get the mapped reads from bam
samtools view -b -F 4 ${input_bam} > ${output_bam}
samtools view -f 0x40 ${output_bam} > ${read1.bam}
samtools view -f 0x80 ${output_bam} > ${read2.bam}

#getting the softclipped reads
extractSoftclippedRetain 
'
