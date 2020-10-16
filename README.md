
# MiscProcessingScripts
Miscellaneous scripts written to help dealing with big data

## ICGC Subsampling of Whole Genome Sequencing BAM
Recently I have been working with trying to use some data from ICGC and have found that the size of the files is much too big to store locally on our workstation computers. To deal with this problem, I wrote a script to use the ICGC [score-client](https://docs.icgc.org/download/guide/) and [extractSoftclipped](https://github.com/dpryan79/SE-MEI).
The script takes in a couple arguments:
* -t --threads are the number of threads that will be passed to samtools view to count and sample the large bam file. *Default is **1** thread*
* -r --reads_wanted are the number of reads that are wanted in the output file. Some of the ICGC files have over 1 billion reads. *Default is **15000000** reads*
* -d --delete is whether or not to delete the full downloaded ICGC file. *Default is **n***
* -s --seed is the seed that samtools view will use when sampling the ICGC bam file. This is important to maintain repeatability  *Default is **333***
* -m --manifest is the manifest file that is downloaded from ICGC. **THIS IS REQUIRED**

ICGC requires application and granted access in order to get raw sequencing data, but there is other data formats that are available. This script has only been validated for aligned whole genome sequencing data files so it is unknown whether the other data formats (RNA-seq) will work.

TODO:
* Add ability to *not* do `extractSoftclippedRetain` from terminal
* Add selection for `fileType` for downloading non-sequencing data
