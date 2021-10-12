#!/bin/sh -

######################################################################
# LINE-1 Pattern Detection
#
# Author: Juan O. Lopez (juano.lopez@upr.edu)
# License: Creative Commons Attribution-ShareAlike 4.0
# http://creativecommons.org/licenses/by-sa/4.0/legalcode
######################################################################

usage() {
	cat << EOF
Usage: L1PD.sh (-f file1.fq file2.fq -r ref.fa | -b file.bam -r ref.fa | -g genome.fa)
               [-t threads] [-p prefix] [-e edit_distance] [-d dist_threshold] [-h]

L1PD.sh must be executed in one of three modes, each specified with a different flag:
	1) FASTQ: using paired-end reads and a reference genome (FASTA).

		L1PD.sh -f file1.fq file2.fq -r ref.fa

	2) BAM/CRAM: using an alignment file (BAM or CRAM) and a reference (FASTA).

		L1PD.sh -b file.bam -r ref.fa

	3) Genome: Using a genome that is provided in FASTA format.

		L1PD.sh -g genome.fa

	Note: FASTQ mode will transition to BAM mode, which will transition to Genome mode.

Additional arguments:
	-t  Number of additional threads
	    [default: 0]
	-e  Max. edit distance
	    [default: 20]
	-d  Distance difference threshold for establishing patterns
	    [default: 700]
	-p  Prefix text to use for all file names generated
	    [default: Filename prefix of FASTQ/BAM/CRAM file]

	-h  Usage
EOF
}
	#Note: FASTQ mode will align the input files to the reference, generate a BAM file,
	      #and proceed to BAM/CRAM mode.
			#BAM/CRAM mode will re-create the subject genome according to the reference,
			#and proceed to Genome mode.
			#Genome mode will execute the L1PD algorithm to detect LINE-1s in the genome.
		#Two paired-end FASTQ files are provided with the -f flag.
		#A reference genome in FASTA format must also be provided with the -r flag.
		#A BAM/CRAM file is provided with the -b flag.
		#A reference genome in FASTA format must also be provided with the -r flag.

isInteger() {
	# https://stackoverflow.com/questions/806906/how-do-i-test-if-a-variable-is-a-number-in-bash#3951175
	case $1 in
		''|*[!0-9]*) return 1 ;;
		*) return 0 ;;
	esac
}

# Initialize variables to default values
threads=0
edit_distance=20
dist_threshold=700
min_amt_kmers=9

###################################
# Process arguments using getopts #
###################################
# https://stackoverflow.com/questions/7529856/
# retrieving-multiple-arguments-for-a-single-option-using-getopts-in-bash#22696122
fastq_mode=0 bam_mode=0 genome_mode=0
while getopts hr:p:t:e:d:m:f:b:g: opt ; do
	case $opt in
		r) reference=$OPTARG ;;
		p) prefix=$OPTARG ;;
		t) threads=$OPTARG ;;
		e) edit_distance=$OPTARG ;;
		d) dist_threshold=$OPTARG ;;
		m) min_amt_kmers=$OPTARG ;;
		h) usage
			exit 0
			;;
		b) bam_mode=1
			file=$OPTARG
			;;
		g) genome_mode=1
			file=$OPTARG
			;;
		f) 
			fastq_mode=1
			file1=$OPTARG
			# getopts doesn't accept two arguments for a flag, so we must
			# manually check for the second required argument in FASTQ mode.
			# OPTIND has the index of the next argument, if any.
			if [ $OPTIND -le $# ] && case $(eval "echo \${$OPTIND}") in -*) false;; *) true; esac; then
				eval "file2=\${$OPTIND}"
				# "shift" OPTIND since we already used that argument
				OPTIND=$((OPTIND+1))
			else
				usage
				exit 1
			fi
			;;
	esac
done

######################
# Validate arguments #
######################
val=$(expr $fastq_mode + $bam_mode + $genome_mode)
if [ $val -ne 1 ] ; then
	echo "ERROR: You must specify one of the three possible modes!"
	exit 1
fi
if ([ $bam_mode -eq 1 ] || [ $fastq_mode -eq 1 ]) && [ "$reference" = "" ] ; then
	echo "ERROR: Reference file in FASTA format is missing!"
	exit 1
fi
if ! isInteger $threads ; then
	echo "ERROR: Amount of threads must be a positive integer!"
	exit 1
fi
if ! isInteger $edit_distance || [ $edit_distance -lt 1 ] ; then
	echo "ERROR: Edit distance must be a positive integer!"
	exit 1
fi
if ! isInteger $dist_threshold || [ $dist_threshold -lt 1 ] ; then
	echo "ERROR: Distance threshold must be a positive integer!"
	exit 1
fi
if ! isInteger $min_amt_kmers || [ $min_amt_kmers -lt 1 ] ; then
	echo "ERROR: Min. amount of k-mers must be a positive integer!"
	exit 1
fi

if [ -z $prefix ] ; then
	# If no prefix was specified, use the filename prefix (up to first
	# period) as the filename (with different extensions) for the
	# resulting files.  If prefix has a space, only first word is used.
	# Note: Use quotes around the file name to allow spaces.
	if [ $fastq_mode -eq 1 ] ; then
		algnBase=$(basename $file1)
	else
		algnBase=$(basename $file)
	fi
	prefix=${algnBase%%.*}
fi

# Testing
#echo "reference=$reference, prefix=$prefix, threads=$threads, file1=$file1, file2=$file2, ed=$edit_distance, m=$min_amt_kmers, th=$dist_threshold, file=$file"
#exit 0

################
# L1PD scripts #
################
dir_path=$(dirname $0)
if [ $fastq_mode -eq 1 ] ; then
	file=$prefix.bam
	bam_mode=1
	$dir_path/L1PD_files/L1PD_FASTQ.sh -f $file1 $file2 -r $reference -p $prefix -t $threads
	if [ $? -ne 0 ] ; then
		exit 2
	fi
fi
if [ $bam_mode -eq 1 ] ; then
	$dir_path/L1PD_files/L1PD_BAM_CRAM.sh -f $file -r $reference -p $prefix -t $threads
	file=$prefix.fa
	if [ $? -ne 0 ] ; then
		exit 2
	fi
fi
# Genome mode is always executed
$dir_path/L1PD_files/L1PD_Genome.sh -f $file -p $prefix -d $dist_threshold -e $edit_distance -m $min_amt_kmers
