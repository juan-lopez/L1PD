#!/bin/sh -

########################################################################
# L1PD_BAM_CRAM.sh
#
# Generate a consensus genome using the BAM/CRAM file
# and the corresponding reference genome, received as arguments.
#
# Author: Juan O. Lopez (juano.lopez@upr.edu)
# License: Creative Commons Attribution-ShareAlike 4.0
# http://creativecommons.org/licenses/by-sa/4.0/legalcode
########################################################################

usage() {
	cat << EOF
Usage:
	L1PD_BAM_CRAM.sh -f file.bam|file.cram -r reference.fasta [-p prefix] [-t threads]

Arguments:
	-f  Input file.  Either a single BAM/CRAM file. REQUIRED.
	-r  Reference genome in FASTA format. REQUIRED.
	-t  Number of additional threads.
	    [default: 0]
	-p  Filename prefix to use for all files generated.
	    [default: Filename prefix of BAM/CRAM file]

	-h  Usage
EOF
}

# https://stackoverflow.com/questions/2172352/in-bash-how-can-i-check-if-a-string-begins-with-some-value#18558871
beginswith() { case $2 in "$1"*) true;; *) false;; esac; }
endswith() { case $2 in *"$1") true;; *) false;; esac; }

while getopts hr:p:t:f: opt ; do
	case $opt in
		r) reference=$OPTARG ;;
		p) prefix=$OPTARG ;;
		t) threads=$OPTARG ;;
		f) file=$OPTARG ;;
		h) usage
			exit 0
			;;
		?) usage
			exit 1
			;;
	esac
done

# Reference and file are required arguments
if [ -z "$reference" ] || [ -z "$file" ] ; then
	usage
	exit 1
fi
# Assign default values as necessary
if [ -z $prefix ] ; then
	# If no prefix was specified, use the filename prefix (up to first
	# period) as the filename (with different extensions) for the
	# resulting files.  If prefix has a space, only first word is used.
	algnBase=$(basename $file) # Use quotes around $file to allow spaces
	prefix=${algnBase%%.*}
fi
#prefix=$3
# If threads argument not received, use 0 (bcftools' default)
#threads=${4:-0}
: ${threads:=0}

## If the file starts with ftp:// then we download it.
##if [[ $file = 'ftp://'* ]] ; then
#if beginswith "ftp://" "$file" ; then
	#echo "Downloading file with wget..."
	#wget $file
	#file=$algnBase
#fi

# Generate a Binary Variant Call Format (BCF) file
echo "BCFTOOLS MPILEUP | BCFTOOLS CALL"
cmd="bcftools mpileup -Ou -f $reference $file "
cmd="$cmd | bcftools call -mv --threads $threads -Ob -o ${prefix}.bcf"
time sh -c "$cmd"

# Normalize indels (30s)
echo "BCFTOOLS NORM"
time bcftools norm -f $reference ${prefix}.bcf --threads $threads -Ob -o ${prefix}.norm.bcf

# Index the BCF file; necessary to generate the consensus (2s)
echo "BCFTOOLS INDEX"
time bcftools index ${prefix}.norm.bcf

# Create the consensus (12s)
echo "BCFTOOLS CONSENSUS"
time bcftools consensus -f $reference ${prefix}.norm.bcf > ${prefix}.fa

# Now that the consensus is ready, it can be used to calculate the amount of patterns.
# See L1PD_Genome.sh
