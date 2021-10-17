#!/bin/sh -

########################################################################
# L1PD_Genome.sh
#
# Execute the L1PD algorithm on a genome, specified as argument.
#
# Author: Juan O. Lopez (juano.lopez@upr.edu)
# License: Creative Commons Attribution-ShareAlike 4.0
# http://creativecommons.org/licenses/by-sa/4.0/legalcode
########################################################################

# "Constants"
PROBESFASTA="L1PD_50mers.fasta"

# Initialize variables to default values
edit_distance=20
dist_threshold=700
min_amt_kmers=9

while getopts p:e:d:m:f:h opt ; do
	case $opt in
		p) prefix=$OPTARG ;;
		e) edit_distance=$OPTARG ;;
		d) dist_threshold=$OPTARG ;;
		m) min_amt_kmers=$OPTARG ;;
		f) file=$OPTARG ;;
		h) usage
			exit 0
			;;
	esac
done

if [ -z $file ] ; then
	echo "ERROR: Genome in FASTA format is missing!"
	exit 1
fi
if [ -z $prefix ] ; then
	# If no prefix was specified, use the filename prefix (up to first
	# period) as the filename (with different extensions) for the
	# resulting files.  If prefix has a space, only first word is used.
	# Note: Use quotes around the file name to allow spaces.
	algnBase=$(basename $file)
	prefix=${algnBase%%.*}
fi

dir_path=$(dirname $0)
# If .index file doesn't exist, then index the genome
if [ ! -f $file.index ] ; then
	mrfast --index $file
fi
mrfast --search $file --seq $dir_path/$PROBESFASTA -e $edit_distance -o probes_$prefix.sam

# Apply the L1PD algorithm to generate GFF3 and histogram
$dir_path/L1PD.py probes_$prefix.sam $dir_path/$PROBESFASTA -t $dist_threshold -m $min_amt_kmers > $prefix.gff3
