#!/bin/bash
#SBATCH --job-name=probe_generator
#SBATCH --error=probe_generator_error.txt
#SBATCH --output=probe_generator_output.txt
#SBATCH --workdir=/work/jlopez2/jlopez1/L1PD
#SBATCH --mem=4000

######################################################################
# Probe Generator
#
# Authors: Juan O. Lopez    (juano.lopez@upr.edu)
#          Emanuel Martinez (emanuel.martinez8@upr.edu)
#          Javier Quinones  (javier.quinones3@upr.edu)
# License: Creative Commons Attribution-ShareAlike 4.0
# http://creativecommons.org/licenses/by-sa/4.0/legalcode
######################################################################

usage() {
	cat << EOF
Usage:

Arguments:
	-A|--aligned	Named of the file containing the aligned sequences
	-a|--aligner    Aligner software to be used.  Options are: clustalo, muscle, mafft, tcoffee
	                NOTE: The program must be installed and available in your system.
	                [Default: clustalo]
	-c|--component  LINE1 component for which probes are being generated (ORF1 or ORF2).
	                Will be used as prefix for k-mer names and as probe filtering criterion.
	                [Default: ORF1]
	-g|--genome     Reference genome in FASTA format, used for mapping the k-mers.
	-h|--help       Help (shows usage)
	-i|--input      Input file with original sequences
	-j|--join_kmers Final k-mer size to attain after joining k-mers.
	                Consecutive k-mers will be joined when the percentage of mismatches between them
					does not exceed --merge_pct argument.
	-k|--kmersize   Size of k-mer sequences to extract.
	                [Default: 100]
	-m|--metadata   Metadata in CSV format of full-length intact L1s (FLI-L1).
	-o|--output     Output file containing probes.
	                [Default: probes.fasta]
	-p|--merge_pct  Percentage of mismatches allowed between two k-mers in order to merge them
	                into a larger k-mer.
	-q|--sequence	Sequence type, currently supporting Line1 and Alu
			[Default: Line1]
	-r|--identity   Percentage range of identity.
	                [Default: 95]
	-S|--skip_aln   Skip alignment; to be used when input file is already aligned.
	-s|--skip_par	Skip parsing of Alu sequences (provide input sequences and metadata instead).
	-v|--verbose    Adds verbosity and outputs all result details.
	                [Default: non-verbose; outputs just the probes]
	-t|--threads 	Number of threads used by aligner software.
EOF
}


isInteger() {
	# https://stackoverflow.com/questions/806906/how-do-i-test-if-a-variable-is-a-number-in-bash#3951175
	case $1 in
		''|*[!0-9]*) return 1 ;;
		*) return 0 ;;
	esac
}

#############################
# Initialize default values #
#############################
aligned_file=""
aligner="clustalo"
component="ORF1"
identity_range=95
join_kmers=0
kmer_size=100
merge_pct=0
output_file="probes.fasta"
skip_aln=false
skip_par=false
seq_type="Line1"
threads=1

#####################
# Process arguments #
#####################

while [[ $# -gt 0 ]]; do
	case "$1" in
		-A|--aligned)
			aligned_file="$2"
			shift 2
			;;
		-a|--aligner)
			aligner="$2"
			shift 2
			;;
		-c|--component)
			component="$2"
			shift 2
			;;
		-g|--genome)
			ref_genome="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		-i|--input)
			input_file="$2"
			shift 2
			;;
		-j|--join_kmers)
			join_kmers="$2"
			shift 2
			;;
		-k|--kmersize)
			kmer_size="$2"
			shift 2
			;;
		-m|--metadata)
			metadata="$2"
			shift 2
			;;
		-o|--output)
			output_file="$2"
			shift 2
			;;
		-p|--merge_pct)
			merge_pct="$2"
			shift 2
			;;
        -q|--sequence)
			seq_type="$2"
			shift 2
			;;
		-r|--identity)
			identity_range="$2"
			shift 2
			;;
		-S|--skip_aln)
			skip_aln=true
			shift 1
			;;
		-s|--skip_par)
			skip_par=true
			shift 1
			;;
		-t|--threads)
			threads="$2"
			shift 2
			;;
		-v|--verbose)
			verbose="-v"
			shift 1
			;;
		*)
			echo "Invalid option: $1"
			exit 1
			;;
	esac
done

######################
# Validate arguments #
######################

if ! isInteger $kmer_size || [ $kmer_size -lt 1 ] ; then
	echo "ERROR: K-mer size must be a positive integer!"
	exit 1
fi
if ! isInteger $threads || [ $threads -lt 1 ] ; then
	echo "ERROR: Number of threads must be a positive integer!"
	exit 1
fi
if ! isInteger $identity_range || [ $identity_range -lt 1 ] ||
	[ $identity_range -gt 100 ] ; then
	echo "ERROR: Identity range must be a percentage between 1 and 100!"
	exit 1
fi
if ! isInteger $merge_pct || [ $merge_pct -lt 0 ] ||
	[ $merge_pct -gt 100 ] ; then
	echo "ERROR: Merge percentage must be a percentage between 0 and 100!"
	exit 1
fi

if [ "$skip_aln" = true ] ; then
	aligned_file=$input_file
fi

if [ "$ref_genome" = "" ] ; then
	echo "ERROR: Reference genome not provided!"
	exit 1
fi

if [ "$metadata" = "" ] && [ seq_type = "Line1" ] ; then
	echo "ERROR: Metadata CSV not provided!"
	exit 1
fi

if [ "$aligned_file" = "" ] ; then
	# If aligned filename is not provided, we create our own by using
	# the name of the input file and adding "_aligned" before the extension.
	# In this case, the aligned file will be created in pwd.
	extension="${input_file##*.}"
	nameWoExt="${input_file%.*}" # This includes the path
	aligned_file="${nameWoExt}_aligned.${extension}"
fi

###############
# Main script #
###############

# Python scripts are in same directory as this shell script
# They will be executed from that directory (make sure they are executable)

#https://stackoverflow.com/questions/56962129/how-to-get-original-location-of-script-used-for-slurm-job
if [ -n "$SLURM_JOB_ID" ] ; then
	# Need to find dir with script, since SLURM copies scripts to a different directory.
	#command_line=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
	#dir_path=$(dirname "$(echo "$command_line" | awk '{print $1}')")
	# User MUST supply workdir for this to work.
	dir_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/WorkDir=/{print $2;exit}')
else
	# Other L1PD scripts are in the same directory as this script.
	dir_path=$(dirname $0)
fi

# Filter sequences so we are only left with Alu's
# Also performs clustering of similar sequences
if [ "$seq_type" = "Alu" ] && [ "$skip_par" = false ] ; then
	$dir_path/clustering/parse_alu_sequences_cluster.py -i ${input_file} -o "alu_sequences.fasta" -m ${metadata}
	input_file="alu_sequences.fasta"
fi

# Align the sequences
if [ "$skip_aln" = false ] ; then
	$dir_path/align.py -i "$input_file" -o "$aligned_file" -p $aligner -t $threads
	rc=$?
	if [ $rc -ne 0 ] ; then
		exit $rc
	fi
fi

# Extract all k-mers that meet the specified identity range
if [ "$seq_type" = "Line1" ] ; then
	$dir_path/extract_kmers.py --inFASTA "$aligned_file" --outFASTA "${component}_${kmer_size}mers.fasta" -k $kmer_size -r $identity_range -c $component -p $merge_pct -q 1
	if [ ! -s "${component}_${kmer_size}mers.fasta" ]; then
		echo "No kmers were found. Try lowering the identity or the k-mer size."
		exit 1
	fi
else
	$dir_path/extract_kmers.py --inFASTA "$aligned_file" --outFASTA alu_${kmer_size}mers.fasta -k $kmer_size -r $identity_range -c $component -q 0
	if [ ! -s "alu_${kmer_size}mers.fasta" ]; then
		echo "No kmers were found. Try lowering the identity or the k-mer size."
		exit 1
	fi
fi

# Index the genome if necessary
if [ ! -f $ref_genome.index ] ; then
	mrfast --index $ref_genome
	rc=$?
	if [ $rc -ne 0 ] ; then
		echo "ERROR: mrfast was unable to index."
		exit $rc
	fi
fi

# Search for the k-mers within the genome
mkdir -p sam
if [ "$seq_type" = "Line1" ] ; then
	mrfast --search $ref_genome --seq "${component}_${kmer_size}mers.fasta" -o "sam/${component}_${kmer_size}mers.sam"
else
	mrfast --search $ref_genome --seq "alu_${kmer_size}mers.fasta" -o "sam/alu_${kmer_size}mers.sam"
fi
rc=$?
if [ $rc -ne 0 ] ; then
	echo "ERROR: mrfast was unable to search."
	exit $rc
fi

if [ "$seq_type" = "Line1" ] ; then
	$dir_path/kmer_probes.py "sam/${component}_${kmer_size}mers.sam" "${component}_${kmer_size}mers.fasta" "$aligned_file" "$metadata" "$component" $verbose -q 1 > "$output_file"
	# Create a directory to store our orf12 probes do nothing if it already exists
	mkdir -p orf12_probes
	# If this is the orf1 iteration empty the orf12 probe file to ensure a clean output before appending
	if [[ "$component" == "ORF1" ]] ; then
		> orf12_probes/orf12_${kmer_size}mers.fasta
	fi	
	# Append out output to the output of the next orf iteration
	cat "$output_file" >> orf12_probes/orf12_${kmer_size}mers.fasta
else
	$dir_path/kmer_probes.py "sam/alu_${kmer_size}mers.sam" "alu_${kmer_size}mers.fasta" "$aligned_file" "$metadata" "alu" -q 0 > "$output_file"
fi