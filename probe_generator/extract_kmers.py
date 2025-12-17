#!/usr/bin/env python3
"""This script extracts all k-mers from a consensus sequence of an input FASTA file.
The input file is assumed to be aligned, and BioPython is used to generate a
consensus sequence using a particular threshold to obtain "strong" consensus k-mers.
"""

import argparse
from collections import defaultdict
from Bio.Align import AlignInfo
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq, MutableSeq
from merge.algorithm import kmer_frequency, split_kmers, heap_merge_kmers, overlap_kmers, overlapping_kmer_frequency
import heapq

KMER_FREQ_STEP = 50

def position_consensus(sequences, position):
    """
    This function returns the consensus base from a especific position on a aligned object.
    :param sequences: This is a AlignIO object with the aligned FASTA sequences and IDs
    :param position: This is the position where the consensus
    :return: The function returns the consensus at the especifide position
    """
    if position > sequences.get_alignment_length() - 1:
        return None  # the position is out of the sequences length

    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for sequence in sequences:
        base = sequence.seq[position]
        if base in counts:
            counts[base] += 1

    consensus_base = max(counts, key=counts.get)
    return consensus_base


def split_consensus(consensus):
    """
    This function split the consensus into "pieces" according to whenever an X was found
    and stores the relative position of the pieces
    :param consensus: Consensus string
    :return: List of tuples containing relative position and piece
    """

    consensus = list(consensus)
    new_consensus = []

    # Removing consecutive X's from consensus
    i = 0
    while i < len(consensus):
        if consensus[i] != 'X':
            new_consensus.append(consensus[i])
            i += 1

        else:
            new_consensus.append('X')
            while i < len(consensus) and consensus[i] == 'X':
                i += 1

    # Adding an X to the end so the algorithm below doesn't skip the last piece
    if new_consensus[-1] != 'X':
        new_consensus.append('X')

    # Split the consensus into "pieces" according to whenever an X was found
    # Generates a tuple per piece of the form (relative position, piece)
    last_piece_pos = 0
    pieces = []

    i = 0
    while i < len(new_consensus):
        if new_consensus[i] == 'X':
            seq_obj = Seq(''.join(new_consensus[last_piece_pos:i]))
            pieces.append(((last_piece_pos), seq_obj))
            last_piece_pos = i + 1
        i += 1

    return pieces


# def extract_kmers(inName, outName, k, percentage, identityPct, component, j_size, sequence_type):
def extract_kmers(inName, outName, k, percentage, identityPct, component, sequence_type):
    align = AlignIO.read(inName, "fasta")

    identityPct = identityPct / 100  # Convert from % to float
    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.dumb_consensus(threshold=identityPct, ambiguous='X')
    consensus_mutable = MutableSeq(str(consensus))
    consensus_kmers = []  # Save a list of (position,size) tuples for k-mers that satisfy a certain threshold
    consensus_non_kmers = []  # Save a list of (position,size) tuples for sequences of non-base characters such as X and N
    threshold = k  # Threshold through which we consider if a k-mer sequence is acceptable in the consensus
    # threshold = 0
    if percentage != 0:
        percentage = (
                    percentage / 100)  # Percentage of non-base characters we allow in a union of consecutive k-mer sequences split with non-base sequences in between
        pos = 1  # Counter to keep track of k-mer pos through consensus iteration
        counter = 1  # Counter to keep track of k-mer size through consensus iteration

    if percentage != 0:
        for i in range(1, len(consensus)):
            pos += 1

            if consensus[i - 1] in ['A', 'G', 'C', 'T']:
                if consensus[i] in ['A', 'G', 'C', 'T']:
                    counter += 1  # Only add to counter if current and past chars are bases
                else:
                    # Comment out the treshold to retrieve all the bases wip
                    # if counter >= threshold:
                    consensus_kmers.append(
                        (pos - counter, counter))  # If non-base is found, check threshold and append to list
                    counter = 1  # Restart counter for a new sequence of non-bases
            elif consensus[i - 1] not in ['A', 'G', 'C', 'T']:
                if consensus[i] not in ['A', 'G', 'C', 'T']:
                    counter += 1  # Only add to counter if current and past chars are non-bases
                else:
                    consensus_non_kmers.append((pos - counter, counter))  # If a base is found, append to non-bases list
                    counter = 1  # Restart counter for a new sequence of bases

        # Comment out the treshold to retrieve all the bases wip
        # if counter >= threshold:  # Consider any final sequece left after loop is over
        if consensus[-1] in ['A', 'G', 'C', 'T']:
            consensus_kmers.append((pos - counter + 1, counter))
        if consensus[-1] not in ['A', 'G', 'C', 'T']:
            consensus_non_kmers.append((pos - counter + 1, counter))

        print("BASE SEQUENCES")
        print("(pos, size)")
        for seq in consensus_kmers:
            print(seq)
        print("(pos, size)")
        print("NON-BASE SEQUENCES")
        for seq in consensus_non_kmers:
            print(seq)

        # Merge multiple gaps into contiguous base regions

        # New heap merge
        final_kmers = heap_merge_kmers(consensus_kmers, consensus_non_kmers, percentage)

        # Filter by k-mer threshold
        final_kmers = [(start, length) for start, length in final_kmers if length >= k]

        # Print merged sequences for debugging / verification
        print("\nMERGED SEQUENCES:")
        print("BASE SEQUENCES")
        print("(pos, size)")
        for start, length in final_kmers:
            seq = consensus_mutable[start: start + length]
            print(f"Start={start}, Length={length}, Seq={str(seq)}")

        final_overlapping_kmers = final_kmers
        print("Non overlapping kmers generated")
        # TODO add param for exploratory frequency display
        for start, length in kmer_frequency(final_kmers, KMER_FREQ_STEP).items():
            print(f'{start} : {length}')

        final_kmers = split_kmers(final_kmers, k, component)
        # print(final_kmers)
        # print(len(final_kmers))

        print("Overlapping kmers generated")
        for start, length in overlapping_kmer_frequency(final_overlapping_kmers, KMER_FREQ_STEP).items():
            print(f'{start} : {length}')

        final_overlapping_kmers = overlap_kmers(final_overlapping_kmers, k, component)
        # print(final_overlapping_kmers)
        # print(len(final_overlapping_kmers))



    # Split the consensus into "pieces" according to whenever an X was found
    if sequence_type:
        pieces = consensus_mutable.split("X")
    else:
        pieces = split_consensus(consensus_mutable)

    lengths_count = {}
    ns_lengths_count = {}  # Dictionary to store lengths of consecutive Ns
    # Enter each individual piece
    with open(outName, "w") as fhOut:
        if percentage == 0:
            count = 1
            for i in range(len(pieces)):
                rel_pos = 0
                if sequence_type:
                    piece = pieces[i]

                else:
                    rel_pos, piece = pieces[i]

                length = len(piece)  # Find length for that piece
                if length == 0:
                    continue

                # Count consecutive Ns and store in dictionary
                ns_length = 0
                for base in piece:
                    if base == 'N' or base == 'n':
                        ns_length += 1
                    else:
                        if ns_length > 0:
                            if ns_length in ns_lengths_count:
                                ns_lengths_count[ns_length] += 1
                            else:
                                ns_lengths_count[ns_length] = 1
                            ns_length = 0

                # la N al final no se cuenta por que no entraba al else
                if ns_length > 0:
                    if ns_length in ns_lengths_count:
                        ns_lengths_count[ns_length] += 1
                    else:
                        ns_lengths_count[ns_length] = 1

                if length in lengths_count:  # verify if the length is already a key in the dictionary
                    lengths_count[length] += 1
                else:
                    lengths_count[length] = 1

                lCount = 0  # Initialize lCount for next part

                if component != "":
                    prefixString = ">" + component + "_" + str(k) + "mers_"
                else:
                    prefixString = ">" + str(k) + "mers_"

                while lCount < length:  # Try to extract all possible k-mers from this piece
                    pi = piece[lCount:]  # Cut substring of piece starting at lCount

                    relative_pos = count if sequence_type else rel_pos + lCount

                    if len(pi) >= k:  # If length of that subpiece is >= k
                        fhOut.write(prefixString + str(
                            relative_pos) + "\n")  # not sure if adding lCount to relative position is the right move
                        count += 1
                        fhOut.write(str(pi[0:k]) + "\n")  # Cut out k-mer from index 0
                    lCount += 1

        else:
            for prefix, start, length in final_kmers:
                fhOut.write(">" + prefix + "\n")
                fhOut.write(str(consensus_mutable[start: start + length]) + "\n")
            with open("overlapping_merged_kmers.fa","w") as output:
                for prefix, start, length in final_overlapping_kmers:
                    output.write(">" + prefix + "\n")
                    output.write(str(consensus_mutable[start: start + length]) + "\n")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate k-mers from aligned columns with specified % of identity")
    parser.add_argument("--inFASTA", help="Input FASTA file with aligned sequences")
    parser.add_argument("--outFASTA", help="Output FASTA file that will store the k-mers")
    parser.add_argument("-k", type=int, default=50, help="k-mer size")
    # This -p is what decides if a gap should be joined say ATCG XXX ATCG if it meets the percentage requirement it
    # will replace XXX with the most common actual base in that column say ATCG GTAC ATCG
    parser.add_argument("-p", type=int, default=0, help="Non-base percentage allowed within a union of k-mer sequences")
    # Identity % is in charge of the consensus between columns i.e the ammount of bases we are willing to accept as -
    # (empty bases) before we decide it is too ambigous and must be replaced as a X. This affects our unions later on.
    parser.add_argument("-r", type=int, default=95, help="Identity %")
    parser.add_argument("-c", help="String prefix for k-mer name in FASTA output")
    parser.add_argument("-q", type=int, default=1, help="Type of sequence: LINE or Other")
    args = parser.parse_args()
    extract_kmers(args.inFASTA, args.outFASTA, args.k, args.p, args.r, args.c, args.q)
