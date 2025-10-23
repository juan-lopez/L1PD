from Bio import SeqIO
import numpy as np
import Levenshtein
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

def cluster():
    sequences = [str(record.seq) for record in SeqIO.parse("alu_sequences_unfiltered.fasta", "fasta")]
    n = len(sequences)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            dist_matrix[i][j] = Levenshtein.distance(sequences[i], sequences[j])

    # Convert to condensed form
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')  # or 'single', 'complete'

    # Choose a distance threshold for clustering
    threshold = 5
    clusters = fcluster(Z, threshold, criterion='distance')

    # Group sequences
    from collections import defaultdict
    cluster_dict = defaultdict(list)

    for idx, cluster_id in enumerate(clusters):
        cluster_dict[cluster_id].append(sequences[idx])

    # Print example
    with open("clusters.txt", "w") as outfile:
        max_len_seqs = 0
        max_cluster_id = 0
        for cluster_id, seqs in cluster_dict.items():
            if len(seqs) > max_len_seqs:
                max_len_seqs = len(seqs)
                max_cluster_id = cluster_id
            outfile.write(f"\nCluster_{cluster_id} {len(seqs)}\n")
            for seq in seqs:
                outfile.write(seq + "\n")
    return max_cluster_id
