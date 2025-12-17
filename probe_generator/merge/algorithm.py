# class structure wip
from collections import defaultdict
import heapq


class ColRange:
    """
    Represents a contiguous region of bases in the consensus.

    position           = start index in consensus
    length             = number of bases
    next               = next ColRange (the next contiguous base region)
    ambiguity_with_next = accumulated ambiguity % if merged with next
    gap_len            = number of ambigous bases X/N between this and next
    ambiguity_so_far    = ambiguity of the entire sequence including all merges performed so far
    all_ambigous_base_len    = length of all the ambigous bases we have merged in to this sequence
    """

    def __init__(self, position, length):
        self.position = position
        self.length = length
        self.next = None
        self.ambiguity_with_next = float("inf")  # default until computed
        self.gap_len = 0
        self.ambiguity_so_far = float("inf")  # default until computed
        self.all_ambigous_base_len = 0

    def update_ambiguity(self):
        """
        Computes relative ambiguity if this block merges with its next block.
        """
        if self.next is None:
            self.ambiguity_with_next = float("inf")
            return

        # DISCUSS
        # total ambiguity implementation
        num = self.gap_len
        den = self.length + self.next.length + self.gap_len
        self.ambiguity_with_next = num / den
        self.all_ambigous_base_len += num
        self.ambiguity_so_far = self.all_ambigous_base_len / den

    def merge(self):
        """
        Merge this ColRange with its next ColRange.
        Update position, length, ambiguity accumulation.
        """
        if self.next is None:
            return

        # Extend length by: bases + non bases + next bases
        self.length = self.length + self.gap_len + self.next.length

        # After merging, the ambiguity accumulated is the non-base count
        # (your model does not accumulate ambiguity recursively beyond this)
        self.gap_len = (self.next.gap_len
                        if self.next is not None else 0)

        # Link to next->next
        self.next = self.next.next

        # Recompute new relative ambiguity
        self.update_ambiguity()

    def __lt__(self, other):
        """
        Required by heapq — it uses < for ordering.
        Compare by relative ambiguity.
        """
        return self.ambiguity_with_next < other.ambiguity_with_next


def build_colranges(base_regions, non_base_regions):
    """
    Convert (start,length) lists into linked ColRange objects.
    Connect each base region to the next base region.
    Attach the non-base length between them.
    """
    colranges = []

    # Create a list of all sequences represented by ColRange objects
    for (start, length) in base_regions:
        colranges.append(ColRange(start, length))

    # Now link them
    for i in range(len(colranges) - 1):
        A = colranges[i]
        B = colranges[i + 1]

        # Lookup the non-base region directly between A and B
        # TODO There should be better way to do this instead of looking through ALL the non base regions
        # Considering the non bases and bases are in order
        for nb_start, nb_len in non_base_regions:
            # We check to make sure we have the gap in between consecutive regions
            if nb_start == A.position + A.length and nb_start + nb_len == B.position:
                # We link sequence A with sequence B and compute the relative ambiguity
                A.next = B
                A.gap_len = nb_len
                A.update_ambiguity()
                break

    return colranges


def kmer_frequency(tuplist, step):
    """
    Generates a frequency table for how many probes ways we can cut with k length
    """
    hm = defaultdict(int)
    for kmerPos, kmerSize in tuplist:
        for k in range(50, 200, step):
            hm[k] += kmerSize // k
    return hm

def overlapping_kmer_frequency(tuplist, step):
    """
        Generates a frequency table dictionary for how many kmers we can cut with k length using sliding window
    """
    hm = defaultdict(int)
    for kmerPos, kmerSize in tuplist:
        for k in range(50, 200, step):
            if kmerSize >= k:
                hm[k] += kmerSize - k + 1
    return hm






def split_kmers(tupList, newKmerSize, orf):
    """
    Cuts the tup list in to individual k-mers of length: newKmerSize
    """
    splitKmers = []
    ctr = 1
    for start, length in tupList:
        while length - newKmerSize > 0:
            kmerName = f'{orf}_{newKmerSize}mers_{ctr}'
            splitKmers.append((kmerName, start, newKmerSize))
            length -= newKmerSize
            start += newKmerSize
            ctr += 1

    return splitKmers

def overlap_kmers(tupList, newKmerSize, orf):
    """Uses sliding window to get all k-mers """
    overlapKmers = []
    ctr = 1
    for start, length in tupList:
        while length - newKmerSize >= 0:
            kmerName = f'{orf}_{newKmerSize}mers_{ctr}'
            overlapKmers.append((kmerName, start, newKmerSize))
            length -= 1
            start += 1
            ctr += 1

    return overlapKmers

def heap_merge_kmers(base_regions, non_base_regions, max_ambig_fraction):
    """
    Bottom-up merging using a min-heap ordered by relative ambiguity.
    """
    # Even though the sorting is done with relative ambiguity we want to keep track of the total ambiguity
    # Using the TOTAL ambiguity is what we compare to P
    # Del sequence completo verificamos cuantos non bases (ambigous bases) hemos incluido

    # Build linked ColRanges
    colranges = build_colranges(base_regions, non_base_regions)

    # Map: position → object, for removing children later
    head_set = set(colranges)

    # Min heap
    heap = [cr for cr in colranges if cr.next]
    heapq.heapify(heap)

    # Perform merges
    while heap:
        cr = heapq.heappop(heap)

        # Skip if next vanished due to earlier merge
        if cr.next is None:
            continue

        # If exceeding ambiguity: stop merging this chain
        # Since we are using a min heap if the relative ambiguity exceeds our
        # max ambiguity fraction then all the col ranges after also exceed the
        # max ambiguity
        if cr.ambiguity_so_far > max_ambig_fraction:
            break

        # We will merge cr and cr.next:
        child = cr.next

        # Remove child from head set
        if child in head_set:
            head_set.remove(child)

        # Merge operation
        cr.merge()

        # If cr still has a next, reinsert it
        if cr.next is not None:
            heapq.heappush(heap, cr)

    # Now head_set contains ONLY top-level merged blocks
    final = []
    for head in head_set:
        final.append((head.position, head.length))

    final.sort(key=lambda x: x[0])
    return final