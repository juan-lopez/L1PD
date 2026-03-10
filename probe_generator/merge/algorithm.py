# class structure wip
from collections import defaultdict
from functools import reduce
import heapq


class ColRange:
    """
    Represents a contiguous region of bases in the consensus.

    start           = start index in consensus
    length             = number of bases
    next               = next ColRange (the next contiguous base region)
    ambiguity_with_next = accumulated ambiguity % if merged with next
    gap_len            = number of ambigous bases X/N between this and next
    ambiguity_so_far    = ambiguity of the entire sequence including all merges performed so far
    all_ambigous_base_len    = length of all the ambigous bases we have merged in to this sequence
    """


    def __init__(self, start, length):
        self.start = start
        self.length = length
        self.ambiguity_with_next = float("inf")  # default until computed
        self.gap_len = 0
        self.prev = None
        self.next = None
        self.ambiguity_so_far = 0  # default until computed
        self.all_ambigous_base_len = 0


    def update_ambiguity_so_far(self):

        # Add up all the ambigous bases of the sequences being merged
        self.all_ambigous_base_len += self.gap_len + self.next.all_ambigous_base_len
        denominator = self.length + self.next.length + self.gap_len
        self.ambiguity_so_far = self.all_ambigous_base_len / denominator


    def update_ambiguity_with_next(self):
        """
        Computes ambiguity if this block merges with its next block.
        """
        if not self.next:
            return

        # DISCUSS
        # total ambiguity implementation
        # TODO DISCUSS this should only be theoretically updated
        # self.all_ambigous_base_len += self.next.all_ambigous_base_len + self.gap_len
        all_ambigous_with_next = self.all_ambigous_base_len + self.next.all_ambigous_base_len + self.gap_len
        # Must include gap len in all_ambigous_base_len
        denominator = self.length + self.next.length + self.gap_len

        # Use self.gap here + self.all_ambig_len and dont use it in the merge call
        # The formula should be  self.all_ambigous_base_len (left) + gap_len + self.all_ambigous_base_len (right)
        # self.ambiguity_with_next = self.all_ambigous_base_len  / denominator
        self.ambiguity_with_next = all_ambigous_with_next / denominator



    def merge(self):
        """
        Merge this ColRange with its next ColRange.
        Update start, length, ambiguity accumulation.
        """
        if not self.next:
            return

        # print(f'Merging sequence starting at {self.start} with length {self.length} with the sequence sarting at'
        #       f' {self.next.start} with length {self.next.length} ')
        # if self.prev:
        #     print(
        #         f'ambiguity so far for previous {self.prev.start} (including all previous merges) before the merge is {self.prev.ambiguity_so_far}  and ambiguity with next {self.ambiguity_with_next}')
        #
        # print(
        #     f'ambiguity so far for {self.start} (including all previous merges) before the merge is {self.ambiguity_so_far}  and ambiguity with next {self.ambiguity_with_next}')
        # print(
        #     f'ambiguity so far for {self.next.start} (including all previous merges) before the merge is {self.next.ambiguity_so_far} and ambiguity with next {self.next.ambiguity_with_next}')

        self.update_ambiguity_so_far()
        # print(
        #     f'ambiguity for {self.start} (including all previous merges) after the merge is {self.ambiguity_so_far} ')

        # Extend length by: bases + non bases + next bases
        self.length = self.length + self.gap_len + self.next.length


        # Update the ambigous bases to include the gap
        # self.all_ambigous_base_len += self.gap_len

        # After merging, the ambiguity accumulated is the non-base count
        # (your model does not accumulate ambiguity recursively beyond this)
        self.gap_len = self.next.gap_len

        if self.next.next:
            # Link to next->next
            self.next = self.next.next
        else:
            self.next = None

        # Recompute new relative ambiguity
        self.update_ambiguity_with_next()

        if self.prev:
            self.prev.update_ambiguity_with_next()


    def __lt__(self, other):
        """
        Required by heapq — it uses < for ordering.
        Compare by relative ambiguity.
        """
        return (self.ambiguity_with_next < other.ambiguity_with_next) or (self.ambiguity_with_next == other.ambiguity_with_next) and (self.length > other.length)


def build_colranges(base_regions):
    """
    Convert (start,length) lists into linked ColRange objects.
    Connect each base region to the next base region.
    Attach the non-base length between them.
    """
    if len(base_regions) < 2:
        return []

    col_ranges = []

    # Create a list of all sequences represented by ColRange objects
    for start, length in base_regions:
        col_ranges.append(ColRange(start, length))

    # Now link them
    for i in range(len(col_ranges) - 1):
        curr = col_ranges[i]
        next = col_ranges[i + 1]
        # Calculate initial ambiguity with next for current colRange
        curr.gap_len = next.start - (curr.start + curr.length)
        curr.next = next
        curr.update_ambiguity_with_next()
        curr.prev = col_ranges[i-1] if i-1 >= 0 else None

    col_ranges[-1].prev = col_ranges[-2]
    return col_ranges


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


def heap_merge_kmers(base_regions, max_ambig_fraction):
    """
    Bottom-up merging using a min-heap ordered by relative ambiguity.
    """
    # Even though the sorting is done with relative ambiguity we want to keep track of the total ambiguity
    # Using the TOTAL ambiguity is what we compare to P
    # Del sequence completo verificamos cuantos non bases (ambigous bases) hemos incluido

    # Build linked ColRanges
    colranges = build_colranges(base_regions)

    # Map: start → object, for removing children later
    head_set = set(colranges)

    # Min heap
    heap = [cr for cr in colranges if cr.next]
    heapq.heapify(heap)

    # Perform merges
    while heap:
        # print("Popping a value from the heap")
        cr = heapq.heappop(heap)

        # Skip if next vanished due to earlier merge
        if cr.next is None:
            continue

        # If exceeding ambiguity: stop merging this chain
        # Since we are using a min heap if the relative ambiguity exceeds our
        # max ambiguity fraction then all the col ranges after also exceed the
        # max ambiguity
        if cr.ambiguity_with_next > max_ambig_fraction:
            # print(f"Exiting due to ambiguity with next exceding max ambig fraction {cr.ambiguity_with_next} > {max_ambig_fraction}")
            break

        # We will merge cr and cr.next:
        child = cr.next

        # Remove child from head set
        if child in head_set:
            head_set.remove(child)

        # TODO we should also remove the child from the heap as it should onyl exist within the context of
        # its parent once it has bee merged
        if child in heap:
            # print("Removed the child base sequence as it has been merged")
            heap.remove(child)
            heapq.heapify(heap)

        # Merge operation
        cr.merge()

        # print(f'Merged {cr.start} and {child.start}')

        # TODO Instead of doing this I believe we should update the current colrange to point towards t
        # If cr still has a next, reinsert it
        if cr.next is not None:
            heapq.heappush(heap, cr)

    # Now head_set contains ONLY top-level merged blocks
    final = []
    for head in head_set:
        final.append((head.start, head.length))

    final.sort(key=lambda x: x[0])
    # print(final)
    res = 0
    for col_range in head_set:
        res += col_range.all_ambigous_base_len
    # print(res)
    #print(reduce(lambda x, y: x.all_ambigous_base_len + y.all_ambigous_base_len, list(head_set)))
    return final


# import sys
# base_regions = [tuple(map(int, base_region.split(','))) for base_region in sys.stdin.readline().strip().split('  ')]
# for col_range in build_colranges(base_regions):
#     print(col_range.start)
#     print(col_range.gap_len)
#     print(col_range.ambiguity_with_next)
#     print()
# print(heap_merge_kmers(base_regions, .02))
