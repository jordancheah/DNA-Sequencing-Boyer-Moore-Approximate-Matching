#!/usr/bin/env python

"""
Two algorithms are presented here:
a) Boyer-Moore approximate matching with Pigeon Hole.
b) Index-assisted approximate matching with k-mer index, Pigeon Hole

Adapted and Enhanced by Jordan Cheah, Oct 2015
Acknowledgement - Original Author: Ben Langmead
"""

import string
import boyermoore
import bisect


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    

class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            #print(self.index[i][1])
            i += 1
        #print("------end of index hits--------")
        return hits

        
'''
Implement the pigeonhole principle using Index to find exact matches for the partitions. 
Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches 
(substitutions). We will use an 8-mer index. 

Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, finds all 
approximate occurrences of P within T with up to 2 mismatches. Insertions and deletions are not allowed. 
Don't consider any reverse complements.

How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence, 
occur with up to 2 substitutions in the excerpt of human chromosome 1? (Don't consider reverse complements 
here.)

'''

# Pigeon Hole / Boyer-Moore Approximate Matching
def approximate_match(p, t, n):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches)

# Pigeon Hole / K-Mer index Approximate Matching
def queryIndex_approximate_match(p, t, n, index):
    # 8-mer index, p is 24, and we allow up to two mismatches
    # Divide p into 3 segments, so at least one of the segment will be an exact match

    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    indexhits = 0
    for i in range(n+1):
        # split p into n+1 segments
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = index.query(p[start:end])
        indexhits += 1
        
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    print("indexhits = {}", indexhits)
    return list(all_matches)
        
if __name__ == '__main__':
    
    p = 'AACTTG'
    t = 'CACTTAATTTG'
    print(approximate_match(p, t, 2))       # Expect [0, 5]
    
    t = readGenome('chr1.GRCh38.excerpt.fasta')
    p = "GGCGCGGTGGCTCACGCCTGTAAT"

    print("*** RESULTS: approximate_match using Pigeon Hole / Boyer's Moore ***")
    print(approximate_match(p, t, 2))

    print("*** RESULTS: approximate_match using Pigeon Hole / K-Mer Index ***")
    index = Index(t,8)
    print(queryIndex_approximate_match(p, t, 2, index))


    
    