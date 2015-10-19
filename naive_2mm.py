#!/usr/bin/env python
'''
Generic Naive Approximate Matching which allows n mismatches

Adapted and Enhanced by Jordan Cheah, Oct 2015
Acknowledgement - Original Author: Ben Langmead
'''

import string
  
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences
    
def naive_with_counts(p, t):
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(t) - len(p) + 1):
        match = True
        num_alignments += 1
        for j in range(len(p)):
            num_character_comparisons += 1
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences, num_alignments, num_character_comparisons
    
    
'''
Naive Approximate Matching up to 2 mismatches per occurence.
Reverse complement are not considered here.
'''
def naive_2mm(p, t):
    occurrences = []
    
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatch_count = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatch_count += 1
                if mismatch_count > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

'''
Naive Approximate Matching up to 2 mismatches per occurence + output algorithm performance (num_alignments, num_character_comparisons)
Reverse complement are not considered here.
'''
def naive_2mm_with_count(p, t):
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatch_count = 0
        num_alignments += 1
        for j in range(len(p)):  # loop over characters
            num_character_comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                mismatch_count += 1
                if mismatch_count > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record

    return occurrences, num_alignments, num_character_comparisons


'''
Naive Approximate Matching up to 2 mismatches per occurence + output algorithm performance (num_alignments, num_character_comparisons)
Reverse complement are not considered here.
'''
def naive_nmm_with_count(p, t, n):
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatch_count = 0
        num_alignments += 1
        for j in range(len(p)):  # loop over characters
            num_character_comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                mismatch_count += 1
                if mismatch_count > n:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record

    return occurrences, num_alignments, num_character_comparisons


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
 
        
if __name__ == '__main__':
        
    p = 'word'
    t = 'there would have been a time for such a word'
    occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
    print(occurrences, num_alignments, num_character_comparisons)  # Expect ([40], 41, 46)
    
    p = 'needle'
    t = 'needle need noodle needle'
    occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
    print(occurrences, num_alignments, num_character_comparisons)   # Expect ([0, 19], 20, 35)
    
    # Testing of naive_2mm function
    # Result expected: list [0, 4]
    result = naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT')  
    print(result)
    
    # Testing of generic naive_nmm function with up to n mismatches
    # Result expected: list [0, 4]
    result = naive_nmm_with_count('ACTTTA', 'ACTTACTTGATAAAGT', 2)   
    print(result)
    
    # Human Chromosome 1
    # Source: http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/chr1.GRCh38.excerpt.fasta
    t = readGenome('chr1.GRCh38.excerpt.fasta')    
    p = "GGCGCGGTGGCTCACGCCTGTAAT"
    occurrences, num_alignments, num_character_comparisons = naive_nmm_with_count(p, t, 2)
    print("*** RESULTS: naive_2mm_with_counts ***")
    print(occurrences, num_alignments, num_character_comparisons) 
    # Expect [56922, 84641, 147558, 160162, 160729, 191452, 262042, 273669, 364263, 421221, 429299, 465647, 
    #    551134, 635931, 657496, 681737, 717706, 724927, 747359]
    print(len(occurrences))  # Expect 19

    