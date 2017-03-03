######################################################################################
#
#   File:       chrishitte_hw2.py
#
#   Purpose:    The
#
#   Developer:  Christopher Hitte
#               CSCI 5314, Spring 2017, HW1
#
######################################################################################
#
#   Sample command line arguments to run program:
#
#   usage: chrishitte_hw2.py [-h] -F FILE -L [3-8]
#
#   optional arguments:
#       -h, --help           show this help message and exit
#       -F FILE, -f FILE, --file FILE           fasta filename
#       -L [3-8], -l [3-8], --length [3-8]      size of k-mer
#
######################################################################################
#
#   Time Complexity:
#
#   This algorithm runs in O((S * (Ls-K+1)) + N) time, where S is the number of
#   sequences in the FASTA file, Ls is the length of a given sequence, K is the given
#   k-mer length, and N is the number of unique k-mers in the FASTA file.
#
#   Memory Complexity:
#
#   This algorithm uses 4 data structures to calculate the result:
#       kmer_index          - dictionary that maps k-mer strings to KmerInstances
#       kmer_count_index    - dictionary that maps occurence counts to KmerInstances
#       kmer_seq_index      - dictionary that maps sequence counts to KmerInstances
#       KmerInstance        - class used to represent a k-mer once it is detected in
#                             a FASTA file.
#
#   This will result in Nk * KmerInstances being created where Nk is the number of
#   unique k-mers in the FASTA file. There will also be Nk * key strings in
#   kmer_index. The kmer_count_index will contain Nc * key integers where Nc is the
#   number of unique occurence counts within the set of KmerInstances. The
#   kmer_seq_index will contain Ns * key integers where Ns is the number of unique
#   sequence counts within the set of KmerInstances.
#
######################################################################################
#
#   References: Okeson Jan 2016
#               alexokeson_hw1.py
#               Formatting of the header comment and function comments
#
######################################################################################

import argparse
import numpy
import scipy
from scipy import spatial
import time

WIDTH_CONST = 0

sequence_list = []


class Sequence:
    def __init__(self, name):
        self.name = name
        self.seq_data = ""

#
#
#
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '-f', '--file', type=file, help="fasta filename", required=True)
    parser.add_argument('-W', '-w', '--width', type=int, help="width of output", default=60)

    args = parser.parse_args()

    global WIDTH_CONST
    WIDTH_CONST = args.width

    parse_file(args.file)


#
#
#
def parse_file(file):
    read_sequences(file)

    compute_global_distances()


#
#
#
def read_sequences(file):
    line = file.readline()

    while line != '' and line[0] == '>':
        new_sequence = Sequence(line[1:-1])
        sequence_list.append(new_sequence)

        seq_data = ""
        line = file.readline()
        while line != '' and line[0] != '>':
            seq_data += line[:-1]
            line = file.readline()

        new_sequence.seq_data = seq_data


#
#
#
def compute_global_distances():
    i = 0
    j = 1

    while j < len(sequence_list):
        seq_1 = sequence_list[i]
        seq_2 = sequence_list[j]

        compute_global_distance(seq_1, seq_2)

        i += 2
        j += 2


#
#
#
def compute_global_distance(s_seq_1, s_seq_2):

    seq_1 = ' ' + s_seq_1.seq_data
    seq_2 = ' ' + s_seq_2.seq_data

    '''
    if len(seq_1) > len(seq_2):
        i = len(seq_1) - len(seq_2)
        seq_2 += (i * '-')
    elif len(seq_1) < len(seq_2):
        i = len(seq_2) - len(seq_1)
        seq_1 += (i * '-')
    '''

    matrix = numpy.zeros((len(seq_2), len(seq_1)))

    score = 0
    match = 0
    indel = 1
    mismatch = 1

    i = 0
    while i < len(seq_1):
        matrix[0][i] = i
        i += 1

    i = 0
    while i < len(seq_2):
        matrix[i][0] = i
        i += 1
   
    i = 1
    j = 1
    while i < len(seq_2):
        while j < len(seq_1):

            # calc diagonal
            #if characters are aligned, we say its a match and give it a score of 0 so we can minimize it
            #if characters arent aligned we say its a mismatch so that we increase the distance
            diagonal = matrix[i - 1][j - 1]
            if seq_1[j] != seq_2[i]:
                diagonal += 1
            elif seq_1[j] == seq_2[i]:
                k = 4

            #if its a match we will have a score of the diagonal plus 0 - again so we minimize it
            #if its a mismatch, it will be score of the diagonal plus 1 - increasing the distance
            # calc top
            #same concept for the top, only we specifically use indel here
            top = matrix[i - 1][j] + 1
      
            # calc left
            #same concept for the left, again only using indel here
            left = matrix[i][j - 1] + 1

            matrix[i][j] = min(top, left, diagonal)
            #like the matrix we had, we calculate the minmimum to minmize the distance and then put that value in for the matrix

            j += 1

        i += 1
        j = 1

    process_result(seq_1, seq_2, matrix)


def process_result(seq_1, seq_2, matrix):
    align_1 = ""
    align_2 = ""

    i = len(seq_2) - 1
    j = len(seq_1) - 1
    while i > 0 and j > 0:
        u = matrix[i - 1][j]
        l = matrix[i][j - 1]
        d = matrix[i - 1][j - 1]

        if d <= l and d <= u:
            align_1 = seq_1[j] + align_1
            align_2 = seq_2[i] + align_2
            i -= 1
            j -= 1
        elif l <= u:
            align_1 = seq_1[j] + align_1
            align_2 = '-' + align_2
            j -= 1
        else:
            align_1 = '-' + align_1
            align_2 = seq_2[i] + align_2
            i -= 1

    while j > 0:
        align_1 = seq_1[j] + align_1
        j -= 1

    while i > 0:
        align_2 = seq_2[i] + align_2
        i -= 1

    if len(align_1) > len(align_2):
        align_2 = ((len(align_1) - len(align_2)) * '-') + align_2
    elif len(align_2) > len(align_1):
        align_1 = ((len(align_2) - len(align_1)) * '-') + align_1

    
    print_result(align_1, align_2)





def print_result(align_1, align_2):
    length = len(align_1)

    index_str = ""
    score = 0

    i = 0
    while i < length:
        if align_1[i] == align_2[i] and align_1[i] != '-':
            score += 1
            index_str += '|'
        else:
            index_str += ' '

        i += 1

    if length > WIDTH_CONST:
        print WIDTH_CONST * '='
    else:
        print length * '='

    printed = 0
    while printed + WIDTH_CONST < length:
        print align_1[printed : printed + WIDTH_CONST]
        print index_str[printed : printed + WIDTH_CONST]
        print align_2[printed : printed + WIDTH_CONST]
        print WIDTH_CONST * '='

        printed += WIDTH_CONST

    print align_1[printed : length]
    print index_str[printed : length]
    print align_2[printed : length]
    print "{0}({1}/{2})".format((length - printed) * '=', score, length)
    distMatrix(align_1, align_2)



def distMatrix(align_1, align_2):
    i = 0
    j = 0

    #says here that we should use a 'normal scoring for global' for the distance matrix here - https://piazza.com/class/iynxmxbzwq35lz?cid=26

    score = 0
    match = 0
    mismatch = 1
    indel = 1
    #we create an array according to the lengths of the alignments
    scoreMatrix = numpy.zeros((len(align_1), len(align_2)))

    #we fill the array with the placeholder values
    while i < len(align_1):
        scoreMatrix[0][i] = i
        i += 1

    i = 0
    while i < len(align_2):
        scoreMatrix[i][0] = i
        i += 1

    i = 1
    j = 1
    while i < len(align_2):
        while j < len(align_1):
            #if there is a match, the score is equivalent to the match 
            if align_1[j] == align_2[i]:
                score = match
            else:
                #if its not there is a mismatch
                score = mismatch
            #like above, we need to take the minimum - we want the smallest possible distance
            #the hw3 file on moodle says something about have everything be a value of less than 1
            #but he has us using a different scoring metric in that piazza post
            #note: this is somewhat similar to the APSP problem 
            diagScore = scoreMatrix[i-1][j-1] + score
            topScore = scoreMatrix[i-1][j] + indel
            leftScore = scoreMatrix[i][j-1] + indel
            scoreMatrix[i][j] = min(diagScore, topScore, leftScore)
            j+=1
        i+=1
        j = 1
    print scoreMatrix

if __name__ == "__main__":
    main()
