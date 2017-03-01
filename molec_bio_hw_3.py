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

    args = parser.parse_args()

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

    while j <= len(sequence_list):
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
            diagonal = matrix[i - 1][j - 1]
            if seq_1[j] != seq_2[i]:
                diagonal += 1
            else:
                k = 4

            # calc top
            top = matrix[i - 1][j] + 1

            # calc left
            left = matrix[i][j - 1] + 1

            matrix[i][j] = min(top, left, diagonal)

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
    index_str = ""
    score = 0

    i = 0
    while i < len(align_1):
        if align_1[i] == align_2[i] and align_1[i] != '-':
            score += 1
            index_str += '|'
        else:
            index_str += ' '

        i += 1

    print align_1
    print index_str
    print align_2
    print "{0}({1}/{2})".format(len(align_1) * '=', score, len(align_1))
    print



if __name__ == "__main__":
    main()
