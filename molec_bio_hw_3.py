######################################################################################
#
#   File:       molec_bio_hw_3.py
#
#   Purpose:    The
#
#   Developers:  Christopher Hitte
#               CSCI 5314, Spring 2017, HW3
#               William Diment
#               CSCI 4314, Spring 2017, HW3    
#
######################################################################################
#
#   Sample command line arguments to run program:
#
#   usage: molec_bio_hw_3.py -f FILE 
#
#   optional arguments:
#       -F FILE, -f FILE, --file FILE           fasta filename
#       -W WIDTH                    width of an output
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

WIDTH_CONST = 0

sequence_list = []
distance_matrix = None


class Sequence:
    def __init__(self, id, name):
        self.id = id
        self.name = name
        self.seq_data = ""

#we handle all results from the alignment in this class - we put in the sequences associated with the alignemnts,
#the alignments themselves, and then the length/distance associated with the alignment so that we can easily access this information across
#functions
class GlobalAlignmentResult:
    def __init__(self, seq_1, seq_2, distance, length, align_1, align_2):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        self.distance = distance
        self.length = length
        self.align_1 = align_1
        self.align_2 = align_2


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
    #we parse the file
    parse_file(args.file)
    #we calculate the alignments and build the distance matrix
    compute_global_distances()
    #we print the distance matrix!
    print_distance_matrix()


#
#
#
def parse_file(file):
    read_sequences(file)
    #we read the sequences from the file and initialize a global distance matrix that we will use 
    #the matrix is initialized to the size of the sequence lists, as we will need a N x M matrix to account for
    #each of the sequences
    global distance_matrix
    distance_matrix = numpy.zeros((len(sequence_list), len(sequence_list)))


#
#
#
def read_sequences(file):
    id = 0

    #we read the lines from the file
    line = file.readline()
    while line != '' and line[0] == '>':
        #we get the title of the sequence here and append it to a global sequence list
        new_sequence = Sequence(id, line[1:-1])
        sequence_list.append(new_sequence)

        #we look for every character between > and the end of the line so that we get every base in there
        #this is then passed to the class that we have declared above for ease of access
        seq_data = ""
        line = file.readline()
        while line != '' and line[0] != '>':
            seq_data += line[:-1]
            line = file.readline()

        new_sequence.seq_data = seq_data
        #we ID each sequence to associate it with a given number
        id += 1


#
#
#here is our runner function for computing both the alignments and the global distance matrix
def compute_global_distances():
    i = 0
    j = 0

    while i < len(sequence_list):
        seq_1 = sequence_list[i]

        while j < len(sequence_list):

            if j != i:
                seq_2 = sequence_list[j]
                #we return the matrix calculated below, having associated the sequence from the sequence list with the ID in the class list to get the sequence
                matrix = compute_global_distance_matrix(seq_1, seq_2)

                show_alignment = False
                #as an alignment pair is calculated it gets printed right away, the permission bit is flipped here
                if i % 2 == 0 and j - i == 1:
                    show_alignment = True
                #we then calculate the alignment from the matrix above
                calculate_global_distance(seq_1, seq_2, matrix, show_alignment)

            j += 1

        i += 1
        j = i


#
#
#
def compute_global_distance_matrix(seq_1, seq_2):
    #here we begin building our scoring matrix to calculate the alignment, initialize the data from the sequences
    #that we have in the sequence class
    seq_data_1 = ' ' + seq_1.seq_data
    seq_data_2 = ' ' + seq_2.seq_data
    #we initialize a matrix to the size of the sequences
    matrix = numpy.zeros((len(seq_data_2), len(seq_data_1)))

    #for the top row and column, we initialize their values according to the standard scoring system for glbobal alignment
    i = 0
    while i < len(seq_data_1):
        matrix[0][i] = i
        i += 1

    i = 0
    while i < len(seq_data_2):
        matrix[i][0] = i
        i += 1
   
   #we begin to iterate through the matrix to calculate our alignment. as is standard, we try to look for diagonal matches if possible
    i = 1
    j = 1
    while i < len(seq_data_2):
        while j < len(seq_data_1):

            # calc diagonal
            #if characters are aligned, we say its a match and give it a score of 0 so we can minimize it
            #if characters arent aligned we say its a mismatch so that we increase the distance
            diagonal = matrix[i - 1][j - 1]
            if seq_data_1[j] != seq_data_2[i]:
                diagonal += 1

            #if its a match we will have a score of the diagonal plus 0 - again so we minimize it
            #if its a mismatch, it will be score of the diagonal plus 1 - increasing the distance
            # calc top
            #same concept for the top, only we specifically use indel here
            top = matrix[i - 1][j] + 1
      
            # calc left
            #same concept for the left, again only using indel here
            left = matrix[i][j - 1] + 1

             #we calculate the minmimum to minmize the distance and then put that value in for the matrix
            matrix[i][j] = min(top, left, diagonal)

            j += 1

        i += 1
        j = 1

    return matrix


def calculate_global_distance(seq_1, seq_2, matrix, show_alignment=False):
    #again, we associate the sequences with the sequences in the Sequence class
    seq_data_1 = seq_1.seq_data
    seq_data_2 = seq_2.seq_data
    #we initialize our alignment strings
    align_1 = ""
    align_2 = ""
    #we get the length of the sequences here to properly iterate through the matrix. we will use a backtrace through the matrix
    i = len(seq_data_2) - 1
    j = len(seq_data_1) - 1
    while i > 0 and j > 0:
        #we take the values from three directions, above, left, diagonal. We are looking for the minimum here, and will choose from these values below
        u = matrix[i - 1][j]
        l = matrix[i][j - 1]
        d = matrix[i - 1][j - 1]
        #if the diagonal is at least less than or equal to both the left and top values, we choose the diagonal for our value as it represents a match between the two sequences
        if d <= l and d <= u:
            align_1 = seq_data_1[j] + align_1
            align_2 = seq_data_2[i] + align_2
            i -= 1
            j -= 1
        #if the left side is greater than the diagonal but less then the top, we choose the left side for our value, and it represents a gap in the second sequence
        elif l <= u:
            align_1 = seq_data_1[j] + align_1
            align_2 = '-' + align_2
            j -= 1
        #if the top is smaller then the left and the right, we choose the top and it also represents a gap in the first sequence
        else:
            align_1 = '-' + align_1
            align_2 = seq_data_2[i] + align_2
            i -= 1
    #if we reach the end of the first sequence we append the end of the sequence to the alignment here
    while j >= 0:
        align_1 = seq_data_1[j] + align_1
        j -= 1
    #if we reach the end of the second sequence we append the end of the sequence to the alignment here
    while i >= 0:
        align_2 = seq_data_2[i] + align_2
        i -= 1

    #we compare the lengths of the alignments here - if one is longer than the other, we correct it so that it will align properly
    if len(align_1) > len(align_2):
        align_2 = ((len(align_1) - len(align_2)) * '-') + align_2
    elif len(align_2) > len(align_1):
        align_1 = ((len(align_2) - len(align_1)) * '-') + align_1
    
    distance = 0
    length = len(align_1)

    #if the two sequences do not align, there is either a mismatch or an indel, so we increase the 'distance' by 1
    i = 0
    while i < length:
        if align_1[i] != align_2[i]:
            distance += 1

        i += 1
    #we build the distance matrix here. the formula is the mismatches/length, portrayed here as distance/length
    #e compare each sequence to every other sequence, using only the top half of the matrix
    #we do not need to use the full matrix, as the top half/bottom half of the matrix are simply mirrored across the diagonal
    distance_matrix[seq_1.id][seq_2.id] = round(float(distance)/length, 5)
    distance_matrix[seq_2.id][seq_1.id] = round(float(distance) / length, 5)

    if show_alignment:
        print_alignment(GlobalAlignmentResult(seq_1, seq_2, distance, length, align_1, align_2))


def print_alignment(align_result):
    #here we print the sequence names used in the alignments
    print "{0} vs. {1}".format(align_result.seq_1.name, align_result.seq_2.name)

    #here are the various quantifiers associated with each alignment
    distance = align_result.distance
    length = align_result.length
    align_1 = align_result.align_1
    align_2 = align_result.align_2

    index_str = ""
    #we go through here and signify where there is base in the two alignments that matches up, disregarding gaps in the alignment
    #we signify that it matches with a '|' character as is standard - we stick it into an 'index string' to match up with the alignment strings themselves
    i = 0
    while i < length:
        if align_1[i] == align_2[i] and align_1[i] != '-':
            index_str += '|'
        else:
            index_str += ' '

        i += 1

    if length > WIDTH_CONST:
        print WIDTH_CONST * '='
    else:
        print length * '='

    #here we build up the length of the of the print sequence that we want to use, according to the length of the alignment and the width constant we created at the start
    #the width constant makes it look nicer!
    printed = 0
    while printed + WIDTH_CONST < length:
        print align_1[printed : printed + WIDTH_CONST]
        print index_str[printed : printed + WIDTH_CONST]
        print align_2[printed : printed + WIDTH_CONST]
        print WIDTH_CONST * '='

        printed += WIDTH_CONST

    #here we print the alignments, along with the index string showing where there are matches
    #we then print the number of mismatches, along with 
    print align_1[printed : length]
    print index_str[printed : length]
    print align_2[printed : length]
    print "{0}({1}/{2})".format((length - printed) * '=', distance, length)
    print


#here is where we actually print the distance matrix as the function states
def print_distance_matrix():
    print "Distance Matrix for all Sequences:"
    print

    line = ""
    #we build up each lnie in the distance matrix here, and format it according to the specifications that it be within a 0-1 decimal range
    for row in distance_matrix:
        for column in row:
            line += "{0:.5f} ".format(round(column, 5))
        #we print the line, and then clear it for the nextl ine
        print line
        line = ""


if __name__ == "__main__":
    main()
