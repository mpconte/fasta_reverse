"""
    The program will take in two arguments: an input file and an output folder.

    Using python3 with no loops like for and while (use map filter reduce):

    1) Read an input FASTA file of nucleotide sequences (using built-in libraries)
    2) Reverse complement each sequence and output to another fasta file in the output directory (using built-in libraries)
    3) Plot the nucleotide composition of each reverse complement sequence (so a Trellis plot where each facet represents a sequence) and save as an html to the output directory (using Altair)
    4) Report a statistical test between the first two reverse complement sequences to see if they have significantly different distributions of nucleotides and print the results to stderr (using SciPy)
"""

import os
from sys import stderr, argv
from scipy import stats


NUCLEOTIDE_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

# sequence dictionary mapping sequence number to number of each nucleotide
sequences = {}

# global variable to keep track of current sequence number
current_sequence = -1


def reverse_nucleotide(line):
    """
    Perform reverse complement of nucleotide sequence and adds
    each sequence number to a dict

    @param line: nucletide sequence as a string
    """
    global current_sequence
    if line[0] == ">":
        current_sequence = int(line.strip()[1:])
        sequences[current_sequence] = {
            'A': 0, 'T': 0, 'G': 0, 'C': 0
        }
        return line
    try:
        sequence = ''.join(map(lambda x: NUCLEOTIDE_DICT[x], list(line.strip()[::-1])))
        sequences[current_sequence]['A'] += sequence.count('A')
        sequences[current_sequence]['T'] += sequence.count('T')
        sequences[current_sequence]['G'] += sequence.count('G')
        sequences[current_sequence]['C'] += sequence.count('C')
        return sequence + '\n'
    except KeyError:
        print("Invalid nucleotide")
        exit(1)


def reverse_asta_file(input_file, output_dir):
    """
    Perform reverse complement of given fasta file and write results to
    to given output directort

    @param input_file: input fasta file
    @param output_dir: output directory to write reverse complement of fasta file

    """
    # Open fasta file for reading
    with open(input_file, 'r') as fasta_input:
        # Reverse complement each sequence 
        fasta_lines = fasta_input.readlines()        
        reversed_fasta = map(reverse_nucleotide, fasta_lines)

        # Write result to file in the output directory
        with open(os.path.join(output_dir, "output.fasta"), 'w') as output:
            output.writelines(reversed_fasta)


if __name__ == "__main__":
    # Reverse complement fasta file as an argument and write result to the current directory as a second argument
    reverse_asta_file(argv[1], argv[2])

    # Calculate and print T test of nucleotide distributions of first two sequences
    sequence_nums = sorted(sequences.keys)
    first_seq = sequence_nums[0]
    sec_seq = sequence_nums[1]
    t_check = stats.ttest_ind([sequences[first_seq]['A'], sequences[first_seq]['T'], sequences[first_seq]['G'], sequences[first_seq]['C']],
                              [sequences[sec_seq]['A'], sequences[sec_seq]['T'], sequences[sec_seq]['G'], sequences[sec_seq]['C']])
    print('Nucleotide T test: {}'.format(t_check), file=stderr)
