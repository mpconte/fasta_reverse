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
from functools import reduce
import itertools
from scipy import stats
import pandas as pd
import altair as alt
import numpy as np


# Dictionary mapping nucleotides to their complement
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
        sequences[current_sequence] = dict(zip(NUCLEOTIDE_DICT.keys(), [0] * len(NUCLEOTIDE_DICT.keys())))
        return line
    try:
        sequence = ''.join(map(lambda x: NUCLEOTIDE_DICT[x], list(line.strip()[::-1])))
        keys = NUCLEOTIDE_DICT.keys()
        sequences[current_sequence] = dict(zip(keys,
                                               list(map(lambda nucleotide: sequences[current_sequence][nucleotide] + sequence.count(nucleotide),
                                                        keys))))
        return sequence + '\n'
    except KeyError:
        print("Invalid nucleotide")
        exit(1)


def reverse_asta_file(input_file, output_dir):
    """
    Perform reverse complement of given fasta file and write results to
    to given output directory

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
    # Reverse complement given fasta file and write result to a given directory
    reverse_asta_file(argv[1], argv[2])

    # Plot the nucleotide composition of each reverse complement sequence
    sequence_nums = sorted(sequences.keys())
    nucleotides = list(NUCLEOTIDE_DICT.keys())
    frequencies = list(map(lambda sequence_num:
                           list(map(lambda nucleotide:
                                    sequences[sequence_num][nucleotide] /
                                    reduce(lambda a, b: a+b, sequences[sequence_num].values()) * 100.,
                                    nucleotides)),
                           sequence_nums))
    data = pd.DataFrame({
        'sequence': list(itertools.chain(*list(map(lambda sequence_num: np.repeat(sequence_num, 4), sequence_nums)))),
        'nucleotide': nucleotides * len(sequence_nums),
        'frequency': list(itertools.chain(*frequencies))
    })
    chart = alt.Chart(data).mark_point().encode(
        x='nucleotide:N',
        y='frequency:Q',
        row='sequence:Q'
    )
    chart.save("chart.html")

    # Calculate and print T test of nucleotide distributions of first two sequences
    first_seq = sequence_nums[0]
    sec_seq = sequence_nums[1]
    t_check = stats.ttest_ind(list(map(lambda nucleotide: sequences[first_seq][nucleotide], nucleotides)),
                              list(map(lambda nucleotide: sequences[sec_seq][nucleotide], nucleotides)))
    print('Nucleotide T test: {}'.format(t_check), file=stderr)
