import os
import pandas as pd
import altair as alt
import vega_datasets as vds
from scipy import stats
import numpy as np
from sys import stderr, argv


NUCLEUS = {'A':'T','T':'A','G':'C','C':'G'}
sequences = {}
current_sequence = -1

def reverse_nucleotide(line):
    global current_sequence
    if line[0] == ">":
        current_sequence = int(line.strip()[1:])
        sequences[current_sequence] = {
            'A': 0, 'T': 0, 'G': 0, 'C': 0
        }
        return line
    try:
        sequence = ''.join(map(lambda x: NUCLEUS[x] if x in NUCLEUS.keys() else x, list(line.strip()[::-1])))            
        sequences[current_sequence]['A'] += sequence.count('A')
        sequences[current_sequence]['T'] += sequence.count('T')
        sequences[current_sequence]['G'] += sequence.count('G')
        sequences[current_sequence]['C'] += sequence.count('C')
        return sequence + '\n'
    except KeyError:
        print("Invalid nucleotide")
        exit(1)


def calc_distribution(sequence):
    distribution = {}
    total = sequence['A'] + sequence['T'] + sequence['G'] + sequence['C']
    distribution['A'] = sequence['A'] / total
    distribution['T'] = sequence['T'] / total
    distribution['G'] = sequence['G'] / total
    distribution['C'] = sequence['C'] / total
    return distribution

# The program will take in two arguments: an input file and an output folder.

# Using python3 with no loops like for and while (use map filter reduce):

# Read an input FASTA file of nucleotide sequences (using built-in libraries)
# Reverse complement each sequence and output to another fasta file in the output directory (using built-in libraries)
# Plot the nucleotide composition of each reverse complement sequence (so a Trellis plot where each facet represents a sequence) and save as an html to the output directory (using Altair)
# Report a statistical test between the first two reverse complement sequences to see if they have significantly different distributions of nucleotides and print the results to stderr (using SciPy)

def reverse_asta_file(input_file, output_dir): 
    # source = vds.data.cars()    
    # array = source.to_numpy()    
    # chart = alt.Chart(source).mark_point().encode(
    #     x='Horsepower:Q',
    #     y='Miles_per_Gallon:Q',
    #     row='Origin:N'
    # ).interactive()    
    # chart.show()

    # Open fasta file for reading
    with open(input_file, 'r') as fasta_input:
        # Reverse complement each sequence 
        fasta_lines = fasta_input.readlines()        
        reversed_fasta = map(reverse_nucleotide, fasta_lines)
        
        # Write result to file in the output directory
        with open(os.path.join(output_dir, "output.fasta"), 'w') as output:
            output.writelines(reversed_fasta)

        # # Plot trellis plot of nucleotide composition of each sequence
        
        # data = pd.DataFrame( {
        #     'sequence_num': sorted(sequences.keys()),
        #     'A': map(lambda sequence: sequence['A'], sequences.iteritems),
        #     'T': map(lambda sequence: sequence['T'], sequences), 
        #     'G': map(lambda sequence: sequence['G'], sequences),             
        #     'C': map(lambda sequence: sequence['C'], sequences) 
        # })
        # chart = alt.Chart(sequences)
        # chart.mark_point().encode(
        #     x='nucleotide:Q',
        #     y='count:Q',
        #     row='sequence_num:N'
        # )

        # Calculate and print distribution of nucleotides of first two sequences        
        t_check = stats.ttest_ind([sequences[1]['A'], sequences[1]['T'], sequences[1]['G'], sequences[1]['C']], 
                                    [sequences[2]['A'], sequences[2]['T'], sequences[2]['G'], sequences[2]['C']])
        # print('\tA\tT\tG\tC', file=stderr)
        # print('\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(fst_sequence_dist['A'], fst_sequence_dist['T'], fst_sequence_dist['G'], fst_sequence_dist['C']), file=stderr)
        # print('\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(sec_sequence_dist['A'], sec_sequence_dist['T'], sec_sequence_dist['G'], sec_sequence_dist['C']), file=stderr)
        print(t_check, file=stderr)

# Pass file to input fasta file as an argument
if __name__ == "__main__":
    reverse_asta_file(argv[1], os.path.curdir)