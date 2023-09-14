import pandas as pd
import numpy as np
import collections
from itertools import groupby
from Bio import AlignIO, SeqIO
from database import create_table_from_dataframe
import output
import plots
import argparse
import csv
import os


def yield_tuple(f):

    with open(f, 'r') as _f:
        faiter = (x[1] for x in groupby(_f, lambda line: line[0] == ">"))
        for header in faiter:
            header_str = header.__next__()[1:].strip()
            seq = ''.join(s.strip() for s in faiter.__next__())

            yield header_str, seq


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given.
    Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = yield_tuple(file_name)
    for i, (k, v) in enumerate(my_gen):
        # resolve for duplicate sequence headers
        new_key = k.replace(" ", "_") + str(i).zfill(4)
        dct[new_key] = v.upper()

    return dct


def extract_sites_to_list(file_path):
    sites_list = []

    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sites_list.append(int(row['sites']))

    return sites_list


def gethxb2(dict):
    """
    :param dict: a dictionary of your aligned input sequences.
    Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence as a string
    """
    found = False
    hxb2_seq = None
    hxb2_key = None
    for k in dict.keys():
        if "HXB2" in k.upper():
            found = True
            hxb2_key = k
            hxb2_seq = dict[k]
            print(f"Found hxb2 ref. seq. Its full name is: {hxb2_key}")
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")
    return str(hxb2_key), str(hxb2_seq)


def calculate_sequence_length(sequence):
    length = 0
    for char in sequence:
        if char != "-":
            length += 1
    return length


# Calculate the HXB2 aligned positions.
def calc_hxb2pos(hxb2seq, start):

    pos_num = []
    n = start
    s = 0.01

    for i in range(len(hxb2seq)):
        resi = hxb2seq[i]

        if i == 0 and resi == '-':
            raise ValueError(
                "Can't start numbering. HXB2 starts with a gap. "
                "Use a longer HXB2 sequence for numbering")

        if resi != '-':
            pos_num.append(n)
            if i < len(hxb2seq) - 1 and hxb2seq[i + 1] == '-':
                g = n
                pos_num.append(g)
        else:
            g = n + s
            pos_num.append(g)
            s = 0.01
        n += 1

    return pos_num


def calc_hxb2pos2(ref_seq):
    """
      Modified from http://github.com/smallBixtools
      :param ref_seq:
      :return:
      """
    aligned_columns = []
    ref_pos = list(range(1, len(ref_seq)))

    character_count = 0
    added_for_this_char = False
    for curnt_col, seq_char in enumerate(ref_seq):
        curnt_col = curnt_col + 1
        if curnt_col >= len(ref_seq):
            break
        if seq_char != '-':
            character_count += 1
            added_for_this_char = False

        if character_count in ref_pos:
            if not added_for_this_char:
                aligned_columns.append(curnt_col)
                added_for_this_char = True
    if len(aligned_columns) != len(ref_pos):
        while len(aligned_columns) < len(ref_pos):
            aligned_columns.append("NA")

    return ref_pos, aligned_columns


def parse_sequence_alignment(alignment_file):
    """
    Parse the sequence alignment file and extract amino acid sequences.
    """
    sequences = []
    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences


#calculate the amino acid frequencies at each position in the alignment
def alnSiteCompositionDF(aln, characters='ARNDCQEGHILKMFPSTWYV', hxb2_identifier='HXB2'):
    alnRows = aln.get_alignment_length()
    compDict = {char: [0] * alnRows for char in characters}

    for record in aln:
        header = record.id
        seq = record.seq

        # Check if the header contains the HXB2 identifier
        if hxb2_identifier in header.upper():
            continue  # Skip this record and move to the next one

        for aaPos in range(len(seq)):
            aa = seq[aaPos]
            if aa in characters:
                compDict[aa][aaPos] += 1

    return pd.DataFrame.from_dict(compDict)


def makerefcompdf(hxb2seq):

    seq = [aa for aa in hxb2seq if aa != "-"]
    s_seq = ''.join(seq)
    return s_seq


# Normalization of frequency data
def get_residue_frequency(df, position):
    position_data = df.iloc[position -1]
    residue = position_data.idxmax()
    frequency = position_data[residue] / position_data.sum()
    freq_percent = round(position_data[residue] / position_data.sum() *100, 2)
    return frequency, freq_percent


def find_residue_frequencies(clade_c_panel, aligned_pos, query_pos):
    residue_freq = {}

    # Iterate through Clade C reference panel sequences
    total_sequences = len(clade_c_panel)
    for sequence in clade_c_panel:
        residue = sequence.seq[aligned_pos]

        # Update the count for each residue
        residue_freq.setdefault(residue, 0)
        residue_freq[residue] += 1

    # Calculate the frequency of each residue
    for residue, count in residue_freq.items():
        frequency = round(count / total_sequences * 100, 4)
        # residue_freq[residue] = frequency
        residue_freq[residue] = count

    return (query_pos, aligned_pos), residue_freq


def get_alignment_position_and_residue(freq_table, hxb2_seq, target_position):

    alignment_position = 0
    counter = 0
    residue = ''
    for i in range(len(hxb2_seq)):
        if hxb2_seq[i] != '-':
            counter += 1

        if counter == target_position:
            alignment_position = i + 1 # Add 1 to convert from zero-indexed to one-indexed position
            residue = hxb2_seq[i]
            break
    residue_freq = get_residue_frequency(freq_table, alignment_position)
    # print(residue_freq)
    # print(residue_freq[1])
    # print(f"Query Position: {target_position}")
    # print(f"Aligned Position: {alignment_position}")
    # print(f"Reference Residue at position {target_position}: {residue}")
    # print(f"Amino acid frequency at position {target_position}: {residue_freq[0]}- {residue_freq[1]}%\n")
    # print('{} {} {}{}\n'.format(f"Amino acid frequency at position", target_position, residue_freq[1], '%'))
    return target_position, alignment_position, residue
    # return residue_freq, alignment_position


def convert_list_to_dataframe(data_list):

    # Create an empty DataFrame with columns representing each amino acid
    data_df = pd.DataFrame(columns=list(data_list[0][1].keys()))
    # Iterate over the list of tuples and fill in the frequencies
    for position, frequencies in data_list:
        data_df.loc[position] = [frequencies.get(aa, 0.0) for aa in data_df.columns]
    data_df = data_df.reset_index(drop=True)
    data_df = data_df.rename_axis('pos')
    # data_df = data_df.reset_index(drop=True)
    return data_df


def convert_to_frequency_and_mapping(data_structure):
    # Step 1: Determine unique amino acid symbols
    unique_amino_acids = sorted(set(amino_acid for _, amino_acid_counts in data_structure for amino_acid in amino_acid_counts.keys()))

    # Step 2: Create amino acid symbol-to-index mapping
    amino_acid_to_index = {amino_acid: index for index, amino_acid in enumerate(unique_amino_acids)}

    # Determine the number of positions and amino acids
    num_positions = len(data_structure)
    num_unique_amino_acids = len(unique_amino_acids)

    # Step 3: Initialize an empty numpy.ndarray
    frequency_array = np.zeros((num_positions, num_unique_amino_acids), dtype=int)

    # Step 4: Fill the array with frequency counts data
    for row, (position, amino_acid_counts) in enumerate(data_structure):
        for amino_acid, count in amino_acid_counts.items():
            column_index = amino_acid_to_index[amino_acid]
            frequency_array[row, column_index] = count

    return frequency_array, amino_acid_to_index


def countmat_from_dataframe(dataframe):
    # Convert the filtered DataFrame to a numpy.ndarray
    frequency_array = dataframe.to_numpy()

    # Get the column names as the alphabet
    alph = list(dataframe.columns)
    alphabet = ''.join(alph)
    return frequency_array, alphabet


def main(file_path, query_sites):

    # get the hxb2 reference sequence
    hxb2name, hxb2seq = gethxb2(fasta_to_dct(file_path))

    aln_pos2 = calc_hxb2pos2(hxb2seq)
    refpos = aln_pos2[0]
    aln_refpos = aln_pos2[1]
    alignment = AlignIO.read(file_path, 'fasta')
    df_freq = alnSiteCompositionDF(alignment)
    # print(df_freq)
    df_freq.to_csv("../deployment/Dash_app/Clade_Cref.csv", index=True, header=True)
    table_name = 'CladeC_refpanel'
    create_table_from_dataframe(df_freq, table_name)
    # row_df = df_freq.iloc[[446]]
    # print(row_df)
    sites = extract_sites_to_list(query_sites)
    df_data = []
    al = []

    for i in sites:
        aln_pos = get_alignment_position_and_residue(df_freq, hxb2seq, i)
        all_freq = find_residue_frequencies(alignment, aln_pos[1], i)
        al.append(aln_pos[1] -1)
        df_data.append(all_freq)
    filtered_df = df_freq.loc[al]
    print(filtered_df)
    raw_counts_csv = output.generate_outputdir('the_rawcounts.csv', zipped=False)
    filtered_df.to_csv(raw_counts_csv, header=True, index=False)
    count_mat, mapped_alphabet = countmat_from_dataframe(filtered_df)
    # output.generate_report(df_data, 'frequencies_for_sites.txt')
    plots.make_logogram(count_mat, mapped_alphabet, sites,  r'../data/weblogo_config.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script generates amino-acid reference frequencies from '
                                                 '200 reference Clade C panel viral')
    parser.add_argument('-f', type=str, metavar='', required=True, help='Input FASTA file containing clade c sequences')
    parser.add_argument('-i', type=str, metavar='', required=False, help='list of sites of interest in csv format file')
    parser.add_argument('-c', type=str, metavar='', required=False)
    args = parser.parse_args()
    current_directory = os.getcwd()
    file_path = args.f
    sites_of_interest = args.i
    # file_path = os.path.join(current_directory, r'../raw_data', args.f)
    # sites_of_interest = os.path.join(current_directory, r'..\raw_data', args.i)
    main(file_path, sites_of_interest)

    # Assuming you have already defined 'alignment' as your Bio.Align.MultipleSeqAlignment object
    # Save the DataFrame to the SQLite database

    # Read the DataFrame back from the SQLite database
    # retrieved_df = read_dataframe_from_table(table_name)

    # test_f = 'fasta_testoutfile.fasta'
    # logo_from_seq(test_f)
