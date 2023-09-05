import pandas as pd
import numpy as np
import os
import sys
import argparse
import collections
from itertools import groupby
from pprint import pprint


f = '2014_cladeC_HXBR_AY2_TRANSLATED.fasta'


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    with open(fasta_name, 'r') as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header_str = header.__next__()[1:].strip()
            # join all sequence lines into one.
            seq = "".join(s.strip() for s in faiter.__next__())
            yield header_str, seq


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for i, (k, v) in enumerate(my_gen):
        # resolve for duplicate sequence headers
        new_key = k.replace(" ", "_") + str(i).zfill(4)
        dct[new_key] = v.upper()

    return dct


def gethxb2(dict):
    """
    :param dict: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
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


def calc_alignedpos(ref_seq):
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


def posnumcalc(hxb2seq, start):

    pos_num = []
    n = start
    s = 0.01
    m = len(hxb2seq) - 1
    for i, resi in enumerate(hxb2seq):
        if i == 0 and resi == '-':
            print("Can't start numbering. HXB2 starts with a gap. Use a longer HXB2 sequence for the numbering")
        if i != m:
            if resi != '-' and hxb2seq[i+1] != '-':
                pos_num.append(n)
                n += 1
            elif resi != '-' and hxb2seq[i+1] == '-':
                g = n
                pos_num.append(g)
            elif resi == '-' and hxb2seq[i+1] == '-':
                g = n + s
                pos_num.append(g)
                s += 0.01
            elif resi == '-' and hxb2seq[i+1] != '-':
                g = n + s
                pos_num.append(g)
                s = 0.01
                n += 1
        else:
            if resi != '-':
                pos_num.append(n)
            elif resi == '-':
                g = n + s
                pos_num.append(g)
    return pos_num


def calc_new_offset(ref_pos, aligned_pos, offset):

    startpos = offset[0]
    endpos = offset[1]
    df = pd.DataFrame(columns=['ref_positions', 'aligned_positions'])
    df['ref_positions'] = ref_pos
    df['aligned_positions'] = aligned_pos
    start = df.loc[df['ref_positions'] == startpos, 'aligned_positions'].iloc[0]
    end = df.loc[df['ref_positions'] == endpos, 'aligned_positions'].iloc[0]
    new_pos = (start, end)

    return new_pos


ref_seq = gethxb2(fasta_to_dct(f))
pprint(ref_seq)
aln_pos = calc_alignedpos(ref_seq[1])


# pprint(aln_pos[1])
# print(aln_pos)
# df_refpos = calc_new_offset(aln_pos[0], aln_pos[1])
# print(df_refpos)


