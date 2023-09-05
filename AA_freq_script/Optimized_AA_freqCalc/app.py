import argparse
import os
import pathlib
import subprocess
from Bio import AlignIO
from itertools import groupby
import collections
import pandas as pd

__author__ = 'Colin Anthony'

# original_pos: len(HXB2_seq[np.array]
# HXB2_refpos: calc_ref[]
# query_pos: [np.array]
# aa_freq: [np.array]
#norm_aa_freq: [np.array]


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
        # resolve duplicate sequence headers
        new_key = f"{k.replace(' ', '_')}{i:04d}"
        dct[new_key] = v.upper()

    return dct



def get_hxb2(sequence_dict):
    """
    :param sequence_dict: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence as a string
    """
    for k, seq in sequence_dict.items():
        if "HXB2" in k.upper():
            print(f"Found HXB2 reference sequence: {k}")
            return k, seq

    print("We could not find a sequence with 'HXB2' in its name. "
          "Please make sure you have an HXB2 ref. seq. included")
    return None, None


def pos_num_calc(hxb2_seq, start):
    pos_num = []
    n = start
    s = 0.01
    m = len(hxb2_seq) - 1
    for i, resi in enumerate(hxb2_seq):
        if i == 0 and resi == '-':
            print("Can't start numbering. HXB2 starts with a gap. Use a longer HXB2 sequence for the numbering")
        if i != m:
            if resi != '-' and hxb2_seq[i + 1] != '-':
                pos_num.append(n)
                n += 1
            elif resi != '-' and hxb2_seq[i + 1] == '-':
                g = n
                pos_num.append(g)
            elif resi == '-' and hxb2_seq[i + 1] == '-':
                g = n + s
                pos_num.append(g)
                s += 0.01
            elif resi == '-' and hxb2_seq[i + 1] != '-':
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


def aln_site_composition_df(aln, characters="ACDEFGHIKLMNPQRSTVWY"):
  alnRows = aln.get_alignment_length()
  compDict = {char:[0]*alnRows for char in characters}
  for record in aln:
    header = record.id
    seq = record.seq
    for aaPos in range(len(seq)):
      aa = seq[aaPos]
      if aa in characters:
        compDict[aa][aaPos] += 1
  return pd.DataFrame.from_dict(compDict)


def main():

    f = '../2014_cladeC_HXBR_AY2_TRANSLATED.fasta'
    infile = AlignIO.read(f, "fasta")
    dic = fasta_to_dct(f)
    ref_dic = get_hxb2(dic)
    seq = str(ref_dic[0])
    print(type(ref_dic))
    hxb2_aln_pos = pos_num_calc(seq, 1)
    print(hxb2_aln_pos)

    df = aln_site_composition_df(infile, characters="ACDEFGHIKLMNPQRSTVWY")
    print(df)


if __name__ == '__main__':
    main()