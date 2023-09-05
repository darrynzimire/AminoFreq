import pandas as pd
import numpy as np
from src.app import *
from Bio import AlignIO, SeqIO
time_point_1 = ''
time_point_2 = ''


def calculate_freq(file):
    get_tuple = yield_tuple(file)
    seq_dict = fasta_to_dct(file)
    ref_name, ref_seq = gethxb2(seq_dict)
    calc_align_pos = calc_hxb2pos2(ref_name)
    read_alignment = AlignIO.read(file, 'fasta')
    frequencies = alnSiteCompositionDF(read_alignment)
    return frequencies


freq = calculate_freq('../data/V703_T1_AA.fasta')
print(freq)

# freq2 = calculate_freq('../data/V703_C3_AA.fasta')
# print(freq2)