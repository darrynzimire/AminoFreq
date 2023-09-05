import csv
from itertools import groupby
from pprint import pprint
import argparse
import sys
import pathlib
import collections
import subprocess
# external libraries
import pandas as pd
import os


f = '2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED.fasta'
soi = 'sites2.csv'
f2 = 'fasta_testoutfile.fasta'
def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


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

sites_of_interest_df = pd.read_csv(soi, sep=',', engine='python')
# print(sites_of_interest_df)
sites_list = sorted((sites_of_interest_df["sites"].dropna().unique().tolist()))
print(sites_list)
#
x = fasta_to_dct(f)
hbxref = gethxb2(x)
print(hbxref[1])
# pprint(hbxref[1])

pos_num = posnumcalc(hbxref[1], 1)
print(pos_num)

#dictionary with seq header as key:sequence as value
virus_d = fasta_to_dct(f)
# # pprint(virus_d)
#
# # for each sequence, collect the required sites
collected_sites_d = collections.defaultdict(list)
# print(collected_sites_d)
for seq_name, seq in virus_d.items():
    # print(seq_name, seq)
    for site in sites_list:
        idx = pos_num.index(site)
        # print(idx)
        p = collected_sites_d[seq_name].append(seq[idx])
        # print(p)

        with open(f2, 'w') as fh:
            for collected_name, collected_sites_list in collected_sites_d.items():
                sites_to_use = "".join(collected_sites_list)
                fh.write(f">{collected_name}\n{sites_to_use}\n")
sites = ','.join([str(int(x)) for x in sites_list])
# print(sites)

# pprint(x)




# def create_csv_with_sites(filename, sites):
#     with open(filename, 'w', newline='') as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow(['sites'])  # Write the column header
#
#         for site in sites:
#             writer.writerow([site])  # Write each site as a separate row
#
#     print(f"CSV file '{filename}' created successfully.")
#
# # Example usage
# sites = [332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367]
# filename = 'sites.csv'
# create_csv_with_sites(filename, sites)