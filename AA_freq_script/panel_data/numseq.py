from Bio import SeqIO

f = '2014_cladeCpanel_V703_synthseq_refCGHXBR.fasta'
f2 = '2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED.fasta'
f3 = '2014_cladeCpanel_V703_synthseq_HXBR_AY_CODON-ALIGNMENT.fasta'


def count_sequences(fasta_file):
    count = 0
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            count += 1
    return count

# fasta_file 3


sequence_count = count_sequences(f3)
print("Number of sequences:", sequence_count)
