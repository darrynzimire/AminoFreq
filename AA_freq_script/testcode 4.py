import pandas as pd
from Bio import SeqIO, AlignIO
# from Optimized_AA_freqCalc import visualization as vis


def parse_sequence_alignment(alignment_file):
    """
    Parse the sequence alignment file and extract amino acid sequences.
    """
    sequences = []
    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences


def calculate_amino_acid_frequencies(sequences):
    """
    Calculate the frequency of each amino acid at each position in the alignment.
    """
    frequencies = pd.DataFrame(0, columns=range(len(sequences[0])),
                               index=["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S",
                                      "T", "V", "W", "Y"])

    for sequence in sequences:
        for i, amino_acid in enumerate(sequence):
            frequencies.at[amino_acid, i] += 1

    return frequencies


def normalize_frequencies(frequencies):
    """
    Normalize the amino acid frequencies by dividing each count by the total number of sequences.
    """
    total_sequences = frequencies.sum(axis=0)
    normalized_frequencies = frequencies / total_sequences

    return normalized_frequencies


def visualize_frequency_matrix(frequencies):
    """
    Visualize the frequency matrix as a heatmap or color-coded table.
    """
    # Use a suitable visualization library (e.g., seaborn, matplotlib) to create a heatmap or color-coded table.
    # Here's an example using seaborn:
    import seaborn as sns
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    sns.heatmap(frequencies, cmap="coolwarm", linewidths=0.5)
    plt.xlabel("Position")
    plt.ylabel("Amino Acid")
    plt.title("Amino Acid Frequencies")
    plt.show()


pos_of_interest = [279, 280, 458, 459, 462, 369, 371]


def alnSiteCompositionDF(aln, characters="ACDEFGHIKLMNPQRSTVWY"):
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
    # Step 1: Parse sequence alignment
    alignment_file = "2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED.fasta"
    sequences = parse_sequence_alignment(alignment_file)
    # print(sequences)
    # Step 2: Calculate amino acid frequencies
    # frequencies = calculate_amino_acid_frequencies(sequences)
    cladeC_alignment = AlignIO.read(alignment_file, "fasta")
    frequencies = alnSiteCompositionDF(cladeC_alignment)
    print(frequencies)
    # for i in pos_of_interest:
    #     print(frequencies.iloc[i-1].to_frame().T)
    # subset_df = frequencies.iloc[pos_of_interest]  # Select the rows corresponding to the positions
    # print(subset_df)
    # print(df.iloc[278].to_frame().T)
    # Print the column names at the top
    # print("\t".join(frequencies.columns))
    # Print the selected rows

    # Step 3: Normalize frequencies
    # normalized_frequencies = normalize_frequencies(frequencies)
    # print(normalized_frequencies)

    # Step 4: Visualize frequency matrix
    # visualize_frequency_matrix(normalized_frequencies)

    # step 5: bar-chart
    # vis.generate_bar_chart(frequencies)


if __name__ == "__main__":
    main()
