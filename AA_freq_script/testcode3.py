import logomaker
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import seqlogo
from Bio import AlignIO
from scipy import stats
import logomaker as lm

f = 'fasta_testoutfile.fasta'
f2 = '2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED.fasta'
f3 = '2014_cladeC_HXBR_AY2_TRANSLATED.fasta'
f4 = '2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED_for_logo.fasta'


def extract_sequences_from_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        current_sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(str(current_sequence))
                    current_sequence = ''
            else:
                current_sequence += line
        if current_sequence:
            sequences.append(str(current_sequence))
    return sequences


# print(extract_sequences_from_fasta(f4))


def generate_sequence_logo(aligned_sequences):
    # Calculate the frequency matrix from the aligned sequences
    frequency_matrix = lm.alignment_to_matrix(aligned_sequences)
    # print(frequency_matrix)
    # counts_mat = lm.alignment_to_matrix(seqs)
    print(frequency_matrix.head())
    # Create a sequence logo using the frequency matrix
    # logo = logomaker(frequency_matrix)
    lm.Logo(frequency_matrix)

    # Customize the appearance of the sequence logo
    # lm.style_title('Sequence Logo')
    # lm.style_color('classic')
    # lm.style_show_consensus(True)
    # lm.format_xticks(format='digit')

    # Display the sequence logo
    p = lm.sequence_to_matrix(aligned_sequences)
    return p


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


seq = ''.join(extract_sequences_from_fasta(f4))
print(seq)
# visualize_frequency_matrix(generate_sequence_logo(extract_sequences_from_fasta(seq)))

# generate_sequence_logo(extract_sequences_from_fasta(f4))
# t7_alignment = AlignIO.read(f4, "fasta")
# print(t7_alignment)


# def alnSiteCompositionDF(aln, characters="ACDEFGHIKLMNPQRSTVWY"):
#   alnRows = aln.get_alignment_length()
#   compDict = {char:[0]*alnRows for char in characters}
#   for record in aln:
#     header = record.id
#     seq = record.seq
#     for aaPos in range(len(seq)):
#       aa = seq[aaPos]
#       if aa in characters:
#         compDict[aa][aaPos] += 1
#   return pd.DataFrame.from_dict(compDict)


# t7_alignmentSiteCompDF = alnSiteCompositionDF(t7_alignment)
# print(t7_alignmentSiteCompDF)
# generate_sequence_logo(f)

# t7_alignmentSiteFreqDF = t7_alignmentSiteCompDF.div(t7_alignmentSiteCompDF.sum(axis=1), axis=0)
# print(t7_alignmentSiteFreqDF)

# t7_alignmentSiteFreqSeqLogo= seqlogo.Ppm(t7_alignmentSiteFreqDF, alphabet_type="AA")
# print(t7_alignmentSiteFreqSeqLogo)
# seqlogo.seqlogo(t7_alignmentSiteFreqSeqLogo, ic_scale=False, format='svg', size='large')


# Chi-square analysis of MSA
# alignment_file = "path/to/alignment.fasta"
# alignment = AlignIO.read(open(f3), "fasta")


# def compositionMatrix(aln):
# 	compDict = {}
# 	fixedCharacters = ["-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
# 	for record in aln:
# 		header = record.id
# 		seq = record.seq
# 		currentSeqMat = [0]*21
# 		for i in range(len(seq)):
# 			aa = seq[i]
# 			try:
# 				characterPos = fixedCharacters.index(aa)
# 				currentSeqMat[characterPos]+= 1
# 			except:
# 				print("ERROR:", header, "contains character ("+aa+") not in the list:",fixedCharacters)
# 		compDict[header] = currentSeqMat
# 	compDF = pd.DataFrame.from_dict(compDict, orient='index',
# 					   columns=fixedCharacters)
# 	return compDF
#
#
# def chi2test(compDF):
#     seqTotals = compDF.sum(axis=1)
#     gaps = compDF["-"]
#     gapsPerSeq = gaps / seqTotals
#
#     nonGap = compDF.loc[:, 'A':'Y']
#     nonGapTotals = nonGap.sum().to_frame()
#     nonGapSeqTotals = nonGap.sum(axis=1).to_frame()
#     numCharacters = nonGapTotals.sum()
#     expectedFreq = nonGapTotals / numCharacters
#
#     expectedCountArray = np.dot(nonGapSeqTotals, expectedFreq.transpose())
#     expectedCountDF = pd.DataFrame(expectedCountArray, columns=nonGap.columns, index=nonGap.index.values)
#     chi2DF = ((expectedCountDF - nonGap) ** 2) / expectedCountDF
#     chi2Sum = chi2DF.sum(axis=1)
#
#     pValueDf = 1 - stats.chi2.cdf(chi2Sum, 19)
#     outDF = pd.DataFrame({"Gap/Ambiguity": gapsPerSeq, "p-value": pValueDf})
#     outDF.index.name = 'header'
#
#     return outDF
#
#
# matrix = compositionMatrix(alignment)
# result = chi2test(matrix)
# print(result)
#
# inf = '2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED_for_logo.data.txt'
#
# data = []
# with open(inf) as file:
#     for line in file:
#         line = line.strip()
#         if line.startswith('#') or not line:
#             continue
#         values = line.split('\t')
#         data.append([int(values[0])] + [int(count) for count in values[1:-4]])
#
# df = pd.DataFrame(data, columns=['Position'] + list('ACDEFGHIKLMNPQRSTVWY'))
#
# print(df)
# # Calculate the sum of frequencies for each amino acid
# aa_frequencies = df.iloc[:, 1:].sum()
#
# # Plot the amino acid frequencies
# plt.figure(figsize=(10, 6))
# plt.bar(aa_frequencies.index, aa_frequencies.values)
# plt.xlabel('Amino Acid')
# plt.ylabel('Frequency')
# plt.title('Amino Acid Frequencies')
# plt.xticks(rotation=45)
# plt.show()
#
# # Plotting heatmap
#
# #Transpose the DataFrame to have positions as rows and amino acids as columns
# df_transposed = df.iloc[:, 1:].transpose()
#
# # Plot the amino acid frequencies at each position
# plt.figure(figsize=(12, 8))
# sns.heatmap(df_transposed, cmap='YlGnBu', linewidths=0.5)
# plt.xlabel('Position')
# plt.ylabel('Amino Acid')
# plt.title('Amino Acid Frequencies at Each Position')
# plt.xticks(rotation=45)
# plt.show()

# plotting a bar plot
#
# import matplotlib.pyplot as plt
# import numpy as np
#
# # Get the positions from the 'Position' column
# positions = df['Position'].values.astype(int)
#
# # Extract the amino acid columns from the DataFrame
# amino_acids = df.columns[1:]
#
# # Get the frequencies for each amino acid as a 2D numpy array
# frequencies = df[amino_acids].values.astype(float)
#
# # Plot the amino acid frequencies at each position
# plt.figure(figsize=(10, 6))
#
# for i in range(len(amino_acids)):
#     plt.bar(positions, frequencies[:, i], label=amino_acids[i])
#
# plt.xlabel('Position')
# plt.ylabel('Frequency')
# plt.title('Amino Acid Frequencies at Each Position')
# plt.legend()
# plt.xticks(rotation=45)
# plt.show()

# import matplotlib.pyplot as plt
# import numpy as np
#
# # Get the positions from the 'Position' column
# positions = df['Position'].values.astype(int)
#
# # Extract the amino acid columns from the DataFrame
# amino_acids = df.columns[1:]
#
# # Get the frequencies for each amino acid as a 2D numpy array
# frequencies = df[amino_acids].values.astype(float)
#
# # Plot the amino acid frequencies at each position
# plt.figure(figsize=(10, 6))
#
# for i in range(len(amino_acids)):
#     plt.bar(positions, frequencies[:, i], label=amino_acids[i])
#
# plt.xlabel('Position')
# plt.ylabel('Frequency')
# plt.title('Amino Acid Frequencies at Each Position')
#
# # Position the legend next to the graph
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.xticks(rotation=45)
# plt.tight_layout()  # Adjust the layout to prevent overlapping
#
# plt.show()
#
# # write a table
#
# # Create a new DataFrame for the table
# table_data = {'Position': positions}
# for i, aa in enumerate(amino_acids):
#     table_data[aa] = frequencies[:, i]
#
# table_df = pd.DataFrame(table_data)
#
# # Display the table
# print(table_df)
#
#
#
