from collections import Counter


def calculate_amino_acid_frequencies(fasta_file):
    """
    Calculate the amino acid frequencies at all positions in an aligned FASTA file.
    :param fasta_file: Path to the aligned FASTA file.
    :return: A dictionary containing the amino acid frequencies at each position.
    """
    frequencies = {}

    with open(fasta_file, 'r') as file:
        lines = file.readlines()

        # Initialize the frequencies dictionary
        for line in lines:
            if line.startswith('>'):
                seq_len = len(lines[lines.index(line) + 1].strip())
                # print(seq_len)
                frequencies.update({i: Counter() for i in range(seq_len)})

        # Calculate the frequencies
        for i in range(seq_len):
            for line in lines:
                # print(line)
                if not line.startswith('>'):
                    frequencies[i].update(line[i])

    return frequencies


# Usage example
fasta_file = '2014_cladeCpanel_V703_synthseq_HXBR_AY_TRANSLATED.fasta'
frequencies = calculate_amino_acid_frequencies(fasta_file)
print(frequencies)

# Accessing the frequencies
for position, counts in frequencies.items():
    print(f"Position {position + 1}: {counts}")
