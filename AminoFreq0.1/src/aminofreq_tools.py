import pandas as pd


def identify_conserved_regions(data_df, threshold=0.9):
    """
    Identify conserved regions in the HIV-1 envelope protein.

    Args:
        data_df (pd.DataFrame): A DataFrame containing amino acid frequencies at different positions.
        threshold (float, optional): The threshold for amino acid conservation. Defaults to 0.9.

    Returns:
        list: A list of tuples representing the conserved regions (start_position, end_position).
    """
    conserved_regions = []
    positions = data_df.index
    num_sequences = len(data_df)

    for pos in positions:
        frequencies = data_df.loc[pos]
        max_frequency = frequencies.max()
        most_common_aa = frequencies.idxmax()
        conservation_score = max_frequency / num_sequences

        if conservation_score >= threshold:
            if not conserved_regions or pos != conserved_regions[-1][1] + 1:
                conserved_regions.append((pos, pos))
            else:
                conserved_regions[-1] = (conserved_regions[-1][0], pos)

    return conserved_regions


def identify_highly_variable_positions(data_df, threshold=0.1):
    """
    Identify positions in the HIV-1 envelope protein with high variability.

    Args:
        data_df (pd.DataFrame): A DataFrame containing amino acid frequencies at different positions.
        threshold (float, optional): The threshold for identifying highly variable positions. Defaults to 0.1.

    Returns:
        list: A list of positions that show a high frequency of mutations or variability.
    """
    highly_variable_positions = []
    positions = data_df.index

    for pos in positions:
        frequencies = data_df.loc[pos]
        variability = frequencies.std()

        if variability >= threshold:
            highly_variable_positions.append(pos)

    return highly_variable_positions


def compare_to_consensus(frequencies, consensus_sequence):
    """
    Compare the amino acid frequencies at specific positions to the consensus sequence.

    Args:
        frequencies (dict): A dictionary containing the amino acid frequencies at specific positions.
        consensus_sequence (str): The consensus sequence (HXB2) for comparison.

    Returns:
        dict: A dictionary containing the comparison results for each position.
            - 'match': Number of amino acids that match the consensus sequence.
            - 'mismatch': Number of amino acids that mismatch the consensus sequence.
            - 'frequency_match': Frequency of amino acids that match the consensus sequence.
            - 'frequency_mismatch': Frequency of amino acids that mismatch the consensus sequence.
    """
    result = {
        'match': 0,
        'mismatch': 0,
        'frequency_match': 0,
        'frequency_mismatch': 0
    }

    for position, frequency in frequencies.items():
        consensus_aa = consensus_sequence[position - 1]
        if frequency > 0:
            if consensus_aa == '-':
                result['mismatch'] += 1
                result['frequency_mismatch'] += frequency
            elif consensus_aa in frequency:
                result['match'] += 1
                result['frequency_match'] += frequency
            else:
                result['mismatch'] += 1
                result['frequency_mismatch'] += frequency

    return result


def find_unique_amino_acids(frequencies, consensus_sequence):
    """
    Find unique or signature amino acid residues present in the Clade C panel but not in the reference sequence (HXB2).

    Args:
        frequencies (dict): A dictionary containing the amino acid frequencies at specific positions in the Clade C panel.
        consensus_sequence (str): The consensus sequence (HXB2) for comparison.

    Returns:
        dict: A dictionary containing the unique or signature amino acid residues for each position.
            - 'unique_residues': List of amino acid residues unique to the Clade C panel at each position.
            - 'unique_frequencies': List of corresponding frequencies of unique amino acids.
    """
    unique_residues = []
    unique_frequencies = []

    for position, frequency in frequencies.items():
        consensus_aa = consensus_sequence[position - 1]
        if frequency > 0 and consensus_aa != '-' and consensus_aa not in frequency:
            unique_residues.append(consensus_aa)
            unique_frequencies.append(frequency)

    return {
        'unique_residues': unique_residues,
        'unique_frequencies': unique_frequencies
    }


def find_high_gap_frequency_positions(frequencies):
    """
    Find positions in the envelope protein with a high frequency of insertions or deletions (gaps) in the Clade C panel.

    Args:
        frequencies (dict): A dictionary containing the amino acid frequencies at specific positions in the Clade C panel.

    Returns:
        list: A list of positions with high gap frequency.
    """
    high_gap_positions = []

    for position, frequency in frequencies.items():
        if '-' in frequency:
            gap_frequency = frequency['-']
            if gap_frequency > 0.5:  # You can adjust the threshold for considering a position to have a high gap frequency
                high_gap_positions.append(position)

    return high_gap_positions




