from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import numpy as np
import math

def summary_stats():
    pass


def compare_groups(control_freq, vrc01_freq, specific_sites, alpha=0.05):
    p_values = []
    for site in specific_sites:
        control_mutations = list(control_freq[site].keys())
        vrc01_mutations = list(vrc01_freq[site].keys())
        all_mutations = sorted(list(set(control_mutations + vrc01_mutations)))

        # Create contingency table
        contingency_table = []
        for mutation in all_mutations:
            control_count = control_freq[site].get(mutation, 0)
            vrc01_count = vrc01_freq[site].get(mutation, 0)
            contingency_table.append([control_count, vrc01_count])

        # Perform chi-square test
        _, p_value, _, _ = chi2_contingency(contingency_table)

        p_values.append(p_value)

    # Adjust p-values for multiple testing using the Benjamini-Hochberg procedure
    adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]

    significant_sites = [site for site, p_value in zip(specific_sites, adjusted_p_values) if p_value < alpha]

    return adjusted_p_values, significant_sites


def calculate_amino_acid_entropy(frequencies):
    """
    Calculate amino acid entropy based on frequencies.

    Args:
        frequencies (dict): A dictionary of amino acid frequencies.

    Returns:
        float: Amino acid entropy value.
    """
    total_frequency = sum(frequencies.values())
    entropy = 0.0
    for frequency in frequencies.values():
        if frequency > 0:
            p = frequency / total_frequency
            entropy -= p * np.log2(p)
    return entropy


def amino_acid_entropy(frequencies):
    """
    Calculate the amino acid entropy based on frequencies.

    Parameters:
        frequencies (dict): A dictionary containing amino acid frequencies.
                            The keys are amino acid symbols (single-letter codes)
                            and the values are the corresponding frequencies.

    Returns:
        float: The calculated amino acid entropy.

        # Example usage:
        # amino_acid_frequencies = {'A': 20, 'C': 15, 'G': 12, 'T': 18, 'N': 5}
    """
    # Calculate the total number of amino acids
    total_amino_acids = sum(frequencies.values())

    # Calculate the entropy
    entropy = 0
    for frequency in frequencies.values():
        if frequency > 0:
            probability = frequency / total_amino_acids
            entropy -= probability * math.log2(probability)

    return entropy
