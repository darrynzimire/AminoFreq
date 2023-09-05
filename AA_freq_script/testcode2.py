import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
from pprint import pprint
# Example amino acid frequencies for clade C virus sequences at different positions
amino_acid_freqs = {
    'Position1': {'Seq1': {'A': 10, 'C': 5, 'D': 7, 'E': 3}, 'Seq2': {'A': 8, 'C': 6, 'D': 5, 'E': 4}},
    'Position2': {'Seq1': {'A': 6, 'C': 5, 'D': 9, 'E': 2}, 'Seq2': {'A': 7, 'C': 5, 'D': 4, 'E': 6}}
}

# Store results for each position
results = []

# Iterate over each position
for position, freqs in amino_acid_freqs.items():
    # Create a contingency table for the position
    sequences = list(freqs.keys())
    amino_acids = sorted(list(set(freqs[sequences[0]].keys()) | set(freqs[sequences[1]].keys())))
    contingency_table = []
    for amino_acid in amino_acids:
        row = [freqs[sequence].get(amino_acid, 0) for sequence in sequences]
        contingency_table.append(row)
    pprint(contingency_table)
    # Perform chi-square test
    chi2, p_value, _, _ = chi2_contingency(contingency_table)

    # Store results
    results.append({'Position': position, 'Chi2': chi2, 'P-value': p_value})

# Convert results to a pandas DataFrame
results_df = pd.DataFrame(results)

# Set significance threshold
alpha = 0.05

# Identify positions with significant differences
significant_positions = results_df[results_df['P-value'] < alpha]['Position']

# Print significant positions
print("Significant positions:")
print(significant_positions)

# Plot p-values
plt.figure(figsize=(10, 6))
plt.plot(results_df['Position'], results_df['P-value'], 'bo-', label='P-value')
plt.axhline(y=alpha, color='r', linestyle='--', label='Significance Threshold')
plt.xlabel('Position')
plt.ylabel('P-value')
plt.title('P-values for Clade C Virus Sequences')
plt.legend()
plt.show()
