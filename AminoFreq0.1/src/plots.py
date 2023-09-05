import matplotlib.pyplot as plt
import pandas as pd
from io import BytesIO
import matplotlib.cm as cm
import output
from weblogo import *
import json
import sys


def make_logogram(counts_mat, alphabet_map, sites, config_path):

    with open(config_path, 'r') as config_file:
        config = json.load(config_file)
    logodata = LogoData.from_counts(alphabet_map, counts_mat)
    logooptions = LogoOptions()

    for key, value in config.items():
        setattr(logooptions, key, value)
    logooptions.annotate = sites
    # logooptions.logo_end = counts_mat.shape[0]
    color_scheme  = config.get("color_scheme")
    if color_scheme == "hydrophobicity":
        logooptions.color_scheme = hydrophobicity
    if color_scheme == "chemistry":
        logooptions.color_scheme = chemistry
    if color_scheme == "charge":
        logooptions.color_scheme = charge

    format = LogoFormat(logodata, logooptions)
    print(csv(logodata))
    # Create PNG logo
    png = png_formatter(logodata, format)
    # png = weblogo.logo_formatter.png_formatter(logodata, format)

    # Save PNG logo to a file
    _f = output.generate_outputdir('logogram.png', zipped=False)
    with open(_f, 'wb') as png_file:
        png_file.write(png)
# #
# #     # Create PDF logo
#     pdf = pdf_formatter(logodata, format)
# #
# #     # Save PDF logo to a file
# #     _f1 = output.generate_outputdir('position_logo.pdf', zipped=False)
#     with open(_f, 'wb') as pdf_file:
#         pdf_file.write(pdf)


def visualize_stacked_bar_charts(positions, amino_acid_freqs):
    """
    Visualize amino acid frequencies at specific site positions using stacked bar charts.

    Parameters:
        positions (list): A list of site positions.
        amino_acid_freqs (list of dict): A list of dictionaries containing amino acid frequencies.
                                         Each dictionary has keys as amino acid residues and values as frequencies.

    Returns:
        None (displays the plots)
    """
    num_positions = len(positions)

    # Create a colormap with a different color for each residue
    colors = cm.tab20b(range(20))

    # Create a subplot for each position
    fig, axs = plt.subplots(nrows=num_positions, ncols=1, figsize=(10, 6*num_positions))

    # Generate stacked bar charts for each position
    for i, position in enumerate(positions):
        amino_acid_freq = amino_acid_freqs[i]

        # Extract amino acids and their frequencies from the dictionary
        amino_acids = list(amino_acid_freq.keys())
        frequencies = list(amino_acid_freq.values())

        # Create a stacked bar chart for the current position
        axs[i].bar(amino_acids, frequencies, color=colors)

        # Set plot title and labels
        axs[i].set_title(f"Amino Acid Frequencies at Position {position}")
        axs[i].set_xlabel("Amino Acid Residue")
        axs[i].set_ylabel("Frequency (%)")

    # Adjust layout and display the plots
    plt.tight_layout()
    plt.show()


# volcano plot
def generate_volcano_plot(reference_freq, treatment_freq, significance_threshold=0.05):
    """
    Generate a Volcano plot to visualize statistically significant differences in amino acid frequencies.

    :param reference_freq: A dictionary containing amino acid frequencies in the reference arm.
                           Format: {amino_acid: frequency}
    :param treatment_freq: A dictionary containing amino acid frequencies in the treatment arm.
                           Format: {amino_acid: frequency}
    :param significance_threshold: The significance threshold for the statistical test (default: 0.05).
    :return: None (the plot will be displayed or saved based on your configuration)
    """

    # Calculate the log-fold change in frequencies
    log_fold_change = {}
    for aa in reference_freq:
        if aa in treatment_freq and reference_freq[aa] > 0 and treatment_freq[aa] > 0:
            log_fold_change[aa] = abs(treatment_freq[aa] - reference_freq[aa])

    # Separate significant and non-significant points
    significant_points = {aa: log_fold_change[aa] for aa in log_fold_change if log_fold_change[aa] > significance_threshold}
    non_significant_points = {aa: log_fold_change[aa] for aa in log_fold_change if log_fold_change[aa] <= significance_threshold}

    # Create the Volcano plot
    plt.figure(figsize=(10, 6))
    plt.scatter(list(non_significant_points.values()), list(non_significant_points.keys()), color='grey', alpha=0.5)
    plt.scatter(list(significant_points.values()), list(significant_points.keys()), color='red', alpha=0.8)
    plt.axvline(significance_threshold, color='black', linestyle='dashed')
    plt.xlabel('Log Fold Change')
    plt.ylabel('Amino Acid')
    plt.title('Volcano Plot: Significant Differences in Amino Acid Frequencies')
    plt.legend(['Threshold', 'Non-Significant', 'Significant'])
    plt.tight_layout()
    # Display or save the plot (modify this based on your preference)
    plt.show()


# Generate a bar-plot
def bar_plot(position, frequencies):
    residues = list(frequencies.keys())
    frequencies = list(frequencies.values())

    plt.bar(residues, frequencies)
    plt.xlabel('Amino Acid Residue')
    plt.ylabel('Frequency')
    plt.title(f'Amino Acid Frequencies at Position {position}')
    plt.xticks(rotation=90)
    plt.show()
