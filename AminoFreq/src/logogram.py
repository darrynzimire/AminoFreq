import Bio
import sys
import argparse
import numpy as np
from weblogo import LogoData, LogoOptions, LogoFormat, pdf_formatter, png_formatter


def make_logogram(df, positions):
    selected_rows = df.loc[positions]
    counts_matrix = selected_rows.values
    transposed_counts_matrix = np.transpose(counts_matrix)
    residues = ''.join(selected_rows.columns.tolist())
    logodata = LogoData.from_counts(residues, transposed_counts_matrix)
    logooptions = LogoOptions()
    logooptions.title = "Position Logo"
    #Create format
    logoformat = LogoFormat(logodata, logooptions)
    # Create PNG logo
    png = png_formatter(logodata, logoformat)

    # Save PNG logo to a file
    with open('positions_logo.png', 'wb') as png_file:
        png_file.write(png)

    # Create PDF logo
    pdf = pdf_formatter(logodata, logoformat)

    # Save PDF logo to a file
    with open('positions_logo.pdf', 'wb') as pdf_file:
        pdf_file.write(pdf)


