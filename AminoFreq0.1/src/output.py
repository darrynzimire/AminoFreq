import datetime
import os


def generate_outputdir(name, zipped):

    directory = '../amino_frequency_output'
    if not os.path.isdir(directory):
        os.mkdir(directory)

    file_path = os.path.join(directory, name)
    return file_path


def generate_report(output1_data, output2_file):
    _freq = generate_outputdir(output2_file, zipped=False)
    #
    with open(_freq, "w+") as f:

    #     # Add the date (machine date) when the script was run
        date_today = datetime.date.today()
        f.write(f"Date: {date_today}\n")
        f.write("query sites:\n")
        f.write("Reference: HXB2\n")
    #
        for data in output1_data:
            query_position = data[0]
            # print(query_position)
            aligned_position = data[1]
            # print(aligned_position)
    #         reference_residue = data[2]
    #         amino_acid_frequency = data[3]
    #
    #         f.write(f"\nQuery Position: {query_position}\n")
    #         f.write(f"Aligned Position: {aligned_position}\n")
    #         f.write(f"Reference Residue at position {query_position}: {reference_residue} (HXB2 numbering)\n")
    #         f.write(f"Amino acid frequency at position {query_position}: {amino_acid_frequency:.2f}%\n\n")
    #
    #         # Assuming the data is in the format (aligned_position, amino_acid_frequencies)
    #         amino_acid_changes = dict(data[4])  # Convert set to dictionary
    #         total_sequences = data[5]
    #
    #         f.write("Amino acid Change:\tFrequency:\t\tno. of sequences\n")
    #         for change, frequency in amino_acid_changes.items():
    #             percentage = (frequency / total_sequences) * 100
    #             num_sequences_with_change = frequency
    #             f.write(f"{reference_residue}{query_position}{change:<15}{percentage:.1f}%\t\t\t({num_sequences_with_change}/{total_sequences})\n")
    #

def generate_logging():
    pass

