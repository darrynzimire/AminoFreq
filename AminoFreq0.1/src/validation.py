from Bio import AlignIO
import os


def check_fasta_alignment(file_path):
    try:
        alignment = AlignIO.read(file_path, 'fasta')

        # Check sequence lengths
        seq_lengths = [len(seq.seq) for seq in alignment]
        if len(set(seq_lengths)) != 1:
            return False, "Sequences have varying lengths"

        # Check if sequences are aligned
        if any("-" in seq.seq for seq in alignment):
            return False, "Sequences are not aligned"

        return True, "Alignment is valid"

    except Exception as e:
        return False, f"Error: {str(e)}"


def is_csv_file(file_path):
    _, file_extension = os.path.splitext(file_path)
    return file_extension.lower() == '.csv'

