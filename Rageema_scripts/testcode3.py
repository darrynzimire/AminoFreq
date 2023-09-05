import csv
import argparse
import os


def remove_duplicates(input_file):

    os.makedirs(os.getcwd(), exist_ok=True)

    # Read the input CSV file
    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        data = list(reader)

    unique_entries = {}

    for row in data:
        uct_id = row[0]
        test_date = row[1]
        test_value = row[2]
        if uct_id in unique_entries:

            unique_entries[uct_id][0].append(test_date)
            unique_entries[uct_id][1].append(test_value)
        else:

            unique_entries[uct_id] = ([test_date], [test_value])

    with open(os.path.join(os.getcwd(), '{}{}'.format(input_file, '_proc.csv')), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['uct_id', 'test_date', 'test_value'])  # Write the header row

        for uct_id, (test_dates, test_values) in unique_entries.items():
            test_dates_str = ', '.join(test_dates)
            test_values_str = ', '.join(test_values)
            writer.writerow([uct_id, test_dates_str, test_values_str])


if __name__ == '__main__':
    # Parse the command-line arguments
    parser = argparse.ArgumentParser(description='Remove duplicates from a CSV file and concatenate dates and test values.')
    parser.add_argument('--input_file', help='Input CSV file')
    # parser.add_argument('--output_file', help='Output CSV file')
    parser.add_argument('--output_dir', help='Directory to write the output CSV file')
    args = parser.parse_args()
    i = args.input_file
    remove_duplicates(i)
