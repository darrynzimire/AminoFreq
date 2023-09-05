import csv
import argparse
import os


def remove_duplicates(input_file, output_file, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the input CSV file
    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        data = list(reader)

    # Create a dictionary to store unique entries and their corresponding dates
    unique_entries = {}

    # Process the data and remove duplicates
    for row in data:
        uct_id = row[0]
        test_date = row[1]

        # Check if the UCT ID already exists in the dictionary
        if uct_id in unique_entries:
            # Append the test date to the existing entry
            unique_entries[uct_id].append(test_date)
        else:
            # Add a new entry with the test date
            unique_entries[uct_id] = [test_date]

    # Write the unique entries with concatenated dates to the output CSV file
    with open(os.path.join(output_dir, output_file), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['uct_id', 'test_date'])  # Write the header row

        for uct_id, test_dates in unique_entries.items():
            # Concatenate the test dates with commas
            test_dates_str = ', '.join(test_dates)
            writer.writerow([uct_id, test_dates_str])


def main():
    remove_duplicates(args.input_file, args.output_file, args.output_dir)


if __name__ == '__main__':
    # Parse the command-line arguments
    print(os.getcwd())
    parser = argparse.ArgumentParser(description='Remove duplicates from a CSV file and concatenate dates.')
    parser.add_argument('--input_file',     type=str,   help='Input CSV file')
    parser.add_argument('--output_file',    type=str,   help='Output CSV file')
    parser.add_argument('--output_dir',     type=str,   help='Directory to write the output CSV file')
    args = parser.parse_args()
    i = args.input_file
    o = args.output_file
    od = args.output_dir

    # Call the function to remove duplicates and write the output
    # remove_duplicates(args.input_file, args.output_file, args.output_dir)
    remove_duplicates(i, o, od)


