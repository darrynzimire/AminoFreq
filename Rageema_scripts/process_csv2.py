import argparse
import sys
import pandas as pd
import warnings

parser = argparse.ArgumentParser(description='csv_data_wrangling')
parser.add_argument('-f1', type=str, required=True)
parser.add_argument('-f2', type=str, required=True)

args = parser.parse_args()
f1 = args.f1
f2 = args.f2


def get_input(file):
    try:
        df = pd.read_csv(file, sep=',')
    except Exception:
        print(f'The input file {file} is not a valid csv')
        print('script terminated')
        sys.exit(1)

    return df


def find_matching_rows(f1, f2):
    file1 = get_input(f1)
    file2 = get_input(f2)
    file1['virus1_viral_load_test_date'] = pd.to_datetime(file1['virus1_viral_load_test_date'])
    file2['virus2_sample_collection_date'] = pd.to_datetime(file2['virus2_sample_collection_date'])

    # Create an empty DataFrame to store the matched rows
    matched_rows = pd.DataFrame(columns=file1.columns)

    # Iterate through each sample_id1 in file 2
    for sample_id1, date2 in zip(file2['sample_id'], file2['virus2_sample_collection_date']):
        # Find the matching rows in file 1 where sample_id1 matches sample_id
        matching_rows = file1[file1['sample_id'] == sample_id1]

        # If there are multiple matching rows, select the row with the closest date
        if len(matching_rows) > 1:
            closest_idx = (matching_rows['virus1_viral_load_test_date'] - date2).abs().idxmin()
            closest_row = matching_rows.loc[closest_idx]
            matched_rows = pd.concat([matched_rows, closest_row])
        elif len(matching_rows) == 1:
            matched_rows = pd.concat([matched_rows, matching_rows])

    return matched_rows


def main():
    # Suppress FutureWarning
    warnings.filterwarnings("ignore", category=FutureWarning)
    matched_df = find_matching_rows(f1, f2).reset_index(drop=True)
    matched_df.index += 1
    print(matched_df)
    matched_df.to_csv('processed_viraldata.csv', index=False)


if __name__ == '__main__':
    main()
