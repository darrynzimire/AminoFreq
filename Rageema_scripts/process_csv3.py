
import argparse
import sys

import numpy as np
import pandas as pd
import warnings
from pprint import pprint
from dateutil.parser import parse

parser = argparse.ArgumentParser(description='csv_data_wrangling')
parser.add_argument('-f1', type=str, metavar='', required=True, help='Raw unprocessed HIV-1 data in csv format')
parser.add_argument('-f2', type=str, metavar='', required=True, help='SARS-COV2 data in csv format')


args = parser.parse_args()
f1 = args.f1
f2 = args.f2

date_formats = {
    "%Y-%m-%d": ['%Y-%m-%d', '%d-%m-%y', '%Y/%m/%d', '%d/%b/%Y']
}


# Function to convert the date value to the target format
def convert_date(date_str):
    for target_format, formats in date_formats.items():
        for date_format in formats:
            try:
                return pd.to_datetime(date_str, format=date_format).strftime(target_format)
            except ValueError:
                pass
    return date_str  # Return original value if no conversion is possible


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

    # df['date_column'] = df['date_column'].apply(convert_date)
    file1.iloc[:, 1] = file1.iloc[:, 1].apply(convert_date)
    # print('FILE1')
    # print(file1.iloc[:, 1].head(5))
    file2.iloc[:, 1] = file2.iloc[:, 1].apply(convert_date)
    # print('FILE 2')
    # print(file2.iloc[:, 1].head(5))
    file1.iloc[:, 1] = pd.to_datetime(file1.iloc[:, 1],  dayfirst=True).dt.date
    # print(file1.iloc[:, 1].head(5))
    file2.iloc[:, 1] = pd.to_datetime(file2.iloc[:, 1], dayfirst=True).dt.date
    # print(file2.iloc[:, 1].head(5))
    # print(type(file1.iloc[:, 1]))
    # Create an empty DataFrame to store the matched rows
    matched_rows = pd.DataFrame(columns=file1.columns)

    # Iterate through each sample_id1 in file 2
    for sample_id1, date2 in zip(file2['sample_id'], file2.iloc[:, 1]):
        # print(np.array(date2))
        # Find the matching rows in file 1 where sample_id1 matches sample_id
        matching_rows = file1[file1['sample_id'] == sample_id1]
        # If there are multiple matching rows, select the row with the closest date
        if len(matching_rows) > 1:
            closest_idx = (matching_rows.iloc[:, 1] - date2).abs().idxmin()
            # print(closest_idx)
            closest_row = matching_rows.loc[closest_idx]
            # print(closest_row)
            matched_rows = matched_rows._append(closest_row)
        elif len(matching_rows) == 1:
            matched_rows = matched_rows._append(matching_rows)

    merged_df = pd.merge(matched_rows, file2, on='sample_id', how='left')
    merged_df.iloc[:, [1, 3]] = merged_df.iloc[:, [1, 3]].apply(lambda column: pd.to_datetime(column).dt.date)
    # print(merged_df.columns)
    # print(merged_df.iloc[:, 3])
    # Calculate the difference in days and add it as a new column
    # Calculate the number of days between the two columns
    merged_df['days_difference'] = (merged_df.iloc[:, 1] - merged_df.iloc[:, 3]).dt.days
    print(type(merged_df.iloc[:, 3]))
    merged_df.iloc[:, 4] = merged_df.iloc[:, 4].abs()
    print(merged_df['days_difference'])

    # merged_df['no_of_days'] = merged_df.iloc[:, [1, 3]].apply(lambda row: (row.max() - row.min()).days, axis=1)
    # print(merged_df.head(10))

    return merged_df


def main():

    # Suppress FutureWarning
    # warnings.filterwarnings("ignore", category=FutureWarning)
    matched_df = find_matching_rows(f1, f2).reset_index(drop=True)
    matched_df.index +=1
    # print(matched_df)
    matched_df.to_csv('processed_viraldata.csv', index=False)


if __name__ == '__main__':
    main()