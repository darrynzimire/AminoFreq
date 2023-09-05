import argparse
import sys
import pandas as pd
import warnings
from dateutil.parser import parse

parser = argparse.ArgumentParser(description='csv_data_wrangling')
parser.add_argument('-f1', type=str, metavar='', required=True, help='Raw unprocessed HIV-1 data in csv format')
parser.add_argument('-f2', type=str, metavar='', required=True, help='SARS-COV2 data in csv format')


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


def sanity_check():

    # check that input files are valid csv files
    # check that both files contain the same number of rows
    # check that column names are correct
        # make script col name agnostic

    # check that data types in each column is correct/appropriate
    # add option to specify target directory
    # check if target directory have writing permissions
    # include a requirements.txt file to fend against dependencies errors
    # include date formatting options
    # option to name output file
    # use argv: no need for arguments
    pass

# def detect_date_format(date_column):
#     format = ["%Y-%m-%d", "%m-%d-%Y", '%d-%m-%Y']
#     for fmt in format:
#         try:
#             parse(date_column, fuzzy=False, dayfirst=True, yearfirst=False, default=None, ignoretz=True, format=[fmt])
#             return fmt
#         except ValueError:
#             pass
#         return None

# if file1_date_format is None:
#     print('Error: Unable to detect format in column')
#     return None


def detect_date_format(date_column):
    format = ["%Y-%m-%d", "%m-%d-%Y", '%d-%m-%Y', "%m/%d/%Y"]
    for fmt in format:
        print(fmt)
        try:
            parse(str(date_column), fuzzy=False, dayfirst=True, yearfirst=False, default=None, ignoretz=True)
            print(fmt)
            return fmt
        except ValueError:
            return None


def find_matching_rows(f1, f2):

    file1 = get_input(f1)
    file2 = get_input(f2)

    # file1_date_format = detect_date_format(file1.iloc[:, 1])
    # print(file1_date_format)
    # print(str(file1.iloc[:, 1]))
    # print(file1_date_format)
    # if file1_date_format is None:
    #     print(f"Error: Unable to detect the date format in {f1} column of file")
    file1['virus1_viral_load_test_date'] = pd.to_datetime(file1['virus1_viral_load_test_date'], infer_datetime_format=True, dayfirst=True)
    # file2_date_format = detect_date_format(file2.iloc[:, 1])
    # if file2_date_format is None:
    #     print(f"Error: Unable to detect the date format in {f2} column of file")
    file2['virus2_sample_collection_date'] = pd.to_datetime(file2['virus2_sample_collection_date'], infer_datetime_format=True, dayfirst=True)
    # print(file2['virus2_sample_collection_date'])
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
            matched_rows = matched_rows._append(closest_row)
            # matched_rows = pd.concat(closest_row)
        elif len(matching_rows) == 1:
            matched_rows = matched_rows._append(matching_rows)
            # matched_rows = pd.concat(matching_rows)
    # matched_rows['virus2_sample_collection_date'] = pd.to_datetime(file2['virus2_sample_collection_date'])
    merged_df = pd.merge(matched_rows, file2, on='sample_id', how='left')
    merged_df.iloc[:, [1, 3]] = merged_df.iloc[:, [1, 3]].apply(lambda column: pd.to_datetime(column))

    # merged_df['no_of_days'] = (merged_df['virus2_sample_collection_date'] - merged_df['virus2_sample_collection_date']).dt.days
    # Calculate the difference in days and add it as a new column
    merged_df['no_of_days'] = merged_df.iloc[:, [1,3]].apply(lambda row: (row.max() - row.min()).days, axis=1)

    # print(merged_df.info())
    return merged_df
    # print(matching_rows)
    # return matched_rows


def main():

    # Suppress FutureWarning
    warnings.filterwarnings("ignore", category=FutureWarning)
    matched_df = find_matching_rows(f1, f2).reset_index(drop=True)
    matched_df.index +=1
    # print(matched_df)
    matched_df.to_csv('processed_viraldata.csv', index=False)


if __name__ == '__main__':
    main()