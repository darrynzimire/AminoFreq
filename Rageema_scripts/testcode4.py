import pandas as pd

f = 'virus_1_input_data.csv'
f2 = 'virus_2_input_data.csv'

df = pd.read_csv(f, parse_dates=[1],
                 date_parser=lambda x: pd.to_datetime(x, format='%Y/%m/%d'))

print(df)
# df = pd.read_csv(fn, parse_dates=[0],
#                           date_parser=lambda x: pd.to_datetime(x, format='%m/%d/%Y %I:%M:%S %p'))
