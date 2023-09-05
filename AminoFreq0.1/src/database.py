import sqlite3
import pandas as pd


def create_table_from_dataframe(df, table_name, db_file='CladeCref.db'):
    """
    Create a table in the SQLite database from a DataFrame.

    :param df: The DataFrame to be stored in the database.
    :param table_name: Name of the table to be created.
    :param db_file: Path to the SQLite database file. Default is 'database.db'.
    """
    conn = sqlite3.connect(db_file)
    df.to_sql(table_name, conn, if_exists='replace', index=False)
    conn.close()


def read_dataframe_from_table(table_name, db_file='CladeCref.db'):
    """
    Read a DataFrame from the SQLite database.

    :param table_name: Name of the table to read from.
    :param db_file: Path to the SQLite database file. Default is 'database.db'.
    :return: The DataFrame read from the table.
    """
    conn = sqlite3.connect(db_file)
    df = pd.read_sql_query(f'SELECT * FROM {table_name}', conn)
    conn.close()
    return df






# Connect to the database (if it doesn't exist, it will be created)
# conn = sqlite3.connect('CladeC_frequencies.db')
# cursor = conn.cursor()
#
# # Create a table to store the position-specific amino acid frequencies
# cursor.execute('''CREATE TABLE IF NOT EXISTS amino_acids (
#                     position INTEGER,
#                     amino_acid TEXT,
#                     frequency REAL,
#                     PRIMARY KEY (position, amino_acid)
#                 )''')
#
# conn.commit()
#
# # Assuming you have a list of tuples containing the data like this:
# data = [
#     (1, 'A', 0.05),
#     (1, 'C', 0.02),
#     (1, 'G', 0.08),
#     (1, 'T', 0.85),
#     # More data for other positions
# ]
#
# # Insert the data into the database
# cursor.executemany('INSERT INTO amino_acids (position, amino_acid, frequency) VALUES (?, ?, ?)', data)
#
# conn.commit()
#
#
# def get_amino_acid_frequency_for_position(position, amino_acid):
#     cursor.execute('SELECT frequency FROM amino_acids WHERE position = ? AND amino_acid = ?', (position, amino_acid))
#     result = cursor.fetchone()
#     return result[0] if result else None
#
# # Example usage
# position = 1
# amino_acid = 'A'
# frequency = get_amino_acid_frequency_for_position(position, amino_acid)
# print(f"The frequency of amino acid '{amino_acid}' at position {position} is: {frequency}")
#
# conn.close()
