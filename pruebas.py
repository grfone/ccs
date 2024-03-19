import pandas
import bz2
import pickle
import pandas as pd


# Read the CSV files
fingerprints = pandas.read_csv('./resources/METLIN_CCS_vectorfingerprintsVectorized.csv')
fingerprints.dropna(subset=['Molecule Name'], inplace=True)
fingerprints.drop(columns=['dimer line', 'CCS', 'm/z.2', 'pubChem', 'METLIN ID',], inplace=True)
fingerprints.drop(fingerprints[fingerprints['Molecule Name'].str.contains("Tm_")].index, inplace=True)

# descriptors = pandas.read_csv('./resources/METLIN_CCS_descriptors.csv')
# descriptors.dropna(subset=['Molecule Name'], inplace=True)
# descriptors.drop(columns=['dimer line', 'CCS', 'm/z.2', 'pubChem', 'METLIN ID',], inplace=True)
# descriptors.drop(descriptors[descriptors['Molecule Name'].str.contains("Tm_")].index, inplace=True)


def calculate_average_ccs(row):
    # Extract values from columns CCS1, CCS2, CCS3 for the current row
    ccs1 = row['CCS1']
    ccs2 = row['CCS2']
    ccs3 = row['CCS3']

    # Calculate the average and round to two decimals
    average_ccs = round((ccs1 + ccs2 + ccs3) / 3, 2)

    return average_ccs


fingerprints['rt'] = fingerprints.apply(calculate_average_ccs, axis=1)
# descriptors['rt'] = descriptors.apply(calculate_average_ccs, axis=1)

fingerprints['pid'] = range(len(fingerprints['rt']))
# descriptors['pid'] = range(len(descriptors['rt']))

columns_to_drop = ['Molecule Name', 'Molecular Formula', 'Precursor Adduct', 'CCS1', 'CCS2', 'CCS3', 'CCS_AVG',
                   '% CV', 'm/z', 'Adduct', 'm/z.1', 'Dimer', 'Dimer.1', 'inchi', 'smiles', 'InChIKEY']

fingerprints.drop(columns=columns_to_drop, inplace=True)
# descriptors.drop(columns=columns_to_drop, inplace=True)

# fingerprints shape (61863, 2216)
# descriptors shape (61863, 5668)



# Make Adduct binary by creating three columns and removing the original
#-----------------------------------------------------------------------
# Read metlin database from resources
import pandas as pd
metlin_df = pd.read_csv('resources/Metlin.csv')
# Remove excess columns and rows from the Metlin DataFrame.
metlin_df.drop(columns=['dimer line', 'CCS', 'm/z.2'], inplace=True)
metlin_df.drop(metlin_df.index[-5:], inplace=True)
metlin_df.drop(metlin_df[metlin_df['Molecule Name'].str.contains("Tm_")].index, inplace=True)

# Use pandas.get_dummies() to create binary columns
adduct_dummies = pandas.get_dummies(metlin_df['Adduct'])

# Rename the columns to match the specified conditions
adduct_dummies.columns = ['[M+H]', '[M-H]', '[M+Na]']

# Concatenate the original DataFrame with the new binary columns
metlin_df = pandas.concat([metlin_df, adduct_dummies], axis=1)

# Fill NaN values with 0 (for cases where the original 'Adduct' didn't match any condition)
metlin_df = metlin_df.fillna(0)

# Add columns to fingerprints
print(fingerprints.shape)
fingerprints[['[M+H]', '[M-H]', '[M+Na]']] = metlin_df[['[M+H]', '[M-H]', '[M+Na]']]
print(fingerprints.shape)


with bz2.BZ2File('./results/fingerprints.pklz', 'wb') as f:
    pickle.dump(fingerprints, f)


# with bz2.BZ2File('./results/descriptors.pklz', 'wb') as f:
#     pickle.dump(descriptors, f)


# Assuming 'descriptors' is your DataFrame
# Save DataFrame to CSV
# descriptors.to_csv('./results/descriptors.csv', index=False)


"""
smrt = pandas.read_csv('./resources/descriptors_smrt.csv')
print(smrt.head())
for pid in [25345055, 17541371, 16295966]:
    if (smrt['pid'] == pid).any():
        print(f"There is at least one value in the 'pid' column equal to {pid}.")
    else:
        print(f"There is no value in the 'pid' column equal to {pid}.")





import sys
sys.exit()


compuesto86 = pandas.read_csv('./resources/new_file_sent_by_alberto.csv', sep='\t')
compuesto86.drop(compuesto86.columns[0], axis=1, inplace=True)
num_columns = compuesto86.shape[1]
max_columns_per_chunk = 1024
num_chunks = (num_columns + max_columns_per_chunk - 1) // max_columns_per_chunk
for chunk_number in range(num_chunks):
    start_column = chunk_number * max_columns_per_chunk
    end_column = (chunk_number + 1) * max_columns_per_chunk
    chunk_df = compuesto86.iloc[:, start_column:end_column]
    chunk_df.to_csv(f'./results/comparison_chunk_{chunk_number + 1}.csv', index=False)


"""

