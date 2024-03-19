import pandas
import pandas as pd


# Load datasets
descriptors_Metlin = pandas.read_csv('./resources/METLIN_CCS_descriptors.csv')
descriptors_SMRT = pandas.read_csv('./resources/descriptors_smrt.csv')

# Change the name of the column 'pubChem' to 'pid' in Metlin dataset to be called the same as in the SMRT dataset
descriptors_Metlin.rename(columns={'pubChem': 'pid'}, inplace=True)

# Drop duplicates
descriptors_Metlin.drop_duplicates(subset='pid', keep='first', inplace=True)
descriptors_SMRT.drop_duplicates(subset='pid', keep='first', inplace=True)

# Keep only the common columns, so they can have an easy merge
common_columns = descriptors_Metlin.columns.intersection(descriptors_SMRT.columns)
descriptors_Metlin = descriptors_Metlin[common_columns]
descriptors_SMRT = descriptors_SMRT[common_columns]

# Print the shape of the updated DataFrames
print(descriptors_Metlin.shape)
print(descriptors_SMRT.shape)

# Merge both dataframes
merged_descriptors = pd.concat([descriptors_Metlin, descriptors_SMRT], ignore_index=True)
print(merged_descriptors.shape)

# Search for the same compounds in SMRT and Metlin
duplicated_pids = merged_descriptors['pid'].duplicated(keep=False)
common_compounds = merged_descriptors[duplicated_pids]
print(common_compounds.shape)

# Now, you can save 'common_compounds' to a single file:
# common_compounds.to_csv(f'./results/descriptors_common.csv', index=False)
# Or split it into chunks of 1024 columns each, so my LibreOffice can open it
num_columns = common_compounds.shape[1]
max_columns_per_chunk = 1024
num_chunks = (num_columns + max_columns_per_chunk - 1) // max_columns_per_chunk
for chunk_number in range(num_chunks):
    start_column = chunk_number * max_columns_per_chunk
    end_column = (chunk_number + 1) * max_columns_per_chunk
    chunk_df = common_compounds.iloc[:, start_column:end_column]
    chunk_df.to_csv(f'./results/descriptors_common_chunk_{chunk_number + 1}.csv', index=False)
