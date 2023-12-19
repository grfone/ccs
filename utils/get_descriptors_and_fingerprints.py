import numpy
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from collections import defaultdict
import pandas as pd


def get_them_from_smiles_in_metlin(metlin_df):
    """
    Extract molecular information from SMILES strings in the given Metlin DataFrame.

    Parameters:
    - metlin_df (pd.DataFrame): DataFrame containing Metlin data, including the 'smiles' column.
    """
    # Keep only the important columns
    metlin_df = metlin_df[['smiles', 'Adduct', 'm/z', '% CV', 'Dimer', 'Dimer.1', 'CCS_AVG']]
    # Remove molecules with no smiles (there aren't any in metlin)
    metlin_df = metlin_df[~metlin_df['smiles'].str.contains("--")]
    metlin_df = metlin_df.dropna(subset=['smiles'])
    # Create a column with all the 'mol' from the 'smiles'
    metlin_df['mol'] = metlin_df['smiles'].apply(Chem.MolFromSmiles)
    # Drop rows that rdkit can't be "kekulized" by rdkit (there are 7 rows in total)
    metlin_df = metlin_df.dropna(subset=['mol'])

    # Generate the descriptors
    list_of_descriptor_dictionaries = [Chem.Descriptors.CalcMolDescriptors(mol) for mol in metlin_df['mol']]
    # Create a defaultdict object in order to store lists of values for each key
    descriptors_dictionary = defaultdict(list)
    # Loop through each dictionary and append values to the corresponding key
    for dictionary in list_of_descriptor_dictionaries:
        for key, value in dictionary.items():
            descriptors_dictionary[key].append(value)
    # Convert the defaultdict to a regular dictionary
    descriptors_dictionary = dict(descriptors_dictionary)
    # Convert the dictionary to a DataFrame
    descriptors_df = pd.DataFrame.from_dict(descriptors_dictionary)

    # Generate fingerprints
    fingerprints_list = [list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)) for mol in metlin_df['mol']]
    # Create the name of the columns
    my_fingerprint_keys = ['Bit_{}'.format(i+1) for i in range(len(fingerprints_list[0]))]
    # Create the DataFrame
    fingerprints_df = pd.DataFrame(fingerprints_list, columns=my_fingerprint_keys)

    # Resetting the index of each DataFrame (so no 'new' rows appear when merging)
    metlin_df = metlin_df.reset_index(drop=True)
    descriptors_df = descriptors_df.reset_index(drop=True)
    fingerprints_df = fingerprints_df.reset_index(drop=True)
    # Merge metlin_df with descriptors_df and fingerprints_df
    metlin_df = pd.concat([metlin_df, descriptors_df, fingerprints_df], axis=1)
    # Save metlin_df
    metlin_df.to_csv('./results/metlin_df.csv', index=False)
