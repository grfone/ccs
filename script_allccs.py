import os
import pandas as pd

""" PART1 MERGE ALL FILES
# Directory containing all .csv files
directory = "./resources/allccs"

# Get a list of all .csv files in the directory
csv_files = [file for file in os.listdir(directory) if file.endswith(".csv")]

# Sort the list of .csv files alphabetically
csv_files.sort()

# Initialize an empty DataFrame to store merged data
merged_data = pd.DataFrame()

# Iterate over each .csv file and merge its data into the merged_data DataFrame

for file in csv_files:
    # Read the data from the current .csv file
    file_path = os.path.join(directory, file)
    data = pd.read_csv(file_path)

    # Merge the data into the merged_data DataFrame
    merged_data = pd.concat([merged_data, data], ignore_index=True)

# Write the merged data to a new .csv file
merged_file_path = "./results/AllCCS2.csv"
merged_data.to_csv(merged_file_path, index=False)

print("Merged data saved to:", merged_file_path)

"""

""" PART 2: EXTRACT THE EXPERIMENTAL VALUES TO A NEW FILE
 Load the CSV file into a DataFrame
allccs2_df = pd.read_csv("./results/AllCCS2.csv")

# Filter the DataFrame based on the "Type" column
allccs2_experimental_df = allccs2_df[allccs2_df["Type"] == "Experimental CCS"]

# Save the filtered DataFrame to a new CSV file
allccs2_experimental_df.to_csv("./results/AllCCS2_experimental.csv", index=False)

print("Filtered data saved to: ./resources/AllCCS2_experimental.csv")
"""

"""Parte COGER inchies
import requests


def get_inchi_from_smiles(smiles):    
    url = "http://www.chemspider.com/InChI.asmx/SMILESToInChI"
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    data = {'smiles': smiles}

    response = requests.post(url, headers=headers, data=data)

    if response.status_code == 200:

        inchi = response.text.replace('<?xml version="1.0" encoding="utf-8"?>', '').strip()
        inchi = inchi.replace('<string xmlns="http://www.chemspider.com/">', '').replace('</string>', '').strip()
        return inchi
    else:
        return None


# Load the CSV file into a DataFrame
allccs2_experimental_df = pd.read_csv("./results/AllCCS2_experimental.csv")

# Add a new column "InChI" with values obtained from applying get_inchi_from_smiles to "Structure" column
allccs2_experimental_df["InChI"] = allccs2_experimental_df["Structure"].apply(get_inchi_from_smiles)

# Save the DataFrame with the new column as a CSV file
allccs2_experimental_df.to_csv("./results/AllCCS2_experimental_with_inchis.csv", index=False)

print(allccs2_experimental_df.head())
"""

""" PARTE COMPARAR CCSs
# Load the Metlin.csv and AllCCS2_experimental_with_inchis.csv into DataFrames
metlin_df = pd.read_csv("./resources/Metlin.csv")
allccs2_df = pd.read_csv("./results/AllCCS2_experimental_with_inchis.csv")

# Do the averages and save them into CCS_AVG
metlin_df['CCS_AVG'] = metlin_df[['CCS1', 'CCS2', 'CCS3']].mean(axis=1).round(2)
# Select only the desired columns and rename "inchi" to "InChI" and "CCS_AVG" to "CCS-Metlin"
metlin_df = metlin_df[['inchi', 'Adduct', 'CCS_AVG']].rename(columns={'inchi': 'InChI', 'CCS_AVG': 'CCS-Metlin'})

# Select only the desired columns "InChI", "Adduct", and "CCS" and rename "CCS_AVG" to "CCS-Metlin"
allccs2_df = allccs2_df[['InChI', 'Adduct', 'CCS']].rename(columns={'CCS': 'CCS-AllCCS'})
# Remove the last character from each value in the "Adduct" column
allccs2_df['Adduct'] = allccs2_df['Adduct'].str[:-1]



# Perform an inner merge based on the columns "InChI" and "Adduct" for allccs2_df
coincidences_allccs2_df = pd.merge(allccs2_df, metlin_df, how='inner', on=['InChI', 'Adduct'])

# Calculate the absolute difference between the values in the "CCS-Metlin" and "CCS-AllCCS" columns
coincidences_allccs2_df['Abs_error'] = (abs(coincidences_allccs2_df['CCS-AllCCS'] - coincidences_allccs2_df['CCS-Metlin'])).round(2)

# Calculate the absolute difference between the values in the "CCS" and "CCS_AVG" columns
coincidences_allccs2_df['Rel_error(%)'] = (100 * coincidences_allccs2_df['Abs_error'] / coincidences_allccs2_df['CCS-Metlin']).round(2)

# Save the new DataFrame to a new CSV file
coincidences_allccs2_df.to_csv("./results/coincidences_df.csv", index=False)
print("Coincidences saved to: ./results/coincidences_between_metlin_and_allccs.csv")
"""

import requests


# Global variable to keep track of the number of times the function has been called
call_count = 0

def get_inchi_from_smiles(smiles):
    global call_count  # Access the global variable
    global total_length

    # Increment the call_count variable
    call_count += 1

    if call_count % 1000 == 0: print(100*call_count/total_length, "% done")

    url = "http://www.chemspider.com/InChI.asmx/SMILESToInChI"
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    data = {'smiles': smiles}

    response = requests.post(url, headers=headers, data=data)

    if response.status_code == 200:

        inchi = response.text.replace('<?xml version="1.0" encoding="utf-8"?>', '').strip()
        inchi = inchi.replace('<string xmlns="http://www.chemspider.com/">', '').replace('</string>', '').strip()
        return inchi
    else:
        return None


# Load the CSV file into a DataFrame
allccs2_experimental_df = pd.read_csv("./results/AllCCS2.csv")

total_length = len(allccs2_experimental_df)

# Add a new column "InChI" with values obtained from applying get_inchi_from_smiles to "Structure" column
allccs2_experimental_df["InChI"] = allccs2_experimental_df["Structure"].apply(get_inchi_from_smiles)

# Save the DataFrame with the new column as a CSV file
allccs2_experimental_df.to_csv("./results/AllCCS2_with_inchis.csv", index=False)

print(allccs2_experimental_df.head())
