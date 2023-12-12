import pandas

# Read the CSV file into a pandas DataFrame
metlin = pandas.read_csv('mMetlin.csv')

# Remove all rows with the molecule name 'Tm_XXX' where XXX are three numbers
metlin_without_Tms = metlin[~metlin['Molecule Name'].str.contains("Tm_")]

print("Average standard deviation between runs of the same molecule:")
print(metlin_without_Tms[['CCS1', 'CCS2', 'CCS3']].std(axis=1, skipna=True).mean())

print("Average standard deviation between experiments of the same molecule with different adduct:")
print(metlin_without_Tms.groupby('InChIKEY')['CCS_AVG'].std().mean())
