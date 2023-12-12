import pandas

# Read the CSV file into a pandas DataFrame
metlin = pandas.read_csv('mMetlin.csv')

print("Average standard deviation between runs of the same molecule:")
print(metlin[['CCS1', 'CCS2', 'CCS3']].std(axis=1, skipna=True).mean())

print("Average standard deviation between experiments of the same molecule with different adduct:")
print(metlin.groupby('InChIKEY')['CCS_AVG'].std().mean())
