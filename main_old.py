import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def compute_sverage_sd_and_plot_histogram(compounds, name):
    """
    :param compounds: dataframe with the compound data
    :param name: text string representative of the data to be processed
    :return: average value of the standard deviation of the three CSS measurements of each compound.
    """
    css_sd = compounds[['CCS1', 'CCS2', 'CCS3']].std(axis=1, skipna=True)
    plt.figure(figsize=(24, 8))
    plt.hist(css_sd, bins=300)
    plt.title("Standard deviation of each molecule: " + name, fontsize=24)
    plt.savefig("./results/Sd of the CCS" + name + ".png")
    return css_sd.mean()

def plot_css_histogram(compounds,name, row):
    css = compounds[[row]]
    plt.figure(figsize=(24, 8))
    plt.hist(css, bins=300)
    plt.title("Css: " + name , fontsize=24)
    plt.savefig("./results/CCS" + name + row+ ".png")



# Read the CSV file into a pandas DataFrame
metlin = pd.read_csv('resources/Metlin.csv')
normalize = False

# Dataframe correction
def custom_function(row):
    if row['Dimer.1'] == "Dimmer":
        if row['Adduct'] == "[M+H]":
            return (row['m/z']-1.00295)*2 + 1.00295
        elif row['Adduct'] == "[M-H]":
            return (row['m/z']+1.00295)*2 - 1.00295
        elif row['Adduct'] == "[M+Na]":
            return (row['m/z']-22.9892)*2 + 22.9892
    else:
        return row['m/z']


# Apply the custom function to each row
metlin['m/z'] = metlin.apply(custom_function, axis=1)
if normalize:
    for column in ['CCS1', 'CCS2', 'CCS3', 'CCS_AVG']:
        metlin[column] = metlin[column] / metlin['m/z']


print("CSS average value and standard deviations:")
print("All: average CSS:", metlin['CCS_AVG'].mean(),
      "Average standard deviation between runs of the same molecule, all compounds:",
      compute_sverage_sd_and_plot_histogram(metlin, "All"))

met_h = metlin[metlin["Precursor Adduct"].str.contains("M+H", regex=False)].copy(deep=True)
met_na = metlin[metlin["Precursor Adduct"].str.contains("M+Na", regex=False)].copy(deep=True)
met_h_neg = metlin[metlin["Precursor Adduct"].str.contains("M-H", regex=False)].copy(deep=True)
met_h_h_neg = pd.concat([met_h, met_h_neg], axis=0).copy(deep=True)
met_h_na = pd.concat([met_h, met_na], axis=0).copy(deep=True)

print("H+: average CSS:", met_h['CCS_AVG'].mean(), " standard deviation between runs of the same molecule",
      compute_sverage_sd_and_plot_histogram(met_h, "H+"))
print("H-: average CSS:", met_h_neg['CCS_AVG'].mean(), " standard deviation between runs of the same molecule, H-:",
      compute_sverage_sd_and_plot_histogram(met_h_neg, "H-"))
print("Na+: average CSS:", met_na['CCS_AVG'].mean(), "  standard deviation between runs of the same molecule, HNa:",
      compute_sverage_sd_and_plot_histogram(met_na, "Na+"))
print("H+ and H-: average CSS", met_h_h_neg['CCS_AVG'].mean(),
      " standard deviation between runs of the same molecule, H+ and Na+:",
      compute_sverage_sd_and_plot_histogram(met_h_h_neg, "H+ and H-"))

dfs={'metlin' : metlin,
     'met_h' : met_h,
     'met_na' : met_na,
     'met_h_neg' : met_h_neg}
colnames = ['CCS1', 'CCS2', 'CCS3', 'CCS_AVG']

for df_name, df in dfs.items():
    for colname in colnames:
        plot_css_histogram(df,df_name, colname)

css = metlin['CCS1']-metlin['CCS2']
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS1 -CSS2: ", fontsize=24)
plt.savefig("./results/CSS1 _ CSS2" + ".png")
css = css[css >=-1]
css = css[css <= 1]
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS1 -CSS2 0.5: ", fontsize=24)
plt.savefig("./results/CSS1 _ CSS2 1" + ".png")

css = metlin['CCS1']-metlin['CCS3']
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS1 -CSS3: ", fontsize=24)
plt.savefig("./results/CSS1 _ CSS3" + ".png")
css = css[css >=-1]
css = css[css <= 1]
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS1 -CSS3 0.5: ", fontsize=24)
plt.savefig("./results/CSS1 _ CSS3 1" + ".png")

css = metlin['CCS2']-metlin['CCS3']
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS2 -CSS3: ", fontsize=24)
plt.savefig("./results/CSS2 _ CSS3" + ".png")
css = css[css >=-1]
css = css[css <= 1]
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS2 -CSS3 0.5: ", fontsize=24)
plt.savefig("./results/CSS2 _ CSS3 1" + ".png")

css = metlin['CCS_AVG']-metlin['CCS1']
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS_AVG -CSS1: ", fontsize=24)
plt.savefig("./results/CSS_AVG _ CSS1" + ".png")
css = css[css >=-1]
css = css[css <= 1]
plt.figure(figsize=(24, 8))
plt.hist(css, bins=300)
plt.title("CSS_AVG -CSS1 0.5: ", fontsize=24)
plt.savefig("./results/CSS_AVG _ CSS1 1" + ".png")

print("\nAverage standard deviation between experiments of the same molecule with different adducts (H+, H-, Na+):",
      metlin.groupby('InChIKEY')['CCS_AVG'].std().mean())
print("Average standard deviation between experiments of the same molecule with different adducts (H+, H-):",
      met_h_h_neg.groupby('InChIKEY')['CCS_AVG'].std().mean())
print("Average standard deviation between experiments of the same molecule with different adducts (H+, Na+):",
      met_h_na.groupby('InChIKEY')['CCS_AVG'].std().mean())

met_h.rename(columns={'CCS_AVG': 'CCS_AVG_H+'}, inplace=True)
met_h_neg.rename(columns={'CCS_AVG': 'CCS_AVG_H-'}, inplace=True)
met_na.rename(columns={'CCS_AVG': 'CCS_AVG_Na+'}, inplace=True)

metH_Na = pd.merge(met_h, met_na, on="InChIKEY")
metH_Hne = pd.merge(met_h, met_h_neg, on="InChIKEY")

plt.figure()
plt.scatter(metH_Hne["CCS_AVG_H+"], metH_Hne["CCS_AVG_H-"])
plt.title("H+ vs H-")
plt.xlabel("CCS_AVG_H+")
plt.ylabel("CCS_AVG_H-")
plt.savefig("./results/H+ vs H-.png")

plt.figure()
plt.scatter(metH_Na["CCS_AVG_H+"], metH_Na["CCS_AVG_Na+"])
plt.title("H+ vs Na+")
plt.xlabel("CCS_AVG_H+")
plt.ylabel("CCS_AVG_Na+")
plt.savefig("./results/H+ vs Na+.png")

plt.figure()
metH_Na_gt_12 = metH_Na[metH_Na["CCS_AVG_H+"] > 1.2 * metH_Na["CCS_AVG_Na+"]]
plt.scatter(metH_Na_gt_12["CCS_AVG_H+"], metH_Na_gt_12["CCS_AVG_Na+"])
plt.title("H+ vs Na+, H+ > 1.2 Na")
plt.xlabel("CCS_AVG_H+")
plt.ylabel("CCS_AVG_Na+")
plt.savefig("./results/H+ vs Na+, H+ > 1.2 Na.png")

metH_Hne_gt_12 = metH_Hne[metH_Hne["CCS_AVG_H+"] > 1.2 * metH_Hne["CCS_AVG_H-"]]
plt.scatter(metH_Hne_gt_12["CCS_AVG_H+"], metH_Hne_gt_12["CCS_AVG_H-"])
plt.title("H+ vs H-, H+ > 1.2 H-")
plt.xlabel("CCS_AVG_H+")
plt.ylabel("CCS_AVG_H-")
plt.savefig("./results/H+ vs H-, H+ > 1.2 H-.png")

plt.figure()
metH_Na_gt_08 = metH_Na[metH_Na["CCS_AVG_H+"] < 0.8 * metH_Na["CCS_AVG_Na+"]]
plt.scatter(metH_Na_gt_08["CCS_AVG_H+"], metH_Na_gt_08["CCS_AVG_Na+"])
plt.title("H+ vs Na+, H+ > 0.8 Na")
plt.xlabel("CCS_AVG_H+")
plt.ylabel("CCS_AVG_Na+")
plt.savefig("./results/H+ vs Na+, H+ > 0.8 Na.png")

plt.figure()
metH_Hne_gt_08 = metH_Hne[metH_Hne["CCS_AVG_H+"] < 0.8 * metH_Hne["CCS_AVG_H-"]]
plt.scatter(metH_Hne_gt_08["CCS_AVG_H+"], metH_Hne_gt_08["CCS_AVG_H-"])
plt.title("H+ vs H-, H+ > 0.8 H-")
plt.xlabel("CCS_AVG_H+")
plt.ylabel("CCS_AVG_H-")
plt.savefig("./results/H+ vs H-, H+ > 0.8 H-.png")

#code to check CCS_AVG

x = metlin['CCS_AVG'] - metlin[['CCS1', 'CCS2', 'CCS3']].mean(axis=1)
error = [x.abs() > 0.01]
print("Diferences of more than 0.01: ", np.sum(error))
error = [x.abs() > 0.1]
print("Diferences of more than 0.1: ", np.sum(error))
error = [x.abs() > 1]
print("Diferences of more than 1: ", np.sum(error))