import matplotlib.pyplot as plt
import pandas as pd


# Read the CSV file into a pandas DataFrame
metlin = pd.read_csv('resources/Metlin.csv')


def plots_before_correction():
    # Create a mask based on the values in 'Dimer.1'
    # ----------------------------------------------
    mask_monomer = metlin['Dimer.1'] == 'Monomer'

    plt.figure()

    # Plotting
    plt.scatter(metlin.loc[mask_monomer, 'm/z'], metlin.loc[mask_monomer, 'CCS_AVG'], c='blue', s=20, label='Monomer')
    plt.scatter(metlin.loc[~mask_monomer, 'm/z'], metlin.loc[~mask_monomer, 'CCS_AVG'], c='red', s=20, label='Dimer')

    # Add title, labels and legend
    plt.title("CCS_AVG vs m/z(before correction)")
    plt.xlabel('m/z')
    plt.ylabel('CCS_AVG')
    plt.legend()

    # Save the plot
    plt.savefig("./results/Monomer vs Dimer.png")

    # Create a mask based on the values in 'Dimer.1' and 'Adduct'
    # -----------------------------------------------------------
    mask_monomer = metlin['Dimer.1'] == 'Monomer'
    mask_h_minus = metlin['Adduct'] == '[M-H]'
    mask_h_plus = metlin['Adduct'] == '[M+H]'
    mask_na_plus = metlin['Adduct'] == '[M+Na]'

    # Plotting H- adduct
    plt.figure()
    plt.scatter(metlin.loc[mask_h_minus, 'm/z'], metlin.loc[mask_h_minus, 'CCS_AVG'], s=5)
    plt.title("CCS_AVG vs m/z, H- adduct (before correction)")
    plt.xlabel('m/z')
    plt.ylabel('CCS_AVG')
    plt.savefig("./results/H- adduct.png")

    # Plotting H+ adducts
    plt.figure()
    plt.scatter(metlin.loc[mask_h_plus, 'm/z'], metlin.loc[mask_h_plus, 'CCS_AVG'], s=5)
    plt.title("CCS_AVG vs m/z, H+ adduct (before correction)")
    plt.xlabel('m/z')
    plt.ylabel('CCS_AVG')
    plt.savefig("./results/H+ adduct.png")

    # Plotting Na+ adducts
    plt.figure()
    plt.scatter(metlin.loc[mask_na_plus, 'm/z'], metlin.loc[mask_na_plus, 'CCS_AVG'], s=5)
    plt.title("CCS_AVG vs m/z, Na+ adduct (before correction)")
    plt.xlabel('m/z')
    plt.ylabel('CCS_AVG')
    plt.savefig("./results/Na+ adduct.png")


    # Plotting all adducts
    light_blue, mask_monomer_h_minus = '#ADD8E6', (mask_monomer & mask_h_minus)
    medium_blue, mask_monomer_h_plus, = '#6495ED', (mask_monomer & mask_h_plus)
    dark_blue, mask_monomer_na_plus = '#4169E1', (mask_monomer & mask_na_plus)
    light_red, mask_dimer_h_minus = '#FFC0CB', (~mask_monomer & mask_h_minus)
    medium_red, mask_dimer_h_plus = '#FF1493', (~mask_monomer & mask_h_plus)
    dark_red, mask_dimer_na_plus = '#B22222', (~mask_monomer & mask_na_plus)

    plt.figure()

    plt.scatter(metlin.loc[mask_monomer_h_minus, 'm/z'], metlin.loc[mask_monomer_h_minus, 'CCS_AVG'], c=light_blue, s=5,
                label='H- monomers')
    plt.scatter(metlin.loc[mask_monomer_h_plus, 'm/z'], metlin.loc[mask_monomer_h_plus, 'CCS_AVG'], c=medium_blue, s=5,
                label='H+ monomers')
    plt.scatter(metlin.loc[mask_monomer_na_plus, 'm/z'], metlin.loc[mask_monomer_na_plus, 'CCS_AVG'], c=dark_blue, s=5,
                label='Na+ monomers')
    plt.scatter(metlin.loc[mask_dimer_h_minus, 'm/z'], metlin.loc[mask_dimer_h_minus, 'CCS_AVG'], c=light_red, s=5,
                label='H- dimes')
    plt.scatter(metlin.loc[mask_dimer_h_plus, 'm/z'], metlin.loc[mask_dimer_h_plus, 'CCS_AVG'], c=medium_red, s=5,
                label='H+ dimers')
    plt.scatter(metlin.loc[mask_dimer_na_plus, 'm/z'], metlin.loc[mask_dimer_na_plus, 'CCS_AVG'], c=dark_red, s=5,
                label='Na+ dimers')

    # Add title, labels and legend
    plt.title("CCS_AVG vs m/z(before correction)")
    plt.xlabel('m/z')
    plt.ylabel('CCS_AVG')
    plt.legend()

    # Plot the plot
    plt.savefig("./results/Monomer vs Dimer (different adducts).png")


plots_before_correction()


def compute_average_sd_and_plot_histogram(compounds, name):
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

for column in ['CCS1', 'CCS2', 'CCS3', 'CCS_AVG']:
    metlin[column] = metlin[column] / metlin['m/z']


print("CSS average value and standard deviations:")
print("All: average CSS:", metlin['CCS_AVG'].mean(),
      "Average standard deviation between runs of the same molecule, all compounds:",
      compute_average_sd_and_plot_histogram(metlin, "All"))

met_h = metlin[metlin["Precursor Adduct"].str.contains("M+H", regex=False)].copy(deep=True)
met_na = metlin[metlin["Precursor Adduct"].str.contains("M+Na", regex=False)].copy(deep=True)
met_h_neg = metlin[metlin["Precursor Adduct"].str.contains("M-H", regex=False)].copy(deep=True)
met_h_h_neg = pd.concat([met_h, met_h_neg], axis=0).copy(deep=True)
met_h_na = pd.concat([met_h, met_na], axis=0).copy(deep=True)

print("H+: average CSS:", met_h['CCS_AVG'].mean(), " standard deviation between runs of the same molecule",
      compute_average_sd_and_plot_histogram(met_h, "H+"))
print("H-: average CSS:", met_h_neg['CCS_AVG'].mean(), " standard deviation between runs of the same molecule, H-:",
      compute_average_sd_and_plot_histogram(met_h_neg, "H-"))
print("Na+: average CSS:", met_na['CCS_AVG'].mean(), "  standard deviation between runs of the same molecule, HNa:",
      compute_average_sd_and_plot_histogram(met_na, "Na+"))
print("H+ and H-: average CSS", met_h_h_neg['CCS_AVG'].mean(),
      " standard deviation between runs of the same molecule, H+ and Na+:",
      compute_average_sd_and_plot_histogram(met_h_h_neg, "H+ and H-"))


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

plt.figure()
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