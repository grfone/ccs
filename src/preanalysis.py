import os
import matplotlib.pyplot as plt
import pandas as pd
import gdown

from utils.from_copied_link_to_download_link import transform_link


class Preanalysis:
    def __init__(self, normalization_flag):
        """
        Initialize the Preanalysis object.

        Parameters:
            normalization_flag (bool): A flag indicating whether normalization should be performed.
        """
        self.metlin_df = self._load_data()
        self._remove_excess_columns_and_rows()
        self._correct_mz_column()
        self.normalize_ccs_columns() if normalization_flag else None
        os.mkdir('results') if not os.path.exists('./results') else None
        os.mkdir('results/preanalysis_graphs') if not os.path.exists('./results/preanalysis_graphs') else None
        self.do_first_set_of_plots(normalization_flag)
        self.do_second_set_of_plots(normalization_flag)

    @staticmethod
    def _load_data():
        """
        Load the 'Metlin' database.

        Returns:
            pd.DataFrame: The loaded Metlin database as a Pandas DataFrame.
        """
        # Create resources folder if it doesn't exist yet
        if not os.path.exists('./resources'):
            os.mkdir('resources')

        # Download 'Metlin' database
        if not os.path.exists('./resources/Metlin.csv'):
            link_2_metlin_file = 'https://drive.google.com/file/d/1u4g6bvGnVxqM4p0XJs867VeAF4RlYeD3/view?usp=drive_link'
            gdown.download(transform_link(link_2_metlin_file), './resources/Metlin.csv')

        # Read metlin database from resources
        metlin_df = pd.read_csv('resources/Metlin.csv')

        return metlin_df

    def _remove_excess_columns_and_rows(self):
        """
        Remove excess columns and rows from the Metlin DataFrame.
        """
        self.metlin_df.drop(columns=['dimer line', 'CCS', 'm/z.2'], inplace=True)
        self.metlin_df.drop(self.metlin_df.index[-5:], inplace=True)
        self.metlin_df.drop(self.metlin_df[self.metlin_df['Molecule Name'].str.contains("Tm_")].index, inplace=True)

    def _correct_mz_column(self):
        """
        Correct the 'm/z' column in the Metlin DataFrame based on the real weight of the dimer.
        """
        # Dataframe correction function
        def mz_correction_function(row):
            if row['Dimer.1'] == "Dimmer":
                if row['Adduct'] == "[M+H]":
                    return (row['m/z'] - 1.00295) * 2 + 1.00295
                elif row['Adduct'] == "[M-H]":
                    return (row['m/z'] + 1.00295) * 2 - 1.00295
                elif row['Adduct'] == "[M+Na]":
                    return (row['m/z'] - 22.9892) * 2 + 22.9892
            else:
                return row['m/z']

        # Apply the correction function the m/z of all dimers because their m/z is 2M, not just M
        self.metlin_df['m/z'] = self.metlin_df.apply(mz_correction_function, axis=1)

    def normalize_ccs_columns(self):
        """
        Normalize the 'CCS' columns in the Metlin DataFrame if the normalization flag is True.
        """
        # Apply normalization to each row of the columns 'CCS1', 'CCS2', 'CCS3' and 'CCS_AVG if the normalization flag
        # has been used
        for column in ['CCS1', 'CCS2', 'CCS3', 'CCS_AVG']:
            self.metlin_df[column] = self.metlin_df[column] / self.metlin_df['m/z']

    def do_first_set_of_plots(self, normalization_flag):
        """
        Generate and save the first set of plots based on CCS_AVG and m/z values.

        Parameters:
            normalization_flag (bool): A flag indicating whether normalization has been performed.
        """
        # Create suffix for the plot files if normalization flag is True
        suffix = "_normalized" if normalization_flag else ""

        # Create a mask based on the values in 'Dimer.1'
        # ----------------------------------------------
        mask_monomer = self.metlin_df['Dimer.1'] == 'Monomer'

        plt.figure()

        # Plotting
        plt.scatter(self.metlin_df.loc[mask_monomer, 'm/z'], self.metlin_df.loc[mask_monomer, 'CCS_AVG'], c='blue',
                    s=20, label='Monomer')
        plt.scatter(self.metlin_df.loc[~mask_monomer, 'm/z'], self.metlin_df.loc[~mask_monomer, 'CCS_AVG'], c='red',
                    s=20, label='Dimer')

        # Add title, labels and legend
        plt.title("CCS_AVG vs m/z(before correction)")
        plt.xlabel('m/z')
        plt.ylabel('CCS_AVG')
        plt.legend()

        # Save the plot
        plt.savefig(f"./results/preanalysis_graphs/Monomer vs Dimer{suffix}.png")

        # Create a mask based on the values in 'Dimer.1' and 'Adduct'
        # -----------------------------------------------------------
        mask_monomer = self.metlin_df['Dimer.1'] == 'Monomer'
        mask_h_minus = self.metlin_df['Adduct'] == '[M-H]'
        mask_h_plus = self.metlin_df['Adduct'] == '[M+H]'
        mask_na_plus = self.metlin_df['Adduct'] == '[M+Na]'

        # Plotting H- adduct
        plt.figure()
        plt.scatter(self.metlin_df.loc[mask_h_minus, 'm/z'], self.metlin_df.loc[mask_h_minus, 'CCS_AVG'], s=5)
        plt.title("CCS_AVG vs m/z, H- adduct (before correction)")
        plt.xlabel('m/z')
        plt.ylabel('CCS_AVG')
        plt.savefig(f"./results/preanalysis_graphs/H- adduct{suffix}.png")

        # Plotting H+ adducts
        plt.figure()
        plt.scatter(self.metlin_df.loc[mask_h_plus, 'm/z'], self.metlin_df.loc[mask_h_plus, 'CCS_AVG'], s=5)
        plt.title("CCS_AVG vs m/z, H+ adduct (before correction)")
        plt.xlabel('m/z')
        plt.ylabel('CCS_AVG')
        plt.savefig(f"./results/preanalysis_graphs/H+ adduct{suffix}.png")

        # Plotting Na+ adducts
        plt.figure()
        plt.scatter(self.metlin_df.loc[mask_na_plus, 'm/z'], self.metlin_df.loc[mask_na_plus, 'CCS_AVG'], s=5)
        plt.title("CCS_AVG vs m/z, Na+ adduct (before correction)")
        plt.xlabel('m/z')
        plt.ylabel('CCS_AVG')
        plt.savefig(f"./results/preanalysis_graphs/Na+ adduct{suffix}.png")

        # Plotting all adducts
        light_blue, mask_monomer_h_minus = '#ADD8E6', (mask_monomer & mask_h_minus)
        medium_blue, mask_monomer_h_plus, = '#6495ED', (mask_monomer & mask_h_plus)
        dark_blue, mask_monomer_na_plus = '#4169E1', (mask_monomer & mask_na_plus)
        light_red, mask_dimer_h_minus = '#FFC0CB', (~mask_monomer & mask_h_minus)
        medium_red, mask_dimer_h_plus = '#FF1493', (~mask_monomer & mask_h_plus)
        dark_red, mask_dimer_na_plus = '#B22222', (~mask_monomer & mask_na_plus)

        plt.figure()

        plt.scatter(self.metlin_df.loc[mask_monomer_h_minus, 'm/z'],
                    self.metlin_df.loc[mask_monomer_h_minus, 'CCS_AVG'], c=light_blue, s=5, label='H- monomers')
        plt.scatter(self.metlin_df.loc[mask_monomer_h_plus, 'm/z'],
                    self.metlin_df.loc[mask_monomer_h_plus, 'CCS_AVG'], c=medium_blue, s=5, label='H+ monomers')
        plt.scatter(self.metlin_df.loc[mask_monomer_na_plus, 'm/z'],
                    self.metlin_df.loc[mask_monomer_na_plus, 'CCS_AVG'], c=dark_blue, s=5, label='Na+ monomers')
        plt.scatter(self.metlin_df.loc[mask_dimer_h_minus, 'm/z'],
                    self.metlin_df.loc[mask_dimer_h_minus, 'CCS_AVG'], c=light_red, s=5, label='H- dimes')
        plt.scatter(self.metlin_df.loc[mask_dimer_h_plus, 'm/z'],
                    self.metlin_df.loc[mask_dimer_h_plus, 'CCS_AVG'], c=medium_red, s=5, label='H+ dimers')
        plt.scatter(self.metlin_df.loc[mask_dimer_na_plus, 'm/z'],
                    self.metlin_df.loc[mask_dimer_na_plus, 'CCS_AVG'], c=dark_red, s=5, label='Na+ dimers')

        # Add title, labels and legend
        plt.title("CCS_AVG vs m/z(before correction)")
        plt.xlabel('m/z')
        plt.ylabel('CCS_AVG')
        plt.legend()

        # Plot the plot
        plt.savefig(f"./results/preanalysis_graphs/Monomer vs Dimer (different adducts){suffix}.png")

    def do_second_set_of_plots(self, normalization_flag):
        """
        Generate and save the second set of plots based on CCS_AVG and m/z values.

        Parameters:
            normalization_flag (bool): A flag indicating whether normalization has been performed.
        """
        # Create suffix for the plot files if normalization flag is True
        suffix = "_normalized" if normalization_flag else ""

        def compute_average_sd_and_plot_histogram(compounds, name):
            """
            Compute the average value of the standard deviation of the three CCS measurements of each compound
            and plot a histogram.

            Parameters:
                compounds (pd.DataFrame): DataFrame with the compound data.
                name (str): Text string representative of the data to be processed.

            Returns:
                float: Average value of the standard deviation of the three CCS measurements of each compound.
            """
            css_sd = compounds[['CCS1', 'CCS2', 'CCS3']].std(axis=1, skipna=True)
            plt.figure(figsize=(24, 8))
            plt.hist(css_sd, bins=300)
            plt.title("Standard deviation of each molecule: " + name, fontsize=24)
            plt.savefig(f"./results/preanalysis_graphs/Sd of the CCS {name}{suffix}.png")
            return css_sd.mean()

        print("CSS average value and standard deviations:")
        print("All: average CSS:", self.metlin_df['CCS_AVG'].mean(),
              "Average standard deviation between runs of the same molecule, all compounds:",
              compute_average_sd_and_plot_histogram(self.metlin_df, "All"))

        met_h = self.metlin_df[self.metlin_df["Precursor Adduct"].str.contains("M+H", regex=False)].copy(deep=True)
        met_na = self.metlin_df[self.metlin_df["Precursor Adduct"].str.contains("M+Na", regex=False)].copy(deep=True)
        met_h_neg = self.metlin_df[self.metlin_df["Precursor Adduct"].str.contains("M-H", regex=False)].copy(deep=True)
        met_h_h_neg = pd.concat([met_h, met_h_neg], axis=0).copy(deep=True)
        met_h_na = pd.concat([met_h, met_na], axis=0).copy(deep=True)

        print("H+: average CSS:", met_h['CCS_AVG'].mean(), " standard deviation between runs of the same molecule",
              compute_average_sd_and_plot_histogram(met_h, "H+"))
        print("H-: average CSS:", met_h_neg['CCS_AVG'].mean(),
              " standard deviation between runs of the same molecule, H-:",
              compute_average_sd_and_plot_histogram(met_h_neg, "H-"))
        print("Na+: average CSS:", met_na['CCS_AVG'].mean(),
              "  standard deviation between runs of the same molecule, HNa:",
              compute_average_sd_and_plot_histogram(met_na, "Na+"))
        print("H+ and H-: average CSS", met_h_h_neg['CCS_AVG'].mean(),
              " standard deviation between runs of the same molecule, H+ and Na+:",
              compute_average_sd_and_plot_histogram(met_h_h_neg, "H+ and H-"))

        print(
            "\nAverage standard deviation between experiments of the same molecule with different adducts (H+, H-, Na+):",
            self.metlin_df.groupby('InChIKEY')['CCS_AVG'].std().mean())
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
        plt.savefig(f"./results/preanalysis_graphs/H+ vs H-{suffix}.png")

        plt.figure()
        plt.scatter(metH_Na["CCS_AVG_H+"], metH_Na["CCS_AVG_Na+"])
        plt.title("H+ vs Na+")
        plt.xlabel("CCS_AVG_H+")
        plt.ylabel("CCS_AVG_Na+")
        plt.savefig(f"./results/preanalysis_graphs/H+ vs Na+{suffix}.png")

        plt.figure()
        metH_Na_gt_12 = metH_Na[metH_Na["CCS_AVG_H+"] > 1.2 * metH_Na["CCS_AVG_Na+"]]
        plt.scatter(metH_Na_gt_12["CCS_AVG_H+"], metH_Na_gt_12["CCS_AVG_Na+"])
        plt.title("H+ vs Na+, H+ > 1.2 Na")
        plt.xlabel("CCS_AVG_H+")
        plt.ylabel("CCS_AVG_Na+")
        plt.savefig(f"./results/preanalysis_graphs/H+ vs Na+, H+ > 1.2 Na{suffix}.png")

        plt.figure()
        metH_Hne_gt_12 = metH_Hne[metH_Hne["CCS_AVG_H+"] > 1.2 * metH_Hne["CCS_AVG_H-"]]
        plt.scatter(metH_Hne_gt_12["CCS_AVG_H+"], metH_Hne_gt_12["CCS_AVG_H-"])
        plt.title("H+ vs H-, H+ > 1.2 H-")
        plt.xlabel("CCS_AVG_H+")
        plt.ylabel("CCS_AVG_H-")
        plt.savefig(f"./results/preanalysis_graphs/H+ vs H-, H+ > 1.2 H-{suffix}.png")

        plt.figure()
        metH_Na_gt_08 = metH_Na[metH_Na["CCS_AVG_H+"] < 0.8 * metH_Na["CCS_AVG_Na+"]]
        plt.scatter(metH_Na_gt_08["CCS_AVG_H+"], metH_Na_gt_08["CCS_AVG_Na+"])
        plt.title("H+ vs Na+, H+ > 0.8 Na")
        plt.xlabel("CCS_AVG_H+")
        plt.ylabel("CCS_AVG_Na+")
        plt.savefig(f"./results/preanalysis_graphs/H+ vs Na+, H+ > 0.8 Na{suffix}.png")

        plt.figure()
        metH_Hne_gt_08 = metH_Hne[metH_Hne["CCS_AVG_H+"] < 0.8 * metH_Hne["CCS_AVG_H-"]]
        plt.scatter(metH_Hne_gt_08["CCS_AVG_H+"], metH_Hne_gt_08["CCS_AVG_H-"])
        plt.title("H+ vs H-, H+ > 0.8 H-")
        plt.xlabel("CCS_AVG_H+")
        plt.ylabel("CCS_AVG_H-")
        plt.savefig(f"./results/preanalysis_graphs/H+ vs H-, H+ > 0.8 H-{suffix}.png")
