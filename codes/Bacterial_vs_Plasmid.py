"""
Created on 26/04/2023

Author: Lucy Androsiuk
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import fastcluster
from sklearn import preprocessing

# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

# working directories
out_dir = f"{path}/data_calculations"
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"

# working files
all_taxa = f'{path}/Article/Haroon/RedSea_TaxaCounts_AllTaxa.csv'
all_genus = f'{path}/Article/Haroon/RedSea_TaxaRelAbund_Genus.csv'
all_species = f'{path}/Article/Haroon/RedSea_TaxaRelAbund_Species.csv'
all_plasmids = r"../res/all_cov.csv"
sample_names = r'../res/stations.txt'

#### all taxa table
df_taxa = pd.read_csv(all_taxa,header = 0)
df_taxa['St_Depth'] = df_taxa['Station'].astype(str).str.cat(df_taxa['Depth (m)'].astype(str), sep='_')

#### all genus table
df_genus = pd.read_csv(all_genus,header = 0)
df_genus['St_Depth'] = df_genus['Station'].astype(str).str.cat(df_genus['Depth (m)'].astype(str), sep='_')

#### all species table
df_species = pd.read_csv(all_species,header = 0)
df_species['St_Depth'] = df_species['Station'].astype(str).str.cat(df_species['Depth (m)'].astype(str), sep='_')

#### all plasmids table
df_plasmids = pd.read_csv(all_plasmids, header = 0, index_col = 0)
df_plasmids.index = df_plasmids.index.str.replace(r'_l\w+', '', regex = True)
df_plasmids = df_plasmids.drop('94_LNODE_1', axis = 0)
df_plasmids = df_plasmids.transpose()
df_plasmids = df_plasmids.reset_index().rename(columns={'index':'Sample'})
# Define the columns to exclude from division
exclude_columns = ['Sample']
df_plasmids.loc[:, ~df_plasmids.columns.isin(exclude_columns)] = df_plasmids.loc[:, ~df_plasmids.columns.isin(exclude_columns)] / 100
# Extract only the numeric part of the column
df_plasmids['Sample'] = df_plasmids['Sample'].str.extract('(\d+)').astype(int)
df_stat = pd.read_table(sample_names)
merged = df_plasmids.merge(df_stat, on = 'Sample')
merged = merged.drop('Sample', axis = 1)
#print(merged)

def corr_coef(plasmid_df ,bact_df, name):
    #bullshit
    # create a new dataframe with only the plasmid and bacterium columns
    plasmid_bacteria_df = pd.concat([plasmid_df.iloc[:,:-1], bact_df.iloc[:, 17:-1]], axis=1)
    # calculate the correlation coefficients between all pairs of plasmids and bacteria
    corr_matrix = plasmid_bacteria_df.corr()
    print(corr_matrix)
    bacteria_to_remove = list(bact_df.columns.values)[17:-1]
    print(bacteria_to_remove)
    plasmids_to_remove = list(plasmid_df.columns.values)[:-1]
    print(plasmids_to_remove)
    cols_to_remove = [b for b in bacteria_to_remove if b in corr_matrix.columns]
    rows_to_remove = [p for p in plasmids_to_remove if p in corr_matrix.index]
    # subset the correlation matrix to exclude the rows and columns to remove
    corr_matrix_subset = corr_matrix.drop(rows_to_remove, axis = 0).drop(cols_to_remove, axis = 1)

    print(corr_matrix_subset)
    figure = sns.heatmap(corr_matrix_subset, cmap = 'viridis', vmin = -1, vmax = 1, annot = False,
                     cbar_kws = {'label': 'Pearson correlation', "ticks": [-1,0, 1]})
    figure = figure.get_figure()
    svg_name = 'heatmap_corr_' + name + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'heatmap_corr_' + name + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        figure.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        figure.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')

    figure2 = sns.clustermap(data = corr_matrix_subset,
                             metric = "euclidean",
                             method = 'ward',
                             cmap = "coolwarm",vmin = -1, vmax = 1, center=0,
                             xticklabels = False,
                             yticklabels = False,
                             cbar_pos=(0, .3, .03, .3),
                             cbar_kws = {'label': 'Pearson correlation', "ticks": [-1,0, 1]})

    figure2.ax_col_dendrogram.remove()
    figure2.ax_row_dendrogram.remove()
    plt.show()



def contingency(plasmid_df ,bact_df, name):
    # Step 1: Create a new DataFrame containing only the plasmid and bacterium columns
    plasmid_df = plasmid_df.iloc[:,:-1].applymap(lambda x: 1 if x > 0.99 else 0)
    print(bact_df.iloc[:, 17:-2])
    x = bact_df.iloc[:, 17:-2].values  # returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled, columns = bact_df.iloc[:, 17:-2].columns)
    print(df)
    bact_df = df.applymap(lambda x: 0 if x < 0.99 else 1)
    plasmid_bacteria_df = pd.concat([plasmid_df.iloc[:,:-1], bact_df], axis = 1)
    print(plasmid_bacteria_df)
    # Step 2: Create new columns in the DataFrame representing the presence/absence of each bacterium
    for bacterium in plasmid_bacteria_df.columns[len(plasmid_df.columns) - 1:]:
        plasmid_bacteria_df[bacterium + '_presence'] = plasmid_bacteria_df[bacterium].fillna(0)
    print(plasmid_bacteria_df)
    # Step 3: Create a new DataFrame to store the contingency table
    contingency_table = pd.crosstab(index = bact_df.columns,
                                     columns = plasmid_df[:-1].columns)
    print(contingency_table)
    # Step 4: Calculate the counts for each cell in the contingency table
    for bacterium in plasmid_bacteria_df.columns[len(plasmid_df.columns):]:
        contingency_table[bacterium]['plasmid_absent'] = (
                    (plasmid_bacteria_df[plasmid_bacteria_df[bacterium] == 0]).iloc[:,
                    len(plasmid_df.columns):] == 1).sum()
        contingency_table[bacterium]['plasmid_present'] = (
                    (plasmid_bacteria_df[plasmid_bacteria_df[bacterium] == 1]).iloc[:,
                    len(plasmid_df.columns):] == 1).sum()
    print(contingency_table)
    # Step 5: Calculate the odds ratio for each plasmid-bacterium pair
    odds_ratio_table = contingency_table.copy()
    for bacterium in plasmid_bacteria_df.columns[len(plasmid_df.columns):]:
        total_present = contingency_table[bacterium].sum()
        total_absent = contingency_table[bacterium]['plasmid_absent'] + contingency_table[bacterium][
            'plasmid_present'] - total_present
        odds_ratio_table[bacterium]['plasmid_absent'] = contingency_table[bacterium]['plasmid_absent'] / total_absent
        odds_ratio_table[bacterium]['plasmid_present'] = contingency_table[bacterium]['plasmid_present'] / total_present
        odds_ratio_table[bacterium] = odds_ratio_table[bacterium]['plasmid_present'] / odds_ratio_table[bacterium][
            'plasmid_absent']
    print(odds_ratio_table)
#corr_coef(merged, df_taxa, 'taxa')
#corr_coef(merged, df_species, 'species')
contingency(merged, df_genus, 'genus')

