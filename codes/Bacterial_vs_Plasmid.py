"""
Created on 26/04/2023

Author: Lucy Androsiuk
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import *

all_taxa = r'/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/Article/Haroon/RedSea_TaxaCounts_AllTaxa.csv'
all_genus = r'/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/Article/Haroon/RedSea_TaxaRelAbund_Genus.csv'
all_species = r'/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/Article/Haroon/RedSea_TaxaRelAbund_Species.csv'
all_plasmids = r"../res/all_cov.csv"
sample_names = r'../res/stations.txt'
visuals = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/visualisations"

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
    figure = sns.heatmap(corr_matrix, cmap = 'viridis', vmin = -1, vmax = 1, annot = False,
                     cbar_kws = {'label': 'Pearson correlation', "ticks": [-1,0, 1]})
    figure = figure.get_figure()
    svg_name = 'heatmap_corr_' + name + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'heatmap_corr_' + name + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        figure.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        figure.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()





#corr_coef(merged, df_taxa, 'taxa')
#corr_coef(merged, df_species, 'species')
corr_coef(merged, df_genus, 'genus')

