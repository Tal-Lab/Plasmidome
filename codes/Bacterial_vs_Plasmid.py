"""
Created on 26/04/2023

Author: Lucy Androsiuk
"""
import os, time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import scipy.cluster.hierarchy as sch
from scipy import stats
import fastcluster
from sklearn import preprocessing
from itertools import product
from scipy.stats import fisher_exact
from Bio import Entrez
from Bio import SeqIO
from mlxtend.frequent_patterns import apriori
from mlxtend.frequent_patterns import association_rules
import networkx as nx
from joblib import Parallel, delayed
import multiprocessing as mp
from multiprocessing import Pool
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm


Entrez.email = "androsiu@post.bgu.ac.il"
Entrez.sleep_between_tries = 30
Entrez.api_key = 'ca6eb8a4fb72d1093bbe60949306958cdc08'

# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

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
physical = r"../res/station_physical.xlsx"
lineage = r"../res/rankedlineage.dmp"

#### all taxa table
df_taxa = pd.read_csv(all_taxa,header = 0)
df_taxa['St_Depth'] = df_taxa['Station'].astype(str).str.cat(df_taxa['Depth (m)'].astype(str), sep='_')

#### all genus table
df_genus = pd.read_csv(all_genus,header = 0)
df_genus['St_Depth'] = df_genus['Station'].astype(str).str.cat(df_genus['Depth (m)'].astype(str), sep='_')
df_genus = df_genus.drop(df_genus.columns[-2], axis = 1)
x = df_genus.iloc[:, 17:-1].values  # returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
df_genus_upd = pd.DataFrame(x_scaled, columns = df_genus.iloc[:, 17:-1].columns)
df_genus_upd['St_Depth'] = df_genus['St_Depth']

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

df_lineage = pd.read_table(lineage, sep = "\t|\t", header = 0, engine = 'python')
drop_idx = list(range(1, df_lineage.shape[1], 2))
drop_cols = [j for i, j in enumerate(df_lineage.columns) if i in drop_idx]  # <--
df_lineage = df_lineage.drop(drop_cols, axis = 1)
cols = ['tax_id', 'tax_name', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']
df_lineage.columns = cols

def get_family_for_genus(genus):
    #(genus)
    # Search for the genus using rankedlineage.dmp
    if genus in df_lineage['genus'].unique():
        df_gen = df_lineage.loc[df_lineage['genus'] == genus]
        family = df_gen['family'].iloc[0]
        phylum = df_gen['phylum'].iloc[0]
    else:
        handle = Entrez.esearch(db = 'taxonomy', term = genus)
        record = Entrez.read(handle)
        if record['Count'] != '0':
            id_list = record['IdList']
            handle = Entrez.efetch(db = 'taxonomy', id = id_list[0], retmode = 'xml')
            record = Entrez.read(handle)
            handle.close()
            lineage = record[0]['LineageEx']
            if [item['ScientificName'] for item in lineage if item['Rank'] == 'family']:
                family = [item['ScientificName'] for item in lineage if item['Rank'] == 'family'][0]
            else:
                family = 'Unknown'
            if [item['ScientificName'] for item in lineage if item['Rank'] == 'phylum']:
                phylum = [item['ScientificName'] for item in lineage if item['Rank'] == 'phylum'][0]
            else:
                phylum = 'Unknown'
        else:
            family = 'Unknown'
            phylum = 'Unknown'

    #print(family, phylum)
    # Return the taxonomy information
    return family, phylum

def corr_coef(plasmid_df ,bact_df, name):
    #bullshit
    # create a new dataframe with only the plasmid and bacterium columns
    plasmid_bacteria_df = pd.concat([plasmid_df.iloc[:,:-1], bact_df.iloc[:, 17:-1]], axis=1)
    # calculate the correlation coefficients between all pairs of plasmids and bacteria
    corr_matrix = plasmid_bacteria_df.corr()
    #print(corr_matrix)
    bacteria_to_remove = list(bact_df.columns.values)[17:-1]
    #print(bacteria_to_remove)
    plasmids_to_remove = list(plasmid_df.columns.values)[:-1]
    #print(plasmids_to_remove)
    cols_to_remove = [b for b in bacteria_to_remove if b in corr_matrix.columns]
    rows_to_remove = [p for p in plasmids_to_remove if p in corr_matrix.index]
    # subset the correlation matrix to exclude the rows and columns to remove
    corr_matrix_subset = corr_matrix.drop(rows_to_remove, axis = 0).drop(cols_to_remove, axis = 1)

    #print(corr_matrix_subset)
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

def Physical(vers):
    df = pd.read_excel(physical, index_col = 0, header = 0)
    df = df.dropna()
    df_for_corr = df.reset_index().drop(['Sample', 'Station'], axis=1)
    df_for_corr = df_for_corr.apply(pd.to_numeric)
    df_for_corr = df_for_corr[['Latitude', 'Salinity', 'Chlorophyll', 'Turbidity', 'Temp.', 'Oxygen', 'Depth', 'Nitrate', 'Phosphate', 'Silicate']]
    #print(df_for_corr.columns)
    #mask = np.triu(np.ones_like(df_for_corr.corr(), dtype=bool))
    #f, ax = plt.subplots(figsize=(9, 6))
    #sns.heatmap(df_for_corr.corr(), mask=mask,annot_kws={"size": 10}, fmt='.2f', vmin=-1, vmax=1, annot=True,cmap='coolwarm')
    sns.heatmap(df_for_corr.corr(), annot_kws={"size": 10}, fmt='.2f', vmin=-1, vmax=1, annot=True, cmap='coolwarm')
    #plt.show()
    svg_name = 'heatmap_phys_' + str(vers) + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'heatmap_phys_' + str(vers) + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format='svg', dpi=gcf().dpi, bbox_inches='tight')
        plt.savefig(png_file, format='png', dpi=gcf().dpi, bbox_inches='tight')
    #plt.show()
    df['St_Depth'] = df['Station'].astype(int).astype(str)+ '_'+df['Depth'].astype(str)
    # print(df['Temp.'].sort_values(ascending = False))
    # temperature ranges should be specified not in the code!
    conditions = [
        (df['Temp.'] <= 23),
        (df['Temp.'] > 23) & (df['Temp.'] <= 25),
        (df['Temp.'] > 25) & (df['Temp.'] <= 29),
        (df['Temp.'] > 29)
    ]
    temps = ['21-23', '23-25', '25-29', '29-32']
    df['Temperature'] = np.select(conditions, temps)
    # columns should be specified not in the code!
    df_phys = df[['St_Depth','Latitude','Depth','Temp.','Temperature', 'Salinity', 'Oxygen', 'Chlorophyll', 'Turbidity', 'Nitrate', 'Phosphate', 'Silicate']]
    #print(df_phys)
    return df_phys, df_for_corr

def Station_Order(order, init_order):
    counter_list = list(enumerate(init_order, 0))
    #print("prinying counter list")
    #print(counter_list)
    station_order = order
    #print("Printing station order")
    #print(station_order)
    init_df = pd.DataFrame(counter_list, columns = ['Number', 'ID']).set_index('Number').reindex(station_order)
    #print(init_df)
    final_order = init_df['ID'].to_list()
    #print(final_order)
    return final_order

def Clust_map (vers, df, name, cl, pl):
    #coverage = data_plas(df)
    df = df.set_index('St_Depth')
    # print(coverage.index)
    parameters = Physical(1)[0]
    parameters = parameters.set_index(parameters['St_Depth'])

    figure1 = sns.clustermap(data = df,
                             metric = "euclidean",
                             method = 'ward',
                             cmap = "viridis",
                             xticklabels = False,
                             yticklabels = True,
                             cbar_kws = {'label': 'AP', "orientation": "horizontal"},
                             cbar_pos = (.2, .88, .71, .03)
                             )
    figure1.ax_col_dendrogram.remove()
    figure1.ax_row_dendrogram.remove()
    # station clusters indices correspond to indices of original df

    L_st = figure1.dendrogram_row.linkage
    clusters_st = sch.fcluster(L_st, cl, 'distance')
    #print(clusters_st)
    rows = []
    for i, cluster in enumerate(clusters_st):
        rows.append([df.index[i], cluster])
    cluster_st_df = pd.DataFrame(rows, columns = ["St_Depth", "Sampling points clusters"])
    #print(cluster_st_df)
    cluster_st_df = cluster_st_df.set_index(parameters['St_Depth'])
    stat_cluster = dict(
        zip(cluster_st_df['Sampling points clusters'].unique(),
            sns.color_palette("colorblind", cluster_st_df['Sampling points clusters'].nunique())))
    #up_dict = {5: (1.0, 0.76, 0.04)}
    #stat_cluster.update(up_dict)
    cluster_st = cluster_st_df['Sampling points clusters'].map(stat_cluster)
    #print(cluster_st)
    stations = cluster_st_df.index.values.tolist()

    # plasmid clusters indices correspond to indices of original df

    L_pl = figure1.dendrogram_col.linkage
    clusters_pl = sch.fcluster(L_pl, pl, 'distance')
    cols = []
    for i, cluster in enumerate(clusters_pl):
        cols.append([df.columns[i], cluster])
    cluster_pl_df = pd.DataFrame(cols, columns = ["Plasmids", "Plasmid candidates clusters"])
    cluster_pl_df = cluster_pl_df.set_index('Plasmids')
    plas_cluster = dict(
        zip(cluster_pl_df['Plasmid candidates clusters'].unique(),
            sns.color_palette("Greys", cluster_pl_df['Plasmid candidates clusters'].nunique())))
    cluster_pl = cluster_pl_df['Plasmid candidates clusters'].map(plas_cluster)
    plasmids = cluster_pl_df.index.values.tolist()

    empty = 0 * len(stations)
    for_df = {'St_Depth': stations, 'Empty': empty}
    space_df = pd.DataFrame(for_df)
    space_df = space_df.set_index(parameters['St_Depth'])
    space_df.columns = ['St_Depth', ' ']
    stat_space = dict(zip(space_df[' '].unique(), "white"))
    space_st = space_df[' '].map(stat_space)
    row_colors = pd.concat([cluster_st, space_st], axis = 1)
    sns.set(font_scale = 3.2)
    figure2 = sns.clustermap(data = df,
                             metric = "euclidean",
                             method = 'ward',
                             row_colors = row_colors,
                             col_colors = cluster_pl,
                             figsize = (30, 30),
                             cmap = sns.color_palette("Blues", as_cmap = True),
                             linewidths = 0.0,
                             xticklabels = False,
                             yticklabels = True,
                             rasterized = True,
                             cbar_kws = {"ticks": [0, 1]})
    # figure2.ax_heatmap.tick_params(axis = 'both', which = 'major', pad = 50)
    figure2.fig.subplots_adjust(left = -.01, right = 0.8)
    figure2.cax.set_title("Coverage (%)", pad = 30.0, fontdict = {'fontsize': 36})
    figure2.ax_cbar.set_position((0.92, .36, .03, .3))
    figure2.ax_col_dendrogram.remove()
    figure2.ax_row_dendrogram.remove()

    # GAIW sample labels
    for tick_label in figure2.ax_heatmap.axes.get_yticklabels():
        tick_text = tick_label.get_text()
        if tick_text == "12_47":
            tick_label.set_color('red')
        elif tick_text == "34_50":
            tick_label.set_color('red')
        elif tick_text == "34_100":
            tick_label.set_color('red')

    plt.setp(figure1.ax_heatmap.yaxis.get_majorticklabels())
    plt.setp(figure1.ax_heatmap.xaxis.get_majorticklabels())
    svg_name = name + str(vers) + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = name + str(vers) + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        plt.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()  # Push new figure on stack
    station_order = df.index.values.tolist()
    station_reorder = figure2.dendrogram_row.reordered_ind
    pl_oder = df.columns.values.tolist()
    pl_reorder = figure2.dendrogram_col.reordered_ind
    return station_order, station_reorder, pl_oder, pl_reorder, cluster_st_df, cluster_pl_df, df,cluster_st, cluster_pl

def Clust_bact(vers, df,name, bact,order_init,order_new, order):
    df = df.set_index('St_Depth')
    station_order = Station_Order(order_new, order_init)
    stat_cluster = dict(
        zip(order['Sampling points clusters'].unique(),
            sns.color_palette("colorblind", order['Sampling points clusters'].nunique())))

    # up_dict = {5: (1.0, 0.76, 0.04)}
    # stat_cluster.update(up_dict)
    cluster_st = order['Sampling points clusters'].map(stat_cluster)
    #print(cluster_st)
    df = df.reindex(station_order)
    stations = df.index.values.tolist()
    empty = 0 * len(stations)
    for_df = {'St_Depth': stations, 'Empty': empty}
    space_df = pd.DataFrame(for_df)
    #space_df = space_df.set_index('St_Depth')
    space_df.columns = ['St_Depth', ' ']
    stat_space = dict(zip(space_df[' '].unique(), "white"))
    space_st = space_df[' '].map(stat_space)
    space_df = space_df.set_index('St_Depth')
    row_colors = pd.concat([cluster_st, space_st], axis = 1)
    sns.set(font_scale = 3.2)
    figure = sns.clustermap(data = df,
                             metric = "euclidean",
                             method = 'ward',
                             row_cluster = False,
                             row_colors = row_colors,
                             figsize = (30, 30),
                             cmap = sns.color_palette("Blues", as_cmap = True),
                             linewidths = 0.0,
                             xticklabels = False,
                             yticklabels = True,
                             rasterized = True,
                             cbar_kws = {"ticks": [0, 100]})
    # figure2.ax_heatmap.tick_params(axis = 'both', which = 'major', pad = 50)
    figure.fig.subplots_adjust(left = -.01, right = 0.8)
    figure.cax.set_title("Coverage (%)", pad = 30.0, fontdict = {'fontsize': 36})
    figure.ax_cbar.set_position((0.92, .36, .03, .3))
    figure.ax_col_dendrogram.remove()
    figure.ax_row_dendrogram.remove()

    L_bact = figure.dendrogram_col.linkage
    clusters_bact = sch.fcluster(L_bact, bact, 'distance')
    cols = []
    for i, cluster in enumerate(clusters_bact):
        cols.append([df.columns[i], cluster])
    clusters_bact_df = pd.DataFrame(cols, columns = ["Bacteria", "Bacterial genus clusters"])
    clusters_bact_df = clusters_bact_df.set_index('Bacteria')
    bact_cluster = dict(
        zip(clusters_bact_df['Bacterial genus clusters'].unique(),
            sns.color_palette("Greys", clusters_bact_df['Bacterial genus clusters'].nunique())))
    cluster_bact = clusters_bact_df['Bacterial genus clusters'].map(bact_cluster)
    bacteria = clusters_bact_df.index.values.tolist()

    empty = 0 * len(stations)
    for_df = {'St_Depth': stations, 'Empty': empty}
    space_df = pd.DataFrame(for_df)
    # space_df = space_df.set_index('St_Depth')
    space_df.columns = ['St_Depth', ' ']
    stat_space = dict(zip(space_df[' '].unique(), "white"))
    space_st = space_df[' '].map(stat_space)
    space_df = space_df.set_index('St_Depth')
    row_colors = pd.concat([cluster_st, space_st], axis = 1)
    sns.set(font_scale = 3.2)
    df = df.T
    figure2 = sns.clustermap(data = df,
                             metric = "euclidean",
                             method = 'ward',
                             col_cluster = False,
                             row_colors = cluster_bact,
                             col_colors = row_colors,
                             figsize = (30, 30),
                             cmap = sns.color_palette("Blues", as_cmap = True),
                             linewidths = 0.0,
                             xticklabels = True,
                             yticklabels = False,
                             rasterized = True,
                             cbar_kws = {"ticks": [0, 1]})
    # figure2.ax_heatmap.tick_params(axis = 'both', which = 'major', pad = 50)
    figure2.fig.subplots_adjust(left = -.01, right = 0.8)
    figure2.cax.set_title("Coverage (%)", pad = 30.0, fontdict = {'fontsize': 36})
    figure2.ax_cbar.set_position((0.92, .36, .03, .3))
    figure2.ax_col_dendrogram.remove()
    figure2.ax_row_dendrogram.remove()

    # GAIW sample labels
    for tick_label in figure2.ax_heatmap.axes.get_xticklabels():
        tick_text = tick_label.get_text()
        if tick_text == "12_47":
            tick_label.set_color('red')
        elif tick_text == "34_50":
            tick_label.set_color('red')
        elif tick_text == "34_100":
            tick_label.set_color('red')

    svg_name = name + str(vers) + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = name + str(vers) + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        plt.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()  # Push new figure on stack
    bact_order = df.index.values.tolist()
    bact_reorder = figure2.dendrogram_row.reordered_ind
    return bact_order, bact_reorder, clusters_bact_df, df,cluster_bact

def contingency(plasmid_df ,bact_df, name):
    # Step 1: Create a new DataFrame containing only the plasmid and bacterium columns
    plas_df = plasmid_df.iloc[:,:-1].applymap(lambda x: 1 if x > 0.99 else 0)
    plas_df['St_Depth'] = plasmid_df['St_Depth']
    df = bact_df.iloc[:, :-1].applymap(lambda x: 0 if x < 0.99 else 1)
    df['St_Depth'] = bact_df['St_Depth']
    #bact_df = df.applymap(lambda x: 0 if x < 0.99 else 1)
    merged_df = pd.merge(plas_df, df, on='St_Depth')

    # generate all the possible pairs of plasmids and bacteria
    pairs = list(product(plasmid_df.columns[:-1], df.columns[:-1]))
    # calculate the contingency table and odds ratio for each pair
    results = []
    for plasmid, bacterium in pairs:
        contingency_table = pd.crosstab(merged_df[plasmid], merged_df[bacterium])
        odds_ratio, p_value = fisher_exact(contingency_table)
        results.append((plasmid, bacterium, odds_ratio, p_value))

    # create a new dataframe with the results
    result_df = pd.DataFrame(results, columns = ['plasmid', 'bacterium', 'odds_ratio', 'p_value'])
    #print(result_df)
    # print(df[['NewName', 'station_name']])
    out_file = f'{tables}/odds_ratio.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        result_df.to_csv(out_file, index = False)

def vis_contingency():
    df = pd.read_csv(f'{tables}/odds_ratio.csv', header = 0, index_col = None)
    df = df.sort_values(by = ['p_value'],ascending = False)
    print(df.head(20))
    df = df.sort_values(by = ['p_value'],ascending = True)
    print(df.head(20))

def generate_rules(frequent_itemset):
    return association_rules(frequent_itemset, metric="lift", min_threshold=1)

def association_rules():
    plas_df = merged.iloc[:, :-1].applymap(lambda x: 1 if x > 0.99 else 0)
    plas_df['St_Depth'] = merged['St_Depth']
    df = df_genus_upd.iloc[:, :-1].applymap(lambda x: 0 if x < 0.99 else 1)
    df['St_Depth'] = df_genus_upd['St_Depth']
    # bact_df = df.applymap(lambda x: 0 if x < 0.99 else 1)
    merged_df = pd.merge(plas_df, df, on = 'St_Depth')
    merged_df = merged_df.drop('St_Depth', axis =1)

    # find frequent itemsets using Apriori algorithm
    frequent_itemsets = apriori(merged_df, min_support = 0.2, use_colnames = True)
    #print(frequent_itemsets)
    # generate association rules from frequent itemsets in parallel
    num_cores = 4  # adjust as needed
    rules_list = Parallel(n_jobs = num_cores)(delayed(generate_rules)(itemset) for itemset in frequent_itemsets)

    # combine the list of rules into a single dataframe
    rules = pd.concat(rules_list)

    # sort the rules by decreasing lift
    rules = rules.sort_values('lift', ascending = False)

    # print the top 10 rules by lift
    print(rules.head(10))

def association_rules2():
    odds_ratio_df= pd.read_csv(f'{tables}/odds_ratio.csv', header = 0, index_col = None)

    # Drop rows with NaN values
    odds_ratio_df.dropna(inplace = True)

    # Convert odds ratio values to binary based on a threshold
    odds_ratio_df['presence'] = (odds_ratio_df['odds_ratio'] >= 1.0).astype(int)

    # Keep only relevant columns
    binary_matrix = odds_ratio_df[['plasmid', 'bacterium', 'presence']]
    # create a pivot table
    pivot_table = binary_matrix.pivot_table(index = 'bacterium', columns = 'plasmid', values = 'presence', fill_value = 0)

    # Find frequent itemsets using Apriori algorithm
    frequent_itemsets = apriori(pivot_table, min_support = 0.2, use_colnames = True)

    # Generate association rules from frequent itemsets
    rules = association_rules(frequent_itemsets, metric = "lift", min_threshold = 1)

    # Sort the rules by decreasing lift
    rules = rules.sort_values('lift', ascending = False)

    # Print the top 10 rules by lift
    print(rules.head(10))

def corr_coef2(plasmid_df ,bact_df, order_pl_init, order_pl_new, order_bact_init, order_bact_new, clusters_bact_df, cluster_pl_df, name):
    # create a new dataframe with only the plasmid and bacterium columns
    plasmid_df = plasmid_df.iloc[:, :-1]
    plasmid_order = Station_Order(order_pl_new, order_pl_init)
    plasmid_df = plasmid_df[plasmid_order]
    bact_df = bact_df.iloc[:, :-1]
    bact_order = Station_Order(order_bact_new, order_bact_init)
    bact_df = bact_df[bact_order]

    plasmid_bacteria_df = pd.concat([plasmid_df, bact_df], axis=1)
    # calculate the correlation coefficients between all pairs of plasmids and bacteria
    corr_matrix = plasmid_bacteria_df.corr()
    #print(corr_matrix)
    bacteria_to_remove = list(bact_df.columns.values)
    #print(bacteria_to_remove)
    plasmids_to_remove = list(plasmid_df.columns.values)
    #print(plasmids_to_remove)
    cols_to_remove = [b for b in bacteria_to_remove if b in corr_matrix.columns]
    rows_to_remove = [p for p in plasmids_to_remove if p in corr_matrix.index]
    # subset the correlation matrix to exclude the rows and columns to remove
    corr_matrix_subset = corr_matrix.drop(rows_to_remove, axis = 0).drop(cols_to_remove, axis = 1)
    sns.set(font_scale = 1.5)
    #print(corr_matrix_subset)
    out_file = f'{tables}/correlation_matrix.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        corr_matrix_subset.to_csv(out_file, index = True)
    # get row_colors
    bact_cluster = dict(
        zip(clusters_bact_df['Bacterial genus clusters'].unique(),
            sns.color_palette("Greys", clusters_bact_df['Bacterial genus clusters'].nunique())))
    cluster_bact = clusters_bact_df['Bacterial genus clusters'].map(bact_cluster)
    # get col_colors
    plas_cluster = dict(
        zip(cluster_pl_df['Plasmid candidates clusters'].unique(),
            sns.color_palette("Greys", cluster_pl_df['Plasmid candidates clusters'].nunique())))
    cluster_pl = cluster_pl_df['Plasmid candidates clusters'].map(plas_cluster)

    figure = sns.clustermap(data = corr_matrix_subset,
                             metric = "euclidean",
                             method = 'ward',
                             cmap = "coolwarm",vmin = -1, vmax = 1, center=0,
                             row_cluster = False,
                             col_cluster = False,
                             row_colors = cluster_bact,
                             col_colors = cluster_pl,
                             xticklabels = False,
                             yticklabels = False,
                             cbar_pos=(0, .4, .03, .2),
                             cbar_kws = {'label': 'Pearson correlation', "ticks": [-1,0, 1]})

    figure.ax_col_dendrogram.remove()
    figure.ax_row_dendrogram.remove()
    svg_name = 'heatmap_corr_' + name + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'heatmap_corr_' + name + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        figure.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        figure.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')

    ### calculating distribution of max correlations for plasmids
    df_dist = corr_matrix_subset
    #print(df_dist.max())
    out_file = f'{tables}/max_correlation_plasmids.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        df_dist.max().to_csv(out_file, index = True)
    sns.set_style("white")
    dist = sns.displot(df_dist, x = df_dist.max(), bins = 20)
    svg_name = 'dist_max_plasmid2_' + name + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'dist_max_plasmid2_' + name + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        dist.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        dist.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')

    ### calculating correlations for clusters
    # getting plasmid candidates clusters
    cluster_pl_df = cluster_pl_df.reset_index()
    cluster_pl_df['Plasmid candidates clusters'] = 'C' + cluster_pl_df['Plasmid candidates clusters'].astype(str)
    plasmid_df = plasmid_df.T
    plasmid_df = plasmid_df.reset_index().rename(columns = {'index': 'Plasmids'})
    #print(plasmid_df)
    out_pl = cluster_pl_df.merge(plasmid_df, on = 'Plasmids')
    out_pl.sort_values('Plasmid candidates clusters', inplace = True)

    out_group_pl = out_pl.groupby('Plasmid candidates clusters').mean()  # calculating coverage average for each plasmid candidates cluster at each sampling point

    out_group_pl = out_group_pl.T
    print(out_group_pl)
    # getting bacterial genus clusters
    clusters_bact_df = clusters_bact_df.reset_index()
    clusters_bact_df['Bacterial genus clusters'] = 'Cb' + clusters_bact_df['Bacterial genus clusters'].astype(str)
    #clusters_bact_df['Family'] = clusters_bact_df['Bacteria'].apply(get_family_for_genus)
    bact_df = bact_df.T
    bact_df = bact_df.reset_index().rename(columns = {'index': 'Bacteria'})
    #print(bact_df)
    out_bact = clusters_bact_df.merge(bact_df, on = 'Bacteria')
    out_bact.sort_values('Bacterial genus clusters', inplace = True)
    out_group_bact = out_bact.groupby(
        'Bacterial genus clusters').mean()  # calculating presence average for each bacterial genus cluster at each sampling point

    out_group_bact = out_group_bact.T
    print(out_group_bact)
    #concatenating bacterial clsuters and plasmid clsuters dfs
    plasmid_bacteria_cl_df = pd.concat([out_group_pl, out_group_bact], axis=1)
    print(plasmid_bacteria_cl_df)
    # calculate the correlation coefficients between all pairs of plasmids and bacteria
    corr_cl_matrix = plasmid_bacteria_cl_df.corr()
    bact_cl_to_remove = list(out_group_bact.columns.values)
    pl_cl_to_remove = list(out_group_pl.columns.values)
    cols_to_remove = [b for b in bact_cl_to_remove if b in corr_cl_matrix.columns]
    rows_to_remove = [p for p in pl_cl_to_remove if p in corr_cl_matrix.index]
    # subset the correlation matrix to exclude the rows and columns to remove
    corr_cl_matrix_subset = corr_cl_matrix.drop(rows_to_remove, axis = 0).drop(cols_to_remove, axis = 1)
    sns.set(font_scale = 1.5)
    #print(corr_cl_matrix_subset)
    out_file = f'{tables}/correlation_clusters_matrix2.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        corr_cl_matrix_subset.to_csv(out_file, index = True)

    figure3 = sns.clustermap(data = corr_cl_matrix_subset,
                            metric = "euclidean",
                            method = 'ward',
                            cmap = "coolwarm", vmin = -1, vmax = 1, center = 0,
                            row_cluster = False,
                            col_cluster = False,
                            xticklabels = True,
                            yticklabels = True,
                            annot = corr_cl_matrix_subset,
                            cbar_pos = (0, .3, .03, .3),
                            cbar_kws = {'label': 'Pearson correlation', "ticks": [-1, 0, 1]})

    svg_name = 'clusters_corr2_' + name + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'clusters_corr2_' + name + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        figure3.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        figure3.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()
    '''
    results = clusters_bact_df['Bacteria'].apply(get_family_for_genus)
    clusters_bact_df['Family'], clusters_bact_df['Phylum'] = zip(*results)
    clusters_bact_df = clusters_bact_df.fillna('Unknown')
    out_file = f'{tables}/Bacteria_clusters.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        clusters_bact_df.to_csv(out_file, index = True)
    '''
    out_file = f'{tables}/Bacteria_clusters.csv'
    clusters_bact_df_full = pd.read_csv(out_file, index_col = 0, header = 0)
    clusters_bact_df_full = clusters_bact_df_full.replace({ 'Bacterial genus clusters':{'Ca': 'Cb1',
                                                               'Cb': 'Cb2',
                                                               'Cc': 'Cb3',
                                                               'Cd': 'Cb4',
                                                               'Ce': 'Cb5',
                                                               'Cf': 'Cb6',
                                                               'Cg': 'Cb7'}})
    df_merged = pd.merge(corr_matrix_subset, clusters_bact_df_full, left_index = True, right_on = 'Bacteria')

    df_merged = df_merged.set_index(['Bacterial genus clusters', 'Phylum', 'Family','Bacteria'])
    order_cols = df_merged.columns.tolist()
    # convert the column to a categorical data type with the specified order
    cluster_pl_df['Plasmids'] = pd.Categorical(cluster_pl_df['Plasmids'], categories = order_cols, ordered = True)
    # sort the DataFrame by the categorical column
    cluster_pl_df = cluster_pl_df.sort_values('Plasmids')
    # make tuples from CLusters and plasmids for multiindex
    cols_corr_matrix = list(cluster_pl_df[['Plasmid candidates clusters', 'Plasmids' ]].itertuples(index = False, name = None))
    cols_corr_matrix = pd.MultiIndex.from_tuples(cols_corr_matrix)
    df_merged.columns = cols_corr_matrix
    #print(df_merged)
    out_file = f'{tables}/correlation_matrix_clusters2.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        df_merged.to_csv(out_file, index = True)
    return corr_matrix_subset,cluster_bact,cluster_pl

def grid_image():
    pl_order = Station_Order(order_pl_new, order_pl_init)
    bact_order = Station_Order(order_bact_new, order_bact_init)
    station_order = Station_Order(order_st_new, order_st_init)
    data1 = corr_bact
    data1 = data1.reindex(bact_order)
    data1 = data1.reindex(columns = station_order)
    data2 = corr_pl
    data2 = data2.reindex(station_order)
    data2 = data2.reindex(columns = pl_order)
    data3, cl_bact, cl_pl = corr_coef2(merged,df_genus_upd,order_pl_init,order_pl_new, order_bact_init, order_bact_new, order_bacteria, order_plasmids, 'genus')
    # Create a 2x2 grid of subplots using gridspec
    fig = plt.figure(figsize = (16, 12))
    gs = gridspec.GridSpec(2, 4, figure = fig, width_ratios=[5.5,0.75,5.5,0.75])

    sns.set(font_scale = 1)
    sns.set_style("white")
    ax1 = fig.add_subplot(gs[0, 0:1])
    distplot = sns.histplot(data = data3, x = data3.max(), bins = 20, ax = ax1)
    sns.set_style("white")

    # Create the second clustermap and place it in the upper left subplot
    cl_st_col=cluster_st.reindex(station_order)
    cl_pl_col=cluster_pl.reindex(pl_order)
    cl_bact_col=cluster_bact.reindex(bact_order)
    print(cl_st_col)
    print(cl_pl_col)
    print(cl_bact_col)

    ax2 = fig.add_subplot(gs[1, 0:2])
    heatmap1 = sns.heatmap(data1,
                           cmap = sns.color_palette("Blues", as_cmap = True),
                           linewidths = 0.0,
                           xticklabels = False,
                           yticklabels = False,
                           cbar_kws = {'label': 'Coverage',"ticks": [0, 1]},
                           ax = ax2)
    heatmap1.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row colors
    for i, color in enumerate(cl_st_col):
        heatmap1.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                   transform = heatmap1.get_xaxis_transform(), clip_on = False))

    #heatmap1.set_yticklabels(cluster_bact_df['Bacterial genus clusters'], rotation = 0)  # optionally use the groups as the tick labels
    for i, color in enumerate(cl_bact_col):
        heatmap1.add_patch(plt.Rectangle(xy = (-0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap1.get_yaxis_transform(), clip_on = False))
    heatmap1.set_xlabel('Sampling points', labelpad = 35)  # optionally use the groups as the tick labels
    heatmap1.set_ylabel('Bacterial genus', labelpad = 35)  # optionally use the groups as the tick labels
    # Create the second clustermap and place it in the lower left subplot
    ax3 = fig.add_subplot(gs[1, 2:])
    heatmap2 = sns.heatmap(data3,
                           cmap = "coolwarm",
                           vmin = -1,
                           vmax = 1,
                           center = 0,
                           xticklabels = False,
                           yticklabels = False,
                           cbar_kws = {'label': 'Pearson correlation', "ticks": [-1, 0, 1]},
                           ax = ax3)
    heatmap2.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row colors
    heatmap2.set_xlabel('Plasmid candidates', labelpad = 35)
    heatmap2.set_ylabel('Bacterial genus', labelpad = 35) # optionally use the groups as the tick labels
    for i, color in enumerate(cl_bact_col):
        heatmap2.add_patch(plt.Rectangle(xy = ( -0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap2.get_yaxis_transform(), clip_on = False))
    #heatmap2.set_xticklabels(cluster_pl_df['Plasmid candidates clusters'], rotation = 0)  # optionally use the groups as the tick labels
    for i, color in enumerate(cl_pl_col):
        heatmap2.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                         transform = heatmap2.get_xaxis_transform(), clip_on = False))
    # Create the third clustermap and place it in the upper right subplot
    ax4 =  fig.add_subplot(gs[0, 2:])
    heatmap3 = sns.heatmap(data2,
                           cmap = sns.color_palette("Blues", as_cmap = True),
                           linewidths = 0.0,
                           xticklabels = False,
                           yticklabels = False,
                           rasterized = True,
                           cbar_kws = {'label': 'Coverage',"ticks": [0, 1]},
                           ax = ax4)
    heatmap3.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row colors
    heatmap3.set_xlabel('Plasmid candidates', labelpad = 35)# optionally use the groups as the tick labels
    heatmap3.set_ylabel('Sampling points', labelpad = 35)
    for i, color in enumerate(cl_pl_col):
        heatmap3.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                         transform = heatmap3.get_xaxis_transform(), clip_on = False))
    #heatmap3.set_yticklabels(stations_cl_df['Sampling points clusters'], rotation = 0)  # optionally use the groups as the tick labels
    for i, color in enumerate(cl_st_col):
        heatmap3.add_patch(plt.Rectangle(xy = (-0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap3.get_yaxis_transform(), clip_on = False))


    svg_name = 'clusters_grid4.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'clusters_grid4.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        fig.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        fig.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    sns.set_style("white")
    plt.tight_layout()
    plt.show()

def Pearson_corr(df1, df2, df_to_update):
    ''' The function uses stats.pearsonr to calculate Pearson correlation between two vectors'''
    for index_o, row_o in df1.iterrows():
        # iterating clusters and getting coverage of the index_o cluster at each row_o
        for index_p, row_p in df2.iterrows():
            # iterating physical parameters and getting parameter of the index_p environmental condition at each row_p sampling point
            print("Pearson correlation, p-value for Plasmid %s and %s" % (index_o, index_p))
            correlation, p_value= stats.pearsonr(row_o, row_p) # calculating Pearson correlation and p-value for cluster:env.condition in each sampling point
            print(round(correlation,3),round(p_value,3))
            df_to_update[index_p][index_o] = (round(correlation,3),round(p_value,3)) # updating df_pearson dataframe at with calculated correlation values
    return df_to_update

def PCA_calc():
    # Load the datasets
    table1 = Physical(1)[0] # physical properties at the stations
    table1 = table1.reset_index().drop(['Sample', 'Temperature'], axis = 1)
    table1 = table1[table1.columns[:-1]]
    table1 = table1.set_index('St_Depth')


    table2 = corr_coef2(merged,df_genus_upd,order_pl_init,order_pl_new, order_bact_init, order_bact_new,order_bacteria,order_plasmids, 'genus')[0]  # bacterial_plasmid_corr
    X = pd.merge(table1, table2, left_index = True, right_index = True)
    print(table2)

    table3 = merged  # plasmid presence at the stations
    table3 = table3.set_index('St_Depth')

    cluster_pl_df = order_plasmids.reset_index()
    cluster_pl_df['Plasmid candidates clusters'] = 'C' + cluster_pl_df['Plasmid candidates clusters'].astype(str)




    df_pearson_pl = pd.DataFrame(columns = table1.T.index.to_list(),
                                 index = table3.T.index.to_list())  # empty dataframe for pearson clusters:env.conditions
    ### getting Pearson correlation for each cluster-env.condition
    Pearson_corr(table3.T, table1.T, df_pearson_pl)
    # df.assign(**df[['col2', 'col3']].apply(lambda x: x.str[0]))
    df_pearson_2 = df_pearson_pl.assign(**df_pearson_pl[df_pearson_pl.columns.to_list()].apply(lambda x: x.str[0]))
    print(df_pearson_2)

    table2 = table2.T
    data = table2.join(df_pearson_2)
    data_ix = data.reset_index().rename(columns = {'index': 'Plasmids'})
    data_full = data_ix.merge(cluster_pl_df, on = 'Plasmids')
    data_new = data_full.set_index(['Plasmid candidates clusters', 'Plasmids'])
    print(data_full.columns)
    features = data.columns.to_list()
    x = data.loc[:, features].values
    x = StandardScaler().fit_transform(x)  # normalizing the features

    print(x)
    print(x.shape)
    print(np.mean(x),np.std(x))

    normalised_plasmid = pd.DataFrame(x, columns = features, index = data.index)
    normalised_plasmid.to_csv(f'{tables}/PCA_pre.txt')
    print(normalised_plasmid)

    pca_plasmids = PCA(n_components = 2)
    principalComponents_plasmids = pca_plasmids.fit_transform(x)

    principal_plasmids_Df = pd.DataFrame(data = principalComponents_plasmids, columns = ['principal component 1', 'principal component 2'], index = normalised_plasmid.index)
    principal_plasmids_Df.to_csv(f'{tables}/PCA.txt')
    print(principal_plasmids_Df.tail())

    print('Explained variation per principal component: {}'.format(pca_plasmids.explained_variance_ratio_))

    # Get the loadings (correlations) of original variables with Principal Component 1
    loadings = pd.DataFrame(pca_plasmids.components_.T, columns = ['PC1', 'PC2'], index = data.columns)
    loadings_pc1 = loadings['PC1'].abs().sort_values(ascending = False)
    loadings_pc2 = loadings['PC2'].abs().sort_values(ascending = False)

    # Identify the variables with the highest correlation to Principal Component 1
    variables_pc1 = loadings_pc1.head(5)  # Adjust the number (e.g., top 5) based on your preference

    print("Variables with to Principal Component 1:")
    print(loadings_pc1)

    print("Variables with the highest correlation to Principal Component 2:")
    print(loadings_pc2)

    plt.figure()
    plt.figure(figsize = (10, 10))
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 14)
    plt.xlabel('Principal Component - 1', fontsize = 20)
    plt.ylabel('Principal Component - 2', fontsize = 20)
    plt.title("Principal Component Analysis of Plasmidome", fontsize = 20)

    targets = cluster_pl_df['Plasmid candidates clusters'].unique().tolist()
    targets.sort()
    colors = sns.color_palette("colorblind", cluster_pl_df['Plasmid candidates clusters'].nunique())

    for target, color in zip(targets, colors):
        indicesToKeep = data_full.loc[data_full['Plasmid candidates clusters'] == target, 'Plasmids']
        print(indicesToKeep)
        plt.scatter(principal_plasmids_Df.loc[indicesToKeep, 'principal component 1'], principal_plasmids_Df.loc[indicesToKeep, 'principal component 2'], c = color, s = 50)

    plt.legend(targets, prop = {'size': 15})
    plt.show()


    X = sm.add_constant(X)  # Add a constant column for the intercept
    print(table3)
    model = sm.OLS(table3, X).fit()
    print(model.summary())

    '''
    # correlation dfs
    
    # Perform PCA
    pca = PCA()
    pca.fit(data)
    pca_result = pca.fit_transform(data)

    # Calculate the percentage of variation explained by each principal component
    explained_variance_ratio = pca.explained_variance_ratio_

    # Determine the contribution percentages of physical properties and bacterial population
    n_components = pca.n_components_
    contribution_percentage_physical = np.sum(explained_variance_ratio[:n_components // 3]) * 100
    contribution_percentage_bacteria = np.sum(explained_variance_ratio[n_components // 3:2 * n_components // 3]) * 100

    plt.rcParams.update({'font.size': 14})
    # Plotting the explained variance ratios
    plt.figure(figsize = (8, 6))
    component_labels = [f"Component {i + 1}" for i in range(n_components)]
    plt.bar(component_labels, explained_variance_ratio)
    plt.xlabel('Principal Components')
    plt.ylabel('Explained Variance Ratio')
    plt.title('Explained Variance Ratio by Principal Components')
    plt.show()

    # Scree plot
    plt.figure(figsize = (8, 6))
    cumulative_variance_ratio = np.cumsum(explained_variance_ratio)
    plt.plot(range(1, n_components + 1), cumulative_variance_ratio, marker = 'o', linestyle = '--')
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Explained Variance Ratio')
    plt.title('Scree Plot')
    plt.show()

    # Create a scatter plot of the PCA results
    plt.figure(figsize = (8, 6))
    plt.scatter(pca_result[:, 0], pca_result[:, 1])
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA Scatter Plot')
    plt.show()

    # Print the contribution percentages
    print(f"Physical Properties Contribution: {contribution_percentage_physical}%")
    print(f"Bacterial Population Contribution: {contribution_percentage_bacteria}%")
    '''

def grid_image2():
    pl_order = Station_Order(order_pl_new, order_pl_init)
    bact_order = Station_Order(order_bact_new, order_bact_init)
    station_order = Station_Order(order_st_new, order_st_init)
    data1 = corr_bact
    data1 = data1.reindex(bact_order)
    data1 = data1.reindex(columns = station_order)
    data2 = corr_pl
    data2 = data2.reindex(station_order)
    data2 = data2.reindex(columns = pl_order)
    data3, cl_bact, cl_pl = corr_coef2(merged,df_genus_upd,order_pl_init,order_pl_new, order_bact_init, order_bact_new, order_bacteria, order_plasmids, 'genus')
    # Create a 2x2 grid of subplots using gridspec
    fig = plt.figure(figsize = (16, 12))
    gs = gridspec.GridSpec(2, 4, figure = fig, width_ratios=[5.5,0.75,5.5,0.75])

    sns.set_theme(font_scale = 1, style= "white", font = 'Helvetica')
    #sns.set_style()
    ax1 = fig.add_subplot(gs[1, 0:1])
    distplot = sns.histplot(data = data3, x = data3.max(), bins = 20, ax = ax1)
    distplot.set_xlabel('Max. Pearson correlation')

    # Create the second clustermap and place it in the upper left subplot
    cl_st_col=cluster_st.reindex(station_order)
    cl_pl_col=cluster_pl.reindex(pl_order)
    cl_bact_col=cluster_bact.reindex(bact_order)
    print(cl_st_col)
    print(cl_pl_col)
    print(cl_bact_col)

    ax2 = fig.add_subplot(gs[0, 0:2])
    heatmap1 = sns.heatmap(data1,
                           cmap = sns.color_palette("Blues", as_cmap = True),
                           linewidths = 0.0,
                           xticklabels = False,
                           yticklabels = False,
                           cbar_kws = {'label': 'Coverage',"ticks": [0, 1]},
                           ax = ax2)
    heatmap1.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row colors
    for i, color in enumerate(cl_st_col):
        heatmap1.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                   transform = heatmap1.get_xaxis_transform(), clip_on = False))

    #heatmap1.set_yticklabels(cluster_bact_df['Bacterial genus clusters'], rotation = 0)  # optionally use the groups as the tick labels
    for i, color in enumerate(cl_bact_col):
        heatmap1.add_patch(plt.Rectangle(xy = (-0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap1.get_yaxis_transform(), clip_on = False))
    heatmap1.set_xlabel('Sampling points', labelpad = 35)  # optionally use the groups as the tick labels
    heatmap1.set_ylabel('Bacterial genus', labelpad = 35)  # optionally use the groups as the tick labels
    # Create the second clustermap and place it in the lower left subplot
    ax3 = fig.add_subplot(gs[0, 2:])
    heatmap2 = sns.heatmap(data3,
                           cmap = "coolwarm",
                           vmin = -1,
                           vmax = 1,
                           center = 0,
                           xticklabels = False,
                           yticklabels = False,
                           cbar_kws = {'label': 'Pearson correlation', "ticks": [-1, 0, 1]},
                           ax = ax3)
    heatmap2.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row colors
    heatmap2.set_xlabel('Plasmid candidates', labelpad = 35)
    heatmap2.set_ylabel('Bacterial genus', labelpad = 35) # optionally use the groups as the tick labels
    for i, color in enumerate(cl_bact_col):
        heatmap2.add_patch(plt.Rectangle(xy = ( -0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap2.get_yaxis_transform(), clip_on = False))
    #heatmap2.set_xticklabels(cluster_pl_df['Plasmid candidates clusters'], rotation = 0)  # optionally use the groups as the tick labels
    for i, color in enumerate(cl_pl_col):
        heatmap2.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                         transform = heatmap2.get_xaxis_transform(), clip_on = False))
    # Create the third clustermap and place it in the upper right subplot
    ax4 =  fig.add_subplot(gs[1, 2:])
    heatmap3 = sns.heatmap(data2,
                           cmap = sns.color_palette("Blues", as_cmap = True),
                           linewidths = 0.0,
                           xticklabels = False,
                           yticklabels = False,
                           rasterized = True,
                           cbar_kws = {'label': 'Coverage',"ticks": [0, 1]},
                           ax = ax4)
    heatmap3.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row colors
    heatmap3.set_xlabel('Plasmid candidates', labelpad = 35)# optionally use the groups as the tick labels
    heatmap3.set_ylabel('Sampling points', labelpad = 35)
    for i, color in enumerate(cl_pl_col):
        heatmap3.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                         transform = heatmap3.get_xaxis_transform(), clip_on = False))
    #heatmap3.set_yticklabels(stations_cl_df['Sampling points clusters'], rotation = 0)  # optionally use the groups as the tick labels
    for i, color in enumerate(cl_st_col):
        heatmap3.add_patch(plt.Rectangle(xy = (-0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap3.get_yaxis_transform(), clip_on = False))


    svg_name = 'clusters_grid_reord.eps'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'clusters_grid_reord.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        fig.savefig(svg_file, format = 'eps', dpi = gcf().dpi, bbox_inches = 'tight')
        fig.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.tight_layout()
    plt.show()
#corr_coef(merged, df_taxa, 'taxa')
#corr_coef(merged, df_species, 'species')
#contingency(merged, df_genus_upd, 'genus')
#Clust_map(1,merged,'PlPutUnc_HMannot_', 1150, 1200)
#Clust_map(1,merged,'PlPutUnc_HMannot_', 11.5, 12)
order_st_init, order_st_new, order_pl_init, order_pl_new, order_stations ,order_plasmids,corr_pl,cluster_st,cluster_pl= Clust_map(1,merged,'Plasm_for_bact_', 11.5, 12)
order_bact_init, order_bact_new, order_bacteria, corr_bact, cluster_bact = Clust_bact(1,df_genus_upd,'Bact_HMannot_', 8, order_st_init,order_st_new,order_stations)
#corr_coef2(merged,df_genus_upd,order_pl_init,order_pl_new, order_bact_init, order_bact_new,order_bacteria,order_plasmids, 'genus')
#vis_contingency()
#association_rules2()
#grid_image2()
#PCA_calc()