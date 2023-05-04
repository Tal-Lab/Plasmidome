"""
Created on 26/04/2023

Author: Lucy Androsiuk
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import scipy.cluster.hierarchy as sch
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


Entrez.email = "androsiu@post.bgu.ac.il"

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

def get_family_for_genus(genus):
    Entrez.sleep_between_tries = 20
    handle = Entrez.esearch(db="taxonomy", term=genus)
    record = Entrez.read(handle)
    if len(record["IdList"]) == 0:
        return None
    id = record["IdList"][0]
    handle = Entrez.efetch(db="taxonomy", id=id, retmode="xml")
    record = Entrez.read(handle)
    for r in record[0]["LineageEx"]:
        if r["Rank"] == "family":
            return r["ScientificName"]
    return None

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

def Physical(vers):
    df = pd.read_excel(physical, index_col = 0, header = 0)
    df = df.dropna()
    df_for_corr = df.reset_index().drop(['Sample', 'Station'], axis=1)
    df_for_corr = df_for_corr.apply(pd.to_numeric)
    df_for_corr = df_for_corr[['Latitude', 'Salinity', 'Chlorophyll', 'Turbidity', 'Temp.', 'Oxygen', 'Depth', 'Nitrate', 'Phosphate', 'Silicate']]
    print(df_for_corr.columns)
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
    plt.show()
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
    print(clusters_st)
    rows = []
    for i, cluster in enumerate(clusters_st):
        rows.append([df.index[i], cluster])
    cluster_st_df = pd.DataFrame(rows, columns = ["St_Depth", "Sampling points clusters"])
    print(cluster_st_df)
    cluster_st_df = cluster_st_df.set_index(parameters['St_Depth'])
    stat_cluster = dict(
        zip(cluster_st_df['Sampling points clusters'].unique(),
            sns.color_palette("colorblind", cluster_st_df['Sampling points clusters'].nunique())))
    #up_dict = {5: (1.0, 0.76, 0.04)}
    #stat_cluster.update(up_dict)
    cluster_st = cluster_st_df['Sampling points clusters'].map(stat_cluster)
    print(cluster_st)
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
    plt.show()  # Push new figure on stack
    station_order = df.index.values.tolist()
    station_reorder = figure2.dendrogram_row.reordered_ind
    pl_oder = df.columns.values.tolist()
    pl_reorder = figure2.dendrogram_col.reordered_ind
    return station_order, station_reorder, pl_oder, pl_reorder, cluster_st_df, cluster_pl_df

def Clust_bact(vers, df,name, bact,order_init,order_new, order):
    df = df.set_index('St_Depth')
    station_order = Station_Order(order_new, order_init)
    stat_cluster = dict(
        zip(order['Sampling points clusters'].unique(),
            sns.color_palette("colorblind", order['Sampling points clusters'].nunique())))

    # up_dict = {5: (1.0, 0.76, 0.04)}
    # stat_cluster.update(up_dict)
    cluster_st = order['Sampling points clusters'].map(stat_cluster)
    print(cluster_st)
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
    plt.show()  # Push new figure on stack
    bact_order = df.index.values.tolist()
    bact_reorder = figure2.dendrogram_row.reordered_ind
    return bact_order, bact_reorder, clusters_bact_df

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
    print(result_df)
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
    print(frequent_itemsets)
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
    print(corr_matrix)
    bacteria_to_remove = list(bact_df.columns.values)
    print(bacteria_to_remove)
    plasmids_to_remove = list(plasmid_df.columns.values)
    print(plasmids_to_remove)
    cols_to_remove = [b for b in bacteria_to_remove if b in corr_matrix.columns]
    rows_to_remove = [p for p in plasmids_to_remove if p in corr_matrix.index]
    # subset the correlation matrix to exclude the rows and columns to remove
    corr_matrix_subset = corr_matrix.drop(rows_to_remove, axis = 0).drop(cols_to_remove, axis = 1)
    sns.set(font_scale = 1.5)
    print(corr_matrix_subset)

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
    print(df_dist.max())
    out_file = f'{tables}/max_correlation_plasmids.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        df_dist.max().to_csv(out_file, index = False)
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
    print(plasmid_df)
    out_pl = cluster_pl_df.merge(plasmid_df, on = 'Plasmids')
    out_pl.sort_values('Plasmid candidates clusters', inplace = True)
    out_group_pl = out_pl.groupby('Plasmid candidates clusters').mean()  # calculating coverage average for each plasmid candidates cluster at each sampling point

    out_group_pl = out_group_pl.T
    # getting bacterial genus clusters
    clusters_bact_df = clusters_bact_df.reset_index()
    clusters_bact_df['Bacterial genus clusters'] = clusters_bact_df['Bacterial genus clusters'].astype(int).apply(
        lambda x: chr(ord('`') + x))
    clusters_bact_df['Bacterial genus clusters'] = 'C' + clusters_bact_df['Bacterial genus clusters']
    #clusters_bact_df['Family'] = clusters_bact_df['Bacteria'].apply(get_family_for_genus)
    bact_df = bact_df.T
    bact_df = bact_df.reset_index().rename(columns = {'index': 'Bacteria'})
    print(bact_df)
    out_bact = clusters_bact_df.merge(bact_df, on = 'Bacteria')
    out_bact.sort_values('Bacterial genus clusters', inplace = True)
    out_group_bact = out_bact.groupby(
        'Bacterial genus clusters').mean()  # calculating presence average for each bacterial genus cluster at each sampling point

    out_group_bact = out_group_bact.T

    #concatenating bacterial clsuters and plasmid clsuters dfs
    plasmid_bacteria_cl_df = pd.concat([out_group_pl, out_group_bact], axis=1)
    # calculate the correlation coefficients between all pairs of plasmids and bacteria
    corr_cl_matrix = plasmid_bacteria_cl_df.corr()
    print(corr_cl_matrix)
    bact_cl_to_remove = list(out_group_bact.columns.values)
    print(bact_cl_to_remove)
    pl_cl_to_remove = list(out_group_pl.columns.values)
    print(pl_cl_to_remove)
    cols_to_remove = [b for b in bact_cl_to_remove if b in corr_cl_matrix.columns]
    rows_to_remove = [p for p in pl_cl_to_remove if p in corr_cl_matrix.index]
    # subset the correlation matrix to exclude the rows and columns to remove
    corr_cl_matrix_subset = corr_cl_matrix.drop(rows_to_remove, axis = 0).drop(cols_to_remove, axis = 1)
    sns.set(font_scale = 1.5)
    print(corr_cl_matrix_subset)

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

    svg_name = 'clusters_corr_' + name + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'clusters_corr_' + name + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        figure3.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        figure3.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()






#corr_coef(merged, df_taxa, 'taxa')
#corr_coef(merged, df_species, 'species')
#contingency(merged, df_genus_upd, 'genus')
#Clust_map(1,merged,'PlPutUnc_HMannot_', 1150, 1200)
#Clust_map(1,merged,'PlPutUnc_HMannot_', 11.5, 12)
#order_st_init, order_st_new, order_pl_init, order_pl_new, order_stations ,order_plasmids= Clust_map(1,merged,'Plasm_for_bact_', 11.5, 12)
#order_bact_init, order_bact_new, order_bacteria = Clust_bact(1,df_genus_upd,'Bact_HMannot_', 8, order_st_init,order_st_new,order_stations)
#corr_coef2(merged,df_genus_upd,order_pl_init,order_pl_new, order_bact_init, order_bact_new,order_bacteria,order_plasmids, 'genus')
#vis_contingency()
#association_rules2()