# -*- coding: utf-8 -*-
"""
Created on 10/08/2021 8:21

Author: Lucy Androsiuk
"""
# add description
''' This code generates dataframe with plasmid candidates coverage at each sampling station and clusters it. 
    It also uses physical parameters such as station depth and temperature '''


import numpy as np
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
import pylab
from pylab import *
import scipy.cluster.hierarchy as sch
from collections import defaultdict
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import os
from general_analysis import plasmids_by_reads, tables, Station, GetLibSize, CoverageDF, version
from plasmid_detect import Plasmid_class
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

#from CovBAM2 import out_file as datafile

# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
datafile = r"../res/all_cov.csv"
physical = r"../res/station_physical.xlsx"

def data_file():
    ' this function reads coverage file as a dataframe '
    data = pd.read_csv(datafile, sep = ',', index_col = 0, header = 0)
    data.index = data.index.str.replace(r'_l\w+', '',  regex=True)
    data = data.drop('94_LNODE_1', axis=0)
    df=Station()
    data = data.rename(columns = df.set_index('station_name')['St_Depth'].to_dict())
    data = data.rename_axis("Plasmid candidates")
    data = data.rename_axis("Sampling points", axis="columns")
    # transpose data to fit it in clustermap with stations on the y axis and candidates on the x axis
    data_transp = data.transpose()
    return data_transp

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
#Physical(1)

### getting statistics on plasmids by reads
# does it belong here?
def Cov_plasmids():
    df = DF_plasmids_byReads()[0]
    df_count_stat = df.groupby(['NewName']).size().reset_index(name = 'counts')
    df_one_stat = df_count_stat.loc[df_count_stat['counts']==1]
    df_more_stat = df_count_stat.loc[df_count_stat['counts']>1]
    number_more_stat = df_more_stat['NewName'].nunique()
    number_one_stat = df_one_stat['NewName'].nunique()
    #print(number_more_stat)
    #print(df_one_stat['NewName'].unique())
    number_all = df['NewName'].nunique()
    perc_one = (number_one_stat/number_all)*100
    perc_more = (number_more_stat/number_all)*100
    #print(perc_one)
    #print(perc_more)
#Cov_plasmids()

def data_plas(df_pl):
    df_cov = data_file()
    #print(df_pl)
    #df_pl = Plasmid_class()[0]
    final_table_columns = df_pl['Plasmid'].unique()
    df_cov = df_cov.drop(columns=[col for col in df_cov if col not in final_table_columns])
    #print(df_cov.shape)
    return df_cov

def Clust_map2(vers, df, name, cl, pl):
    coverage = data_plas(df)
    #print(coverage.index)
    parameters = Physical(1)[0]
    parameters = parameters.set_index(parameters['St_Depth'])
    # maybe put parameters in loop?
    parameters["Temperature"] = parameters['Temp.'].astype(float).round().astype(int)
    stat_temp = dict(zip(parameters["Temperature"].sort_values().unique(), sns.color_palette("vlag", parameters['Temperature'].nunique())))
    temperature = parameters["Temperature"].map(stat_temp)
    stat_depth = dict(zip(parameters['Depth'].sort_values().unique(), sns.color_palette("PuBu", parameters['Depth'].nunique())))
    depth = parameters['Depth'].map(stat_depth)
    parameters["Salinity"] = parameters['Salinity'].astype(float).round().astype(int)
    stat_salt = dict(zip(parameters['Salinity'].sort_values().unique(), sns.color_palette("YlGn", parameters['Salinity'].nunique())))
    salt = parameters['Salinity'].map(stat_salt)
    parameters["Oxygen"] = parameters['Oxygen'].astype(float).round().astype(int)
    stat_O2 = dict(zip(parameters['Oxygen'].sort_values().unique(), sns.color_palette("ocean_r", parameters['Oxygen'].nunique())))
    O2 = parameters['Oxygen'].map(stat_O2)
    parameters["Nitrate"] = parameters['Nitrate'].astype(float).round().astype(int)
    stat_NO3 = dict(zip(parameters['Nitrate'].sort_values().unique(), sns.color_palette("YlGn", parameters['Nitrate'].nunique())))
    NO3 = parameters['Nitrate'].map(stat_NO3)
    parameters["Phosphate"] = parameters['Phosphate'].astype(float).round(1)
    stat_PO4 = dict(zip(parameters['Phosphate'].sort_values().unique(), sns.color_palette("Purples", parameters['Phosphate'].nunique())))
    PO4 = parameters['Phosphate'].map(stat_PO4)
    parameters["Silicate"] = parameters['Silicate'].astype(float).round().astype(int)
    stat_Si = dict(zip(parameters['Silicate'].sort_values().unique(),sns.color_palette("Oranges", parameters['Silicate'].nunique())))
    Si = parameters['Silicate'].map(stat_Si)
    figure1 = sns.clustermap(data = coverage,
                             metric = "euclidean",
                             method = 'ward',
                             row_colors = [temperature, depth, salt, O2, NO3, PO4, Si],
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
    rows = []
    for i, cluster in enumerate(clusters_st):
        rows.append([coverage.index[i], cluster])
    cluster_st_df = pd.DataFrame(rows, columns=["St_Depth", "Sampling points clusters"])
    cluster_st_df = cluster_st_df.set_index(parameters['St_Depth'])
    stat_cluster = dict(
        zip(cluster_st_df['Sampling points clusters'].unique(), sns.color_palette("colorblind", cluster_st_df['Sampling points clusters'].nunique())))
    cluster_st = cluster_st_df['Sampling points clusters'].map(stat_cluster)
    stations = cluster_st_df.index.values.tolist()
    
    # plasmid clusters indices correspond to indices of original df

    L_pl = figure1.dendrogram_col.linkage
    clusters_pl = sch.fcluster(L_pl, pl, 'distance')
    cols = []
    for i, cluster in enumerate(clusters_pl):
        cols.append([coverage.columns[i], cluster])
    cluster_pl_df = pd.DataFrame(cols, columns = ["Plasmids", "Plasmid candidates clusters"])
    cluster_pl_df = cluster_pl_df.set_index('Plasmids')
    plas_cluster = dict(
        zip(cluster_pl_df['Plasmid candidates clusters'].unique(), sns.color_palette("Greys", cluster_pl_df['Plasmid candidates clusters'].nunique())))
    cluster_pl = cluster_pl_df['Plasmid candidates clusters'].map(plas_cluster)
    plasmids = cluster_pl_df.index.values.tolist()

    empty = 0*len(stations)
    for_df = {'St_Depth':stations, 'Empty': empty}
    space_df = pd.DataFrame(for_df)
    space_df=space_df.set_index(parameters['St_Depth'])
    space_df.columns = ['St_Depth', ' ']
    stat_space = dict(zip(space_df[' '].unique(), "white"))
    space_st = space_df[' '].map(stat_space)
    row_colors = pd.concat([cluster_st, salt, temperature, depth, space_st], axis = 1)
    sns.set(font_scale = 3.2)
    figure2 = sns.clustermap(data = coverage,
                             metric = "euclidean",
                             method = 'ward',
                             row_colors = row_colors,
                             col_colors=cluster_pl,
                             figsize = (30, 30),
                             cmap = sns.color_palette("Blues", as_cmap = True),
                             linewidths = 0.0,
                             xticklabels = False,
                             yticklabels = True,
                             rasterized = True,
                             cbar_kws = {"ticks": [0, 100]})
    #figure2.ax_heatmap.tick_params(axis = 'both', which = 'major', pad = 50)
    figure2.fig.subplots_adjust(left = -.01, right = 0.8)
    figure2.cax.set_title("Coverage(%)", pad=30.0, fontdict={'fontsize':36})
    figure2.ax_cbar.set_position((0.92, .26, .03, .5))
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

    # temperature legend
    tem_legend = []
    for label in parameters["Temperature"].sort_values().unique():
        temp_patch = figure2.ax_row_dendrogram.bar(0, 0, color = stat_temp[label], label = label, linewidth=0)
        tem_legend.append(temp_patch)
    l1 = plt.legend(tem_legend, parameters["Temperature"].sort_values().unique(), loc="upper left", title = 'Temperature', title_fontsize= 'large', ncol=2, bbox_to_anchor=(-0.1, 0.8), bbox_transform=gcf().transFigure, fontsize='large', frameon=False)
    plt.gca().add_artist(l1)
    # depth legend
    # loop?
    dep_legend = []
    for label in parameters['Depth'].sort_values().unique():
        dep_patch = figure2.ax_row_dendrogram.bar(0, 0, color = stat_depth[label], label = label, linewidth=0)
        dep_legend.append(dep_patch)
    l2 = plt.legend(dep_legend, parameters['Depth'].sort_values().unique(), loc= "center left",  title = 'Depth',title_fontsize= 'large', ncol=2, bbox_to_anchor=(-0.1, .5), bbox_transform=gcf().transFigure, fontsize='large', frameon=False)
    salt_legend = []
    for label in parameters['Salinity'].sort_values().unique():
        salt_patch = figure2.ax_row_dendrogram.bar(0, 0, color = stat_salt[label], label = label, linewidth=0)
        salt_legend.append(salt_patch)
    l3 = plt.legend(salt_legend,parameters['Salinity'].sort_values().unique(), loc="lower left", title = 'Salinity',title_fontsize= 'large', ncol=2, bbox_to_anchor=(-0.1, .25), bbox_transform=gcf().transFigure, fontsize='large', frameon=False)
    plt.gca().add_artist(l2)
    plt.setp(figure1.ax_heatmap.yaxis.get_majorticklabels())
    plt.setp(figure1.ax_heatmap.xaxis.get_majorticklabels())
    svg_name=name +str(vers)+'.svg'
    svg_file=f'{visuals}/{svg_name}'
    png_name=name +str(vers)+'.png'
    png_file=f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format = 'svg',dpi=gcf().dpi, bbox_inches='tight')
        plt.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()  # Push new figure on stack
    station_order = coverage.index.values.tolist()
    station_reorder = figure2.dendrogram_row.reordered_ind
    return station_order, station_reorder,cluster_st_df, cluster_pl_df

def Pearson_cor(df1,df2,df_to_update):
    ''' The function uses stats.pearsonr to calculate Pearson correlation between two vectors'''
    for index_o, row_o in df1.iterrows():
        # iterating clusters and getting coverage of the index_o cluster at each row_o
        for index_p, row_p in df2.iterrows():
            # iterating physical parameters and getting parameter of the index_p environmental condition at each row_p sampling point
            print("Pearson correlation, p-value for Cluster %s and %s" % (index_o, index_p))
            correlation, p_value= stats.pearsonr(row_o, row_p) # calculating Pearson correlation and p-value for cluster:env.condition in each sampling point
            print(round(correlation,3),round(p_value,3))
            df_to_update[index_p][index_o] = (round(correlation,3),round(p_value,3)) # updating df_pearson dataframe at with calculated correlation values
    return df_to_update

def Correlation_calculation(class_df, cand_df, name):
    ''' The function calculates correlation between plasmidome clusters and environmental conditions '''
    station_order, station_reorder, cluster_st_df, cluster_pl_df = class_df
    # getting plasmid candidates clusters
    cluster_pl_df = cluster_pl_df.reset_index()
    cluster_pl_df['Plasmid candidates clusters'] = 'C' + cluster_pl_df['Plasmid candidates clusters'].astype(str)
    # getting sampling points clusters
    cluster_st_df = cluster_st_df.reset_index(drop=True)
    cluster_st_df['Sampling points clusters'] = cluster_st_df['Sampling points clusters'].astype(int).apply(lambda x: chr(ord('`')+x))
    cluster_st_df['Sampling points clusters'] = 'C' + cluster_st_df['Sampling points clusters']
    print(cluster_st_df)
    #getting physical parameters dataframes
    df_phys, df_phys_corr = Physical(3)
    df_phys = df_phys[['St_Depth', 'Latitude', 'Salinity', 'Chlorophyll', 'Turbidity', 'Temp.', 'Oxygen', 'Depth', 'Nitrate', 'Phosphate', 'Silicate']]
    df_phys = df_phys.set_index('St_Depth')
    df_phys = df_phys.T # physical parameters: y=parameter; x=sampling point
    df_phys_corr = df_phys_corr.corr() # physical parameters cross-correlation
    col_order = df_phys.columns.tolist()
    ### getting coverage of plasmid candidates in each sampling point
    reads_df = CoverageDF(cand_df['Plasmid'].unique())[1]
    reads_df = reads_df.reset_index()
    reads_df = reads_df.melt(id_vars = ['rname'], var_name = 'Station', value_name = 'coverage')
    station = Station()
    reads_stat = reads_df.merge(station, left_on = 'Station', right_on = 'station_name')
    reads_stat.drop(['station_name', 'Sample'], axis = 1, inplace = True)
    reads_df_st = pd.pivot(reads_stat, values = 'coverage', index = 'rname', columns = 'St_Depth')
    ### getting coverage of plasmid candidates clusters in each sampling point
    out_pl = cluster_pl_df.merge(reads_df_st, left_on = 'Plasmids', right_on = 'rname')
    out_pl.sort_values('Plasmid candidates clusters', inplace = True)
    out_group_pl =out_pl.groupby('Plasmid candidates clusters').mean() # calculating coverage average for each plasmid candidates cluster at each sampling point
    out_group_pl = out_group_pl.reindex(col_order, axis=1) # working dataframe with y=plasmid candidates cluster, x=sampling point, values=coverage.mean
    df_pearson_pl = pd.DataFrame(columns = df_phys.index.to_list(), index = out_group_pl.index.to_list()) # empty dataframe for pearson clusters:env.conditions
    ### getting Pearson correlation for each cluster-env.condition
    Pearson_cor(out_group_pl, df_phys, df_pearson_pl)
    #df.assign(**df[['col2', 'col3']].apply(lambda x: x.str[0]))
    df_pearson_2 = df_pearson_pl.assign(**df_pearson_pl[df_pearson_pl.columns.to_list()].apply(lambda x: x.str[0]))
    # correlation dfs
    df_pl_cluster_corr = out_group_pl.T
    df_pl_cluster_corr = df_pl_cluster_corr.corr()
    df_pl_cluster_corr_P = pd.concat([df_pl_cluster_corr, df_pearson_2], axis=1) #correlation matrix C1-8:C1-8,envs
    df_pearson_2_trans = df_pearson_2.T
    df_pears = pd.concat([df_pearson_2_trans, df_phys_corr], axis=1) #correlation matrix C1-8,envs:C1-8,envs
    df_Pearson = df_pl_cluster_corr_P.append(df_pears)
    ### plotting correlation matrix
    plt.figure(figsize = (15,12))
    sns.set(font_scale = 1.2)
    ax = sns.heatmap(df_Pearson, cmap='coolwarm', vmin=-1, vmax=1, annot=True, cbar_kws = {"ticks": [-1, 1]})
    svg_name = 'heatmap_corr_' + name + str(4) + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'heatmap_corr_' + name + str(4) + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format='svg', dpi=gcf().dpi, bbox_inches='tight')
        plt.savefig(png_file, format='png', dpi=gcf().dpi, bbox_inches='tight')
    plt.show()
    path_file = f'{tables}/Pearson.xlsx'
    with pd.ExcelWriter(path_file, engine="openpyxl", mode = 'a') as writer:
        df_Pearson.to_excel(writer, sheet_name = name)




#Correlation_calculation(Clust_map2(7,Plasmid_class()[2],'PlPutUnc_HMannot_', 1150, 1200),Plasmid_class()[2], 'All')
#Correlation_calculation(Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900),Plasmid_class()[1], 'PlPut')
#Correlation_calculation(Clust_map2(4,Plasmid_class()[0],'Pl_HMannot_', 250, 400),Plasmid_class()[0], 'Pl')
#print(Plasmid_class()[0]['Plasmid'].unique())
#order_pl=Clust_map2(4,Plasmid_class()[0],'Pl_HMannot_', 250, 400)[2]
#order_all=Clust_map2(4,Plasmid_class()[2],'PlPutUnc_HMannot_', 1150, 1500)[2]

#print(Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900)[2])
#Clust_map2(4,Plasmid_class()[0],'Pl_HMannot_', 250, 400)
#Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900)
#Clust_map2(4,Plasmid_class()[2],'PlPutUnc_HMannot_', 1150, 1500)


