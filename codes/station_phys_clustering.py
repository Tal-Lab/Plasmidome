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
import scipy.cluster.hierarchy as sch
from collections import defaultdict
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import os
from general_analysis import plasmids_by_reads, tables, Station, GetLibSize
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
    data = data.rename_axis("Sampling stations", axis="columns")
    # transpose data to fit it in clustermap with stations on the y axis and candidates on the x axis
    data_transp = data.transpose()
    return data_transp

def Physical(vers):
    df = pd.read_excel(physical, index_col = 0, header = 0)
    df = df.dropna()
    df_for_corr = df.reset_index().drop(['Sample', 'Station'], axis=1)
    df_for_corr = df_for_corr.apply(pd.to_numeric)
    mask = np.triu(np.ones_like(df_for_corr.corr(), dtype=bool))
    f, ax = plt.subplots(figsize=(9, 6))
    sns.heatmap(df_for_corr.corr(), mask=mask,annot_kws={"size": 10},
                  fmt='.2f', vmin=-1, vmax=1, annot=True,cmap='coolwarm')
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
    df_phys = df[['St_Depth','Latitude','Depth','Temp.','Temperature', 'Salinity', 'Oxygen', 'Nitrate', 'Phosphate', 'Silicate']]
    #print(df_phys)
    return df_phys
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
    print(number_more_stat)
    print(df_one_stat['NewName'].unique())
    number_all = df['NewName'].nunique()
    perc_one = (number_one_stat/number_all)*100
    perc_more = (number_more_stat/number_all)*100
    print(perc_one)
    print(perc_more)
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
    parameters = Physical(1)
    parameters = parameters.set_index(parameters['St_Depth'])
    # maybe put parameters in loop?
    parameters["Temperature"] = parameters['Temp.'].astype(float).round().astype(int)
    stat_temp = dict(zip(parameters["Temperature"].sort_values().unique(), sns.color_palette("Reds", parameters['Temperature'].nunique())))
    temperature = parameters["Temperature"].map(stat_temp)
    stat_depth = dict(zip(parameters['Depth'].sort_values().unique(), sns.color_palette("PuBu", parameters['Depth'].nunique())))
    depth = parameters['Depth'].map(stat_depth)
    parameters["Salinity"] = parameters['Salinity'].astype(float).round().astype(int)
    stat_salt = dict(zip(parameters['Salinity'].sort_values().unique(), sns.color_palette("Greens", parameters['Salinity'].nunique())))
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
    cluster_st_df = pd.DataFrame(rows, columns=["St_Depth", "Cluster"])
    cluster_st_df = cluster_st_df.set_index(parameters['St_Depth'])
    stat_cluster = dict(
        zip(cluster_st_df['Cluster'].unique(), sns.color_palette("Pastel1", cluster_st_df['Cluster'].nunique())))
    cluster_st = cluster_st_df['Cluster'].map(stat_cluster)
    stations = cluster_st_df.index.values.tolist()
    
    # plasmid clusters indices correspond to indices of original df

    L_pl = figure1.dendrogram_col.linkage
    clusters_pl = sch.fcluster(L_pl, pl, 'distance')
    cols = []
    for i, cluster in enumerate(clusters_pl):
        cols.append([coverage.columns[i], cluster])
    cluster_pl_df = pd.DataFrame(cols, columns = ["Plasmids", "Cluster"])
    cluster_pl_df = cluster_pl_df.set_index('Plasmids')
    plas_cluster = dict(
        zip(cluster_pl_df['Cluster'].unique(), sns.color_palette("Greys", cluster_pl_df['Cluster'].nunique())))
    cluster_pl = cluster_pl_df['Cluster'].map(plas_cluster)
    plasmids = cluster_pl_df.index.values.tolist()

    empty = 0*len(stations)
    for_df = {'St_Depth':stations, 'Empty': empty}
    space_df = pd.DataFrame(for_df)
    space_df=space_df.set_index(parameters['St_Depth'])
    space_df.columns = ['St_Depth', ' ']
    stat_space = dict(zip(space_df[' '].unique(), "white"))
    space_st = space_df[' '].map(stat_space)
    row_colors = pd.concat([cluster_st, salt, temperature, depth, space_st], axis = 1)
    figure2 = sns.clustermap(data = coverage,
                             metric = "euclidean",
                             method = 'ward',
                             row_colors = row_colors,
                             col_colors=cluster_pl,
                             cmap = sns.color_palette("Blues", as_cmap = True),
                             linewidths = 0.0,
                             xticklabels = False,
                             yticklabels = True,
                             rasterized = True,
                             cbar_kws = {"ticks": [0, 100]}
                             )
    # g.set_axis_labels(["Plasmid candidates", "Sampling stations"])
    figure2.fig.subplots_adjust(left = -.01, right = 0.8)
    figure2.cax.set_title("AP")
    figure2.ax_cbar.set_position((0.92, .16, .03, .6))
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
        temp_patch = mpatches.Patch(color = stat_temp[label], label = label)
        tem_legend.append(temp_patch)
    l1 = plt.legend(title = 'Temperature', handles = tem_legend, bbox_to_anchor = (-25, 1.13), loc = "upper right",
                    borderaxespad = 2.0)
    plt.gca().add_artist(l1)
    # depth legend
    # loop?
    dep_legend = []
    for label in parameters['Depth'].sort_values().unique():
        dep_patch = mpatches.Patch(color = stat_depth[label], label = label)
        dep_legend.append(dep_patch)
    l2 = plt.legend(title = 'Depth', handles = dep_legend, bbox_to_anchor = (-25, 0.45), loc = "right",
                    borderaxespad = 2.0)
    salt_legend = []
    for label in parameters['Salinity'].sort_values().unique():
        salt_patch = mpatches.Patch(color = stat_salt[label], label = label)
        salt_legend.append(salt_patch)
    l3 = plt.legend(title = 'Salinity', handles = salt_legend, bbox_to_anchor = (-25, -.11), loc = "lower right",
                    borderaxespad = 2.0)
    plt.gca().add_artist(l2)
    plt.setp(figure1.ax_heatmap.yaxis.get_majorticklabels(), fontsize = 14)
    plt.setp(figure1.ax_heatmap.xaxis.get_majorticklabels(), fontsize = 14)
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


    # calculate correlation vector
    parameters["Temp."] = parameters['Temp.'].astype(float)
    parameters['Depth'] = parameters['Depth'].astype(int)
    parameters['Latitude'] = parameters['Latitude'].astype(float)
    df_correlation = parameters[['Temp.', 'Depth', 'Salinity', 'Oxygen', 'Nitrate', 'Phosphate', 'Silicate', 'Latitude']].merge(
        cluster_st, left_index=True, right_index=True)
    Cluster_array = df_correlation['Cluster'].to_numpy()
    param_arrays = [df_correlation['Temp.'], df_correlation['Depth'], df_correlation['Salinity'],
                    df_correlation['Oxygen'], df_correlation['Latitude'],
                    df_correlation['Nitrate'], df_correlation['Phosphate'], df_correlation['Silicate']]

    '''
    #pearson_dict = {}
    for param_col in param_arrays:
        param_array = param_col.to_numpy()
        model = LinearRegression.fit(Cluster_array, param_array)
        print(LinearRegression.score(Cluster_array, param_array))
        #pearson_param = stats.pearsonr(Cluster_array, param_array)
        #entry = {'r': pearson_param[0].round(3), 'p-value': pearson_param[1].round(3)}
        #pearson_dict[param_col.name] = entry
    #print(pearson_dict)
    #df_params_pearson = pd.DataFrame(pearson_dict).T
    #df_params_pearson.fillna(0, inplace=True)
    #print(df_params_pearson)
    
    # df_correlation = parameters['Depth'].merge(cluster_df,left_index=True, right_index=True)

    slope_t, intercept_t, r_value_t, p_value_t, std_err_t = stats.linregress(df_correlation["Cluster"],
                                                                             df_correlation["Temp."])
    print("Clusters:Temperature,C: slope=%f; intercept=%f; r_value=%f; p_value=%f, str_err=%f" %
          (slope_t.round(3), intercept_t.round(3), r_value_t.round(3), p_value_t.round(3), std_err_t.round(3)))
    slope_d, intercept_d, r_value_d, p_value_d, std_err_d = stats.linregress(df_correlation["Cluster"],
                                                                             df_correlation["Depth"])
    print("Clusters:Depth,m: slope=%f; intercept=%f; r_value=%f; p_value=%f, str_err=%f" % (
        slope_d.round(3), intercept_d.round(3), r_value_d.round(3), p_value_d.round(3), std_err_d.round(3)))
    slope_l, intercept_l, r_value_l, p_value_l, std_err_l = stats.linregress(df_correlation["Cluster"],
                                                                             df_correlation["Latitude"])
    print("Clusters:Latitude,C: slope=%f; intercept=%f; r_value=%f; p_value=%f, str_err=%f" %
          (slope_l.round(3), intercept_l.round(3), r_value_l.round(3), p_value_l.round(3), std_err_l.round(3)))
    '''
    return station_order, station_reorder,cluster_st_df, cluster_pl_df
order_pl=Clust_map2(4,Plasmid_class()[0],'Pl_HMannot_', 250, 400)[2]
order_plput = Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900)[2]
order_all=Clust_map2(4,Plasmid_class()[2],'PlPutUnc_HMannot_', 1150, 1500)[2]

#print(Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900)[2])
#Clust_map2(4,Plasmid_class()[0],'Pl_HMannot_', 250, 400)
#Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900)
#Clust_map2(4,Plasmid_class()[2],'PlPutUnc_HMannot_', 1150, 1500)


