# -*- coding: utf-8 -*-
"""
Created on 19/12/2022 13:56

Author: Lucy

Description: Plasmidome Detection pipeline segments IV-V
"""

import re, os
import pandas as pd
from itertools import combinations
import scipy.cluster.hierarchy as ch
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
from pathlib import Path

### version - change each time - remove irrelevant
version = 1

# dataframe display options (optional)
np.set_printoptions(threshold=np.inf)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option("mode.chained_assignment", None)

# ncbi output structure
catcols = ['qseqid', 'sseqid',  'evalue', 'length', 'pident', 'score', 'qcovs', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq']

### directions
## working directories
visuals = r"../visualisations"
tables = r"../data_calculations"
resource = r"../res"
Path(visuals).mkdir(parents=True, exist_ok=True)
Path(tables).mkdir(parents=True, exist_ok=True)

## working files
# blast all_vs_all output - put it in the data_calculations folder
file = f'{resource}/all_vs_all.csv'
# stations file - needs preparing and should be put in res folder - reqiures all parametes: sample_id, Station number, station location, temperature, depth, any other physical parameters
stations = f'{resource}/stations.txt'

## files for output
plasmids_by_clusters = f'{tables}/Plasmids_byClusters.csv'
heatmap_png = f'{visuals}/HeatmapAP.png'
heatmap_svg = f'{visuals}/HeatmapAP.svg'
# dataframe reader
all_vs_all = pd.read_csv(file, sep = '\t', index_col = None, header = None)
all_vs_all.columns = catcols

def Clean_length(node_name):
    # getting length of each predicted plasmid (node) from its name by pattern
    if re.search('h\_\d+', node_name):
        pos = int(re.search('h\_\d+', node_name).group(0)[2:])
        return pos
    else:
        return node_name

def Clean_Stat(node_name):
    # getting station of each predicted plasmid (node) from its name by pattern
    if re.search("^\d+", node_name):
        pos = str(re.search("^\d+", node_name).group(0)[:])
        return pos
    else:
        return node_name

def Shortname(node_name):
    # getting shorter name of each predicted plasmid
    if re.search('^\d+\w+c', node_name):
        pos = re.search('^\d+\w+c', node_name).group(0)[:-2]
        return pos
    else:
        return node_name

def Get_length(df):
    # extracting length of query and subject sequences
    df['qseqid'] = df['qseqid'].apply(Shortname)
    df['sseqid'] = df['sseqid'].apply(Shortname)
    df['length'] = pd.to_numeric(df['length'], errors = 'ignore')
    df['q_len'] = df['qseqid'].apply(Clean_length)
    df['s_len'] = df['sseqid'].apply(Clean_length)
    df = df.loc[df['q_len'] >= 1000]
    df = df.loc[df['s_len'] >= 1000]
    return df

def overlap( r1, r2 ):
    # calculating overlapping part of the hits
    if r1[0] < r1[1]:
        left = max(r1[0], r2[0])
        right = min(r1[1], r2[1])
        over = right - left
    else:
        left = min(r1[0], r2[0])
        right = max(r1[1], r2[1])
        over = left - right
    return over if over>0 else None

def Pairs(list):
    # getting list of overlap of each hit
    a = []
    for c1, c2 in combinations(list, 2):
        o = overlap(c1, c2)
        if o is not None:
            a.append(o)
    return a

def OverlapSum(l):
    # calculating sum of overlaps
    sum = 0
    for x in l:
        sum += x
    return sum

def Interval(df, x, y):
    # calculates overlaps of specified hit
    df_int = df[[x,y]]
    interval_list = df_int.values.tolist()
    ol_list = Pairs(interval_list)
    sum = OverlapSum(ol_list)
    return int(sum)

def Min_Max(df, x, y):
    # choosing which start and end value to pass
    if df[x].iloc[0] < df[y].iloc[0]:
        z = 'min'
    else:
        z = 'max'
    return z

def DF_size(df,purp):
    # calculating alignment % for each pair
    r,c = df.shape
    total = df['length'].sum()
    q_overlap = Interval(df, 'qstart', 'qend')
    s_overlap = Interval(df, 'sstart', 'send')
    s = total - s_overlap
    q = total - q_overlap
    if df['qseqid'].equals(df['sseqid']):
        df.loc[:, 'align_perc'] = 1
    else:
        if purp == "cluster":
            if r == 1:
                df.loc[:, 'align_perc'] = df['length']/df[['q_len','s_len']].min().min()
            else:
                df.loc[:, 'align_perc'] = min(q,s)/ df[['q_len','s_len']].min().min()
        elif purp == "heatmap":
            if r == 1:
                df.loc[:, 'align_perc'] = df['length'] / df['q_len']
            else:
                df.loc[:, 'align_perc'] = min(q,s) / df['q_len']
        else:
            print("The purpose of function is not known")
    return df

def CleanFiltered(df_name,purp):
    # cleaning data from additional hits
    df = Get_length(df_name)
    grouped = df.groupby(["qseqid", "sseqid"]) # grouping of dataframe by query and subject, to see number of hits in each pair
    keys = grouped.groups.keys()
    appended_df = []
    # iteration over groups
    for name, group in grouped:
        df_new = DF_size(group, purp)
        appended_df.append(df_new)
    appended_df = pd.concat(appended_df)
    to_print = appended_df[['qseqid', 'sseqid', 'evalue', 'pident', 'length',
                            'qstart', 'qend', 'sstart', 'send', 'q_len', 's_len', 'align_perc']]
    print(to_print)
    to_print.to_excel(f'{tables}/Plasmids_AP.xlsx', header = True)
    return appended_df

def GetHits(purp):
    # calculating number of repeating plasmids in each plasmid-group
    df = CleanFiltered(all_vs_all,purp)
    df.loc[:,'hits']=df.groupby(['qseqid'])['qseqid'].transform('count')
    return df

def DFtoPivot(purp):
    df = GetHits(purp)
    # converting df to pivot
    data_pivoted = pd.pivot_table(df, index = ['qseqid'], columns = ['sseqid'], values = ['align_perc'],
                                  fill_value = 0)
    return data_pivoted

def PivotToArray(purp):
    df = DFtoPivot(purp)
    # getting alignment percentage as array
    arrAP = df.values
    return arrAP

def Heatmap(arr):
    val = PivotToArray("heatmap")
    linkage = ch.linkage(arr, 'ward')
    figure1 = sns.clustermap(val,
                             row_linkage = linkage,
                             col_linkage = linkage,
                             linewidths=0,
                             cmap = sns.color_palette("Blues", as_cmap = True),
                             xticklabels=False,
                             yticklabels=False,
                             rasterized = True,
                             cbar_pos=(0.1, .2, .03, .45),
                             cbar_kws = {'label':'AP'})
    figure1.ax_col_dendrogram.remove()
    figure1.ax_row_dendrogram.remove()
    #plt.xlabel('Plasmid Candidates')
    #plt.ylabel('Plasmid Candidates')
    if not os.path.isfile(heatmap_png) and not os.path.isfile(heatmap_svg):
        plt.savefig(heatmap_svg, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        plt.savefig(heatmap_png, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    return linkage

def CandidatesClusters():
    clustering = Heatmap(PivotToArray("cluster"))
    df = DFtoPivot("cluster")
    clusters = ch.fcluster(clustering, t = 1.4, criterion = 'distance')
    clusters = list(clusters)
    labels = list(df.index)
    label_num = list(enumerate(labels))
    new_label_num = list(map((lambda x, y: x + (y,)), label_num, clusters))
    clusters_df = pd.DataFrame(new_label_num, columns = ['Number', 'Plasmid', 'Cluster'])
    clusters_df = clusters_df.drop('Number',axis=1)
    clusters_df.sort_values(by = 'Cluster', ascending = True, inplace = True)
    clusters_df['Count'] = clusters_df.groupby('Cluster')['Cluster'].transform('count')
    clusters_df['Station'] = clusters_df['Plasmid'].apply(Clean_Stat)
    clusters_df['Length'] = clusters_df['Plasmid'].apply(Clean_length)
    if not os.path.isfile(plasmids_by_clusters):
        clusters_df.to_csv(plasmids_by_clusters, index = False)
    return new_label_num

def Histo():
    df = CleanFiltered(all_vs_all,"cluster")
    df = df[['qseqid', 'q_len']]
    df.drop_duplicates(subset = "qseqid",
                         keep = False, inplace = True)
    df['Length_log'] = np.log10(df['q_len'])
    print(df['qseqid'].nunique())
    fig1 = sns.histplot(df, x = "Length_log", stat ='count', bins = 20)
    fig1.set_xlabel("Log{}(length)".format('\u2081\u2080'))
    #fig1.set_xlabel("Length, bp")
    png_name = "HistoLengthCount_" + str(version) + ".png"
    svg_name = "HistoLengthCount_" + str(version) + ".svg"
    plt.savefig(svg_name, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_name, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()

CandidatesClusters()
Histo()