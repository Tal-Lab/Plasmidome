import pandas as pd
import re, os
from Bio import SeqIO
import numpy as np
#import plotly.express as px
#import plotly.graph_objects as go
from sklearn import preprocessing
from scipy.cluster.hierarchy import linkage, fcluster
#import dash_core_components as dcc
#import dash_bio as dashbio
import seaborn as sns
import colorcet as cc
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
import scipy.cluster.hierarchy as sch
from collections import defaultdict
import matplotlib.patches as mpatches
from station_phys_clustering import Clust_map2, Plasmid_class
### Description
# add description

# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

# working directories
out_dir = f"{path}/data_calculations"
visuals = f"{path}/visualisations"

# working files
combined_output = r"../res/unique_plasmids.fasta"
proteins = r"../res/Filtered_ORFs.fasta"
cog_categories = r"../res/Annotations/cog_cats.csv"
eggnog = r"../res/eggnog_FilteredORFs.csv"
stations = r"../res/stations.txt"
reads_coverage = r"../res/all_cov.csv"
plasmids_byreads = f"{out_dir}/AllPlasmids_perStat.csv"
physical = r"../res/station_physical.xlsx"

station_orderCl, station_reorderCl, cluster_st_df, cluster_pl_df = Clust_map2(4,Plasmid_class()[2],'PlPutUnc_HMannot_', 1150, 1500)
#station_orderPlPut, station_reorderPlPut, cluster_st_dfPlPut, cluster_pl_dfPlPut = Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900)
#station_order7, station_reorder7,cluster_st_df7, cluster_pl_df7 = Clust_map2(4,Plasmid_class()[0],'Pl_HMannot_', 250, 400)

def csv_reader(file):
    df = pd.read_csv(file, sep = ',', header = 0, index_col = None)
    df = df.rename(columns={"#query": "Query","COG_category": "COG cat"})
    df['COG cat'].fillna('', inplace = True)
    #print(df.columns)
    return df

def SPLIT(word):
    word_split = list(word)
    #print(word_split)
    return word_split

def SplitColumn():
    dd = csv_reader(eggnog)
    dd['COG cat'] = dd['COG cat'].apply(SPLIT)
    subset = list(dd.columns.difference(['COG cat']))
    # convert list of pd.Series then stack it
    dd = (dd
          .set_index(subset)['COG cat']
          .apply(pd.Series)
          .stack()
          .reset_index()
          .rename(columns = {0: 'COG cat'}))
    return dd

def DefineStation():
    df = pd.read_csv(plasmids_byreads, index_col=None, header = 0)
    df.dropna(inplace = True)
    return df

def Station():
    with open(stations) as adressFile:
        matrix = np.loadtxt(adressFile, dtype = "str")
    df = pd.DataFrame(data = matrix, index = None, columns = matrix[0, 0:])
    df.reset_index()
    df.drop(index = 0, inplace = True)
    df['station_name'] = df['Sample'] + "_coverage"
    return df

def GroupPlaceMapper():
    df1 = DefineStation()
    df2 = Station()
    df1['St_Depth'] = df1['station_name'].map(df2.set_index('station_name')['St_Depth'])
    #df1['station_name'] = df1['Station'] + "_coverage"
    df1['Depth'] = df1['St_Depth'].apply(lambda x: (int(re.search("_\d+", x).group(0)[1:])))
    return df1

def MapToFunc():
    df_plasm = SplitColumn()[['Query', 'COG cat']]
    df_func = csv_reader(cog_categories)
    df_plasm['Functional categories'] = df_plasm['COG cat'].map(df_func.set_index('COG cat')['Functional categories'])
    #print(df_plasm)
    return df_plasm

def MapToPlace():
    df_plasm = MapToFunc()
    df_plasm['qseqid'] = df_plasm['Query'].apply(lambda x: re.search('^.*\_', x).group(0)[:-1])
    df_stations = GroupPlaceMapper()[['NewName', 'St_Depth']]
    new_df = pd.merge(df_plasm, df_stations, left_on = 'qseqid', right_on = 'NewName')
    #print(new_df)
    return new_df

def FreqFuncStat(freq, unknown):
    df = MapToPlace()[['St_Depth', 'COG cat', 'Functional categories']]
    df2 = df[['COG cat', 'Functional categories']].drop_duplicates().fillna('missing')
    df3 = GroupPlaceMapper()[['station_name', 'St_Depth']].drop_duplicates()
    if unknown != 'with':
        indexNames = df[(df['COG cat'] == 'S') | (df['COG cat'] == '-')].index
        df.drop(indexNames, inplace=True)
    df_grouped = df.groupby('St_Depth')['COG cat'].value_counts(normalize = freq).to_frame(name='Function count')
    df_grouped = df_grouped.reset_index()
    df_grouped['Functional categories'] = df_grouped['COG cat'].map(df2.set_index('COG cat')['Functional categories'])
    df_grouped['Functional categories'] = df_grouped[['COG cat', 'Functional categories']].apply(tuple, axis = 1)
    df_grouped['Functional categories'] = df_grouped['Functional categories'].apply(lambda x: ' - '.join(x))
    df_grouped['station_name'] = df_grouped['St_Depth'].map(df3.set_index('St_Depth')['station_name'])
    print(df_grouped)
    return df_grouped

def Station_Order(order):
    init_order = station_orderCl
    #print(init_order)
    counter_list = list(enumerate(init_order, 0))
    #print("prinying counter list")
    #print(counter_list)
    station_order = station_reorderCl
    #print("Printing station order")
    #print(station_order)
    init_df = pd.DataFrame(counter_list, columns = ['Number', 'Sampling Station']).set_index('Number').reindex(station_order)
    #print(init_df)
    final_order = init_df['Sampling Station'].to_list()
    #print(final_order)
    return final_order

def BarChart(freq, unknown):
    'plotting stacked barchart for COGs frequencies in all plasmid candidates by sampling points'
    df = FreqFuncStat(freq, unknown)
    #reorder stations by clusters
    order_Stat = Station_Order(True)
    true_sort1 = [s for s in order_Stat if s in df.St_Depth.unique()]
    print(true_sort1)
    df = df.set_index('St_Depth').loc[true_sort1].reset_index()
    station = df['St_Depth'].to_list()
    ncolors = df['Functional categories'].nunique()
    colors = sns.color_palette(cc.glasbey, n_colors = ncolors)
    labels = df['Functional categories'].unique()
    #prepare figure
    sns.set_theme()  # to make style changable from defaults use this line of code befor using set_style
    plt.figure(figsize = (8, 8))
    sns.set(font_scale = 0.95
            )
    with sns.axes_style("ticks"):
        fig = sns.histplot(df,
                           y=station,
                           weights='Function count',
                           hue='Functional categories',
                           multiple='stack',
                           hue_order = labels[::-1],
                           palette=colors,
                           # Add white borders to the bars.
                           edgecolor='white')
        #fig.set(ylabel='Sampling points', xlabel='Function frequency')
        fig.set_xlabel('Function frequency', fontsize = 'medium')
        fig.set_ylabel('Sampling points', fontsize = 'medium')
        plt.margins(0,0)
        #fig.set_style("white")
        #fig.yaxis.set_label_position("right")
        #fig.yaxis.set_ticks_position("right")
        # Put the legend out of the figure
        fig.legend(labels, title = 'Functional categories (COGs)', bbox_to_anchor = (1.01, 1), ncol = 1, title_fontsize = 'medium',
               fontsize='medium', frameon=False, loc = 2, borderaxespad = 0.)
    #fig.tick_params(axis = 'x', rotation = 90, labelsize = 8)
    #save graph in PNG and vector format
    svg_name = 'barplot_COG_' + unknown + str(3) + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'barplot_COG_' + unknown + str(3) + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format='svg', dpi=gcf().dpi, bbox_inches='tight')
        plt.savefig(png_file, format='png', dpi=gcf().dpi, bbox_inches='tight')
    plt.show()

def BarChart_lim(df_class, cluster, file_name, unknown):
    'plotting stacked barchart for COGs frequencies in plasmids/plasmids+putative plasmids'
    df = MapToFunc()
    df['Plasmid'] = df['Query'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    plasmid_list = df_class['Plasmid'].unique()
    #preparing to reorder plasmids by cluster
    cluster.reset_index(inplace=True)
    cluster.sort_values(by='Plasmid candidates clusters', inplace=True)
    order_pl = cluster['Plasmids'].to_list()

    df_trunc = df.loc[df['Plasmid'].apply(lambda x: x in plasmid_list)]
    df_trunc = df_trunc[['Query', 'Plasmid','COG cat', 'Functional categories']].fillna('missing')
    if unknown != 'with':
        indexNames = df_trunc[(df_trunc['COG cat'] == 'S') | (df_trunc['COG cat'] == '-')].index
        df_trunc.drop(indexNames, inplace=True)
    df_trunc['Functional categories'] = df_trunc[['COG cat', 'Functional categories']].apply(tuple, axis=1)
    df_trunc['Functional categories'] = df_trunc['Functional categories'].apply(lambda x: ' - '.join(x))
    df_grouped = df_trunc.groupby('Plasmid')['Functional categories'].value_counts(normalize=True).to_frame(name='Function count')
    df_grouped = df_grouped.reset_index()
    #reordering plasmids by cluster order
    true_sort1 = [s for s in order_pl if s in df_grouped.Plasmid.unique()]
    df_grouped = df_grouped.set_index('Plasmid').loc[true_sort1].reset_index()
    print(df_grouped['Plasmid'].unique())
    ncolors = df_grouped['Functional categories'].nunique()
    colors = sns.color_palette(cc.glasbey, n_colors=ncolors)
    labels = df_grouped['Functional categories'].unique()
    #figure
    with sns.axes_style("ticks"):
        fig = sns.histplot(df_grouped,
                           x='Plasmid',
                           weights='Function count',
                           hue='Functional categories',
                           hue_order = labels[::-1],
                           multiple='stack',
                           #palette= {'bright',ncolors},
                           palette = colors,
                           # Add white borders to the bars.
                           edgecolor='white',
                           # Shrink the bars a bit so they don't touch.
                           shrink=0.8
                           )
        fig.set(xlabel='Plasmid', ylabel='Function frequency')
        fig.margins(x = 0.1)
        # Put the legend out of the figure
        fig.legend(labels, title='Functional categories (COGs)', bbox_to_anchor = (1.1, 1), ncol = 1, title_fontsize = 16, loc = 2, borderaxespad = 0.)
        fig.tick_params(axis='x', rotation=90, labelsize=8)
    #saving figure
    svg_name = file_name + unknown + str(1) + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = file_name + unknown + str(1) + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        plt.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()

def ProteinsFastaAn():
    records = []
    for record in SeqIO.parse(proteins, "fasta"):
        records.append(record.id)
    plasmids_w_proteins = [(re.search('^.*\_', x).group(0)[:-1]) for x in records]
    plasmids = list(set(plasmids_w_proteins))
    data = {'Plasmids' : plasmids_w_proteins, 'Proteins' : records}
    df = pd.DataFrame.from_dict(data)
    df_grouped = df.groupby('Plasmids')['Proteins'].nunique().to_frame(name = 'Number of Proteins').reset_index()
    min = df_grouped['Number of Proteins'].min()
    max = df_grouped['Number of Proteins'].max()
    mean = df_grouped['Number of Proteins'].mean()
    print("******************** Number of ORFs in plasmids ranges %s - %s *****************" % (str(min), str(max)))
    print("******************* The mean number of ORFs in plasmid is %s *************" % str(mean))
    #fig = px.histogram(df_grouped, x = "Number of Proteins", histnorm = 'percent', nbins = 50)
    #fig.show()
    #fig.write_html("ProteinsPercHisto.html")
    return df

def eggNOGStats():
    df1 = ProteinsFastaAn()
    prodigal = len(df1['Proteins'].to_list())
    print("Number of proteins predicted with PRODIGAL: "+str(prodigal))
    df2 = csv_reader(eggnog)[['Query','COG cat']]
    eggnog_n = len(df2['Query'].to_list())
    print("Number of proteins detected with eggNOG: "+str(eggnog_n))
    df1['COG cat'] = df1['Proteins'].map(df2.set_index('Query')['COG cat']).replace('-', 'missing').fillna('missing')
    #df1['COG cat'] = df1['COG cat']
    print(df1['COG cat'].unique())
    eggnog_perc = eggnog_n/prodigal*100
    formatted_eggnog = "{:.2f}".format(eggnog_perc)
    print("eggNOG mapped only " + str(formatted_eggnog) + "% of predicted genes")
    ### this part is only for detecting proteins with assigned COG categories
    df2['COG cat'] = df2['COG cat'].replace('-', 'missing').fillna('missing')
    df2 = df2[df2['COG cat'] != 'missing']
    print('NUmber of proteins assigned any COG category: %d' % df2['Query'].nunique())
    ###
    plasmids_w_proteins_EN = [(re.search('^.*\_', x).group(0)[:-1]) for x in df2['Query'].to_list()]
    plasmids_EN = list(set(plasmids_w_proteins_EN))
    print("Number of plasmids detected with eggNOG COG cats: "+str(len(plasmids_EN)))
    df3 = df2[df2['COG cat'] == 'S']
    plasmids_w_proteins_S = [(re.search('^.*\_', x).group(0)[:-1]) for x in df3['Query'].to_list()]
    plasmids_S = list(set(plasmids_w_proteins_S))
    print("Number of plasmids detected with eggNOG COG category S - Function unknown: " + str(len(plasmids_S)))
    df4 = df2[df2['COG cat'] != 'S']
    print("Number of proteins detected with eggNOG COG category : " + str(df4['Query'].nunique()))
    eggnog_all = len(df2['Query'].to_list())
    eggnog_S = len(df3['Query'].to_list())
    eggnog_S_perc = eggnog_S / eggnog_all * 100
    formatted_eggnog_S = "{:.2f}".format(eggnog_S_perc)
    print("Number of proteins detected with eggNOG and assigned to COG category S - Function unknown: " + str(eggnog_S))
    print(str(formatted_eggnog_S) + "% of predicted genes are S - Function unknown")
    plasmids_w_proteins_all = df1.Plasmids.unique()
    missing_plasmids = [x for x in plasmids_w_proteins_all if x not in plasmids_EN]
    print(str(len(missing_plasmids)) + " plasmids, whose proteins were not detected with EggNOG: ")
    print(missing_plasmids)
    return df1

def Physical():
    df = pd.read_excel(physical, index_col = 0, header = 0)
    df = df.dropna()
    df['St_Depth'] = df['Station'].astype(int).astype(str)+ '_'+df['Depth'].astype(str)
    df["Temperature"] = df['Temp.'].astype(float).round()
    return df[['St_Depth','Depth','Temperature']]

def Frequency_ofCategory():
    df_station = MapToPlace()[['St_Depth', 'COG cat', 'Functional categories']]
    df_categories = df_station[['COG cat', 'Functional categories']].drop_duplicates()
    df_plasmids = SplitColumn()[['Query', 'COG cat']]
    df_plasmids['Plasmid'] = df_plasmids['Query'].apply(lambda x: (re.search('^.*\_', x).group(0)[:-1]))
    df_plasmids['COG cat'] = df_plasmids['COG cat'].replace('-', 'missing').fillna('missing')
    df_plasmids = df_plasmids[(df_plasmids['COG cat'] != 'missing') & (df_plasmids['COG cat'] != 'S')]
    df_grouped = df_plasmids['COG cat'].value_counts(normalize = True).to_frame(name='Function count').sort_values(by='Function count', ascending = False)
    df_grouped['Function count']=df_grouped['Function count']*100
    result = pd.concat([df_grouped, df_categories.set_index('COG cat')['Functional categories']], axis = 1)
    #df_grouped['Functional categories'] = df_grouped['COG cat'].map(df_categories.set_index('COG cat')['Functional categories'])
    print("******************* COG categories frequency by proteins *****************")
    print(result.head())
    df_station['COG cat'] = df_station['COG cat'].replace('-', 'missing').fillna('missing')
    df_station = df_station[(df_station['COG cat'] != 'missing') & (df_station['COG cat'] != 'S')]
    df_grouped = df_station['COG cat'].value_counts(normalize = True).to_frame(name = 'Function count').sort_values(
        by = 'Function count', ascending = False)
    df_grouped['Function count'] = df_grouped['Function count'] * 100
    result = pd.concat([df_grouped, df_categories.set_index('COG cat')['Functional categories']], axis = 1)
    print("******************* COG categories frequency by stations *****************")
    print(result.head())

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

#Frequency_ofCategory()
#eggNOGStats()
BarChart(True, 'with')
BarChart(True, 'without')
#BarChart_lim(Plasmid_class()[1], cluster_pl_dfPlPut, 'barplot_COG_PlPut', 'with')
#BarChart_lim(Plasmid_class()[1],cluster_pl_dfPlPut, 'barplot_COG_PlPut', 'without')
#BarChart_lim(Plasmid_class()[0], cluster_pl_df7, 'barplot_COG_7pl', 'with')
#BarChart_lim(Plasmid_class()[0], cluster_pl_df7, 'barplot_COG_7pl', 'without')