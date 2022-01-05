import pandas as pd
import re
from Bio import SeqIO
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from sklearn import preprocessing
from scipy.cluster.hierarchy import linkage, fcluster
import dash_core_components as dcc
import dash_bio as dashbio
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
import scipy.cluster.hierarchy as sch
from collections import defaultdict
import matplotlib.patches as mpatches

### Description
# add description

# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

# working directories
out_dir = f"{path}/data_calculations"

# working files
combined_output = r"../res/unique_plasmids.fasta"
proteins = r"../res/Filtered_ORFs.fasta"
cog_categories = r"../res/cog_cats.csv"
eggnog = r"../res/eggnog_FilteredORFs.csv"
groups = r"../res/old_to_new_plasmid.csv"
stations = r"../res/stations.txt"
reads_coverage = r"../res/all_cov.csv"
plasmids_byreads = f"{out_dir}/AllPlasmids_perStat.csv"
physical = r"../res/station_physical.xlsx"

def csv_reader(file):
    df = pd.read_csv(file, sep = ',', header = 0, index_col = None)
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

def Clean_Stat(node_name):
    if re.search("^\d+", node_name):
        pos = str(re.search("^\d+", node_name).group(0)[:])
        return pos
    else:
        return node_name

def DefineStation():
    df = pd.read_csv(groups, index_col=None, header = 0)
    df.dropna(inplace = True)
    df['Station'] = df['Plasmid'].apply(Clean_Stat)
    return df

def DefineStation2():
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
    df1 = DefineStation2()
    df2 = Station()
    df1['St_Depth'] = df1['station_name'].map(df2.set_index('station_name')['St_Depth'])
    #df1['station_name'] = df1['Station'] + "_coverage"
    df1['Depth'] = df1['St_Depth'].apply(lambda x: (int(re.search("_\d+", x).group(0)[1:])))
    return df1

def MapToFunc():
    df_plasm = SplitColumn()[['Query','COG cat']]
    df_func = csv_reader(cog_categories)
    df_plasm['Functional categories'] = df_plasm['COG cat'].map(df_func.set_index('COG cat')['Functional categories'])
    return df_plasm

def MapToPlace():
    df_plasm = MapToFunc()
    df_plasm['qseqid'] = df_plasm['Query'].apply(lambda x: re.search('^.*\_', x).group(0)[:-1])
    df_stations = GroupPlaceMapper()[['NewName', 'St_Depth']]
    new_df = pd.merge(df_plasm, df_stations, left_on = 'qseqid', right_on = 'NewName')
    #print(new_df)
    return new_df

def FuncCountTotal():
    df = MapToPlace()[['St_Depth', 'COG cat', 'Functional categories']]
    df2 = df[['COG cat', 'Functional categories']].drop_duplicates()
    stat_num = df['St_Depth'].nunique()
    #print(df2)
    counts = df.groupby('COG cat')['St_Depth'].nunique().to_frame(name = 'Function Occurance')
    counts.loc[:,'Function Frequency, %'] = (counts['Function Occurance']/stat_num)*100
    counts = counts.reset_index()
    counts['Functional categories'] = counts['COG cat'].map(df2.set_index('COG cat')['Functional categories'])
    counts['Functional categories'] = counts[['COG cat', 'Functional categories']].apply(tuple, axis = 1)
    counts['Functional categories'] = counts['Functional categories'].apply(lambda x: ' - '.join(x))
    print(counts)
    return counts

def FreqFuncStat(freq):
    df = MapToPlace()[['St_Depth', 'COG cat', 'Functional categories']]
    df2 = df[['COG cat', 'Functional categories']].drop_duplicates()
    df3 = GroupPlaceMapper()[['station_name', 'St_Depth']].drop_duplicates()
    indexNames = df[df['COG cat'] == 'S'].index
    df.drop(indexNames, inplace = True)
    df_grouped = df.groupby('St_Depth')['COG cat'].value_counts(normalize = freq).to_frame(name='Function count')
    df_grouped = df_grouped.reset_index()
    df_grouped['Functional categories'] = df_grouped['COG cat'].map(df2.set_index('COG cat')['Functional categories'])
    df_grouped['Functional categories'] = df_grouped[['COG cat', 'Functional categories']].apply(tuple, axis = 1)
    df_grouped['Functional categories'] = df_grouped['Functional categories'].apply(lambda x: ' - '.join(x))
    df_grouped['station_name'] = df_grouped['St_Depth'].map(df3.set_index('St_Depth')['station_name'])
    print(df_grouped)
    return df_grouped

def CoverageDF():
    data = pd.read_csv(reads_coverage, sep = ',', index_col = 0, header = 0)
    columns = list(data.columns)
    arr1 = data[columns].to_numpy()
    arr = arr1.transpose()
    return arr, data

def ClusterMapper ():
    array, df = CoverageDF()
    Z = linkage(array, 'complete', optimal_ordering=True)
    clusters = fcluster(Z, t = 1.2, criterion = 'distance')
    clusters = list(clusters)
    labels = list(df.columns)
    new_label_num = list(zip(clusters, labels))
    new_label_num.sort()
    new_label_num = [x[1] for x in new_label_num]
    #print(new_label_num)
    return new_label_num

def BarChart(freq, y_name):
    df = FreqFuncStat(freq)
    order_Stat = ClusterMapper()
    true_sort1 = [s for s in order_Stat if s in df.station_name.unique()]
    # print(true_sort)
    df = df.set_index('station_name').loc[true_sort1].reset_index()
    station = df['St_Depth'].to_list()
    fig = px.bar(df, x = station, y = 'Function count',
                 color = 'Functional categories',
                 labels = {'Function count' : y_name, 'x' : 'Station'})
    fig.update_layout(font = dict(size = 18))
    fig.show()
    fig.write_html(y_name + '4.html')
    fig.write_image(y_name + '4.svg')

def PieChart():
    df1 = FuncCountTotal()
    fig1 = px.bar(df1, x = 'COG cat', y = 'Function Frequency, %',
                 color = 'Functional categories')
    fig1.update_layout(font = dict(size = 18))
    fig1.show()
    fig1.write_html("AllFunctions_PieChart.html")
    fig1.write_image("AllFunctions_PieChart.svg")

def PlasmidList():
    with open(combined_output) as all_plasdmids:
            non_unique_list = [] #list for plasmids
            for record in SeqIO.parse(all_plasdmids, "fasta"):
                non_unique_list.append(record.id)
    return non_unique_list

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
    fig = px.histogram(df_grouped, x = "Number of Proteins", histnorm = 'percent', nbins = 50)
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
    df1['COG cat'] = df1['Proteins'].map(df2.set_index('Query')['COG cat']).replace('', 'missing').fillna('missing')
    #df1['COG cat'] = df1['COG cat']
    print(df1['COG cat'].unique())
    eggnog_perc = eggnog_n/prodigal*100
    formatted_eggnog = "{:.2f}".format(eggnog_perc)
    print("eggNOG mapped only " + str(formatted_eggnog) + "% of predicted genes")
    ### this part is only for detecting proteins with assigned COG categories
    df2['COG cat'] = df2['COG cat'].replace('', 'missing').fillna('missing')
    df2 = df2[df2['COG cat'] != 'missing']
    #print(df2)
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

def DetectNANs():
    df=eggNOGStats()
    df_grouped = df.groupby(['Plasmids','COG cat']).size()
    df2 = df_grouped.unstack()
    df_cols = df2.columns.to_list()[:-1]
    print(df_grouped)
    df_rows = df2.index.to_list()
    qty_of_nuls=len(df_cols)
    #print(qty_of_nuls)
    df3 = df2[df2[df_cols].isna().sum(axis = 1) > qty_of_nuls]
    #df3 = df2[df2[(df2.isnull().sum(axis = 1) > qty_of_nuls)]]
    #print(df3)

def ClusterMap(name):
    df_plasm = SplitColumn()[['Query','COG cat']]
    df_plasm['qseqid'] = df_plasm['Query'].apply(lambda x: re.search('^.*\_', x).group(0)[:-1])
    df_plasm = df_plasm.drop(['Query'], axis=1)
    df_plasm['Presence'] = 1
    data_piv = df_plasm.pivot_table(index= "qseqid", columns = "COG cat", values = 'Presence').fillna(0)
    columns = list(data_piv.columns.values)
    rows = list(data_piv.index)
    print(len(rows))
    arr = data_piv.loc[rows].values
    clustergram = dashbio.Clustergram(
        data = np.transpose(arr),
        row_labels = columns,
        column_labels = rows,
        center_values = False,
        optimal_leaf_order = True,
        color_map = px.colors.sequential.Plasma,
        color_threshold = {
            'row': 9.9,
            'col': 3.1
        },
        height = 500,
        width = 900,
        line_width = 2,
        hidden_labels = 'column',
        tick_font = {'size': 7})
    dcc.Graph(figure = clustergram)
    clustergram.show()
    #clustergram.write_html(name +".html")
    #clustergram.write_image(name + ".svg")

def Physical():
    df = pd.read_excel(physical, index_col = 0, header = 0)
    df = df.dropna()
    df['St_Depth'] = df['Station'].astype(int).astype(str)+ '_'+df['Depth'].astype(str)
    df["Temperature"] = df['Temp.'].astype(float).round()
    return df[['St_Depth','Depth','Temperature']]

def ClusterMap2(name):
    df_plasm = SplitColumn()[['Query', 'COG cat']]

    df_plasm['Plasmid Candidates'] = df_plasm['Query'].apply(lambda x: re.search('^.*\_', x).group(0)[:-1])
    df_plasm = df_plasm.drop(['Query'], axis = 1)
    df_plasm['Presence'] = 1
    data_piv = df_plasm.pivot_table(index = "Plasmid Candidates", columns = "COG cat", values = 'Presence').fillna(0)
    print(data_piv.index.nunique())
    data_transp = data_piv.transpose()
    figure1 = sns.clustermap(data = data_transp,
                             metric = "euclidean",
                             method = 'ward',
                             linewidths = 0.0,
                             cmap = sns.color_palette("Blues", as_cmap=True),
                             xticklabels = False,
                             yticklabels = True,
                             rasterized = True,
                             cbar_pos = None
                             )
    # g.set_axis_labels(["Plasmid candidates", "Sampling stations"])
    figure1.ax_col_dendrogram.remove()
    figure1.ax_row_dendrogram.remove()
    plt.setp(figure1.ax_heatmap.yaxis.get_majorticklabels(), fontsize = 12)
    svg_name = name + '.svg'
    png_name= name + '.png'
    plt.savefig(svg_name, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_name, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()  # Push new figure on stack

def Frequency_ofCategory():
    df_station = MapToPlace()[['St_Depth', 'COG cat', 'Functional categories']]
    df_categories = df_station[['COG cat', 'Functional categories']].drop_duplicates()
    df_plasmids = SplitColumn()[['Query', 'COG cat']]
    df_plasmids['Plasmid'] = df_plasmids['Query'].apply(lambda x: (re.search('^.*\_', x).group(0)[:-1]))
    df_plasmids['COG cat'] = df_plasmids['COG cat'].replace('', 'missing').fillna('missing')
    df_plasmids = df_plasmids[(df_plasmids['COG cat'] != 'missing') & (df_plasmids['COG cat'] != 'S')]
    df_grouped = df_plasmids['COG cat'].value_counts(normalize = True).to_frame(name='Function count').sort_values(by='Function count', ascending = False)
    df_grouped['Function count']=df_grouped['Function count']*100
    result = pd.concat([df_grouped, df_categories.set_index('COG cat')['Functional categories']], axis = 1)
    #df_grouped['Functional categories'] = df_grouped['COG cat'].map(df_categories.set_index('COG cat')['Functional categories'])
    print("******************* COG categories frequency by proteins *****************")
    print(result.head())
    df_station['COG cat'] = df_station['COG cat'].replace('', 'missing').fillna('missing')
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
#ClusterMap("Clustergram_COG")
#ClusterMap2("Clustergram_COG_3")
#eggNOGStats()
#DetectNANs()
#BarChart(True, 'Function frequency')
#BarChart(False, 'Function count')
#PieChart()
#FreqFuncStat(True)
#csv_reader(cog_categories)
#GroupPlaceMapper()#df['COG cat2'] = df['COG cat'].apply(SPLIT)
#ProteinsFastaAn()