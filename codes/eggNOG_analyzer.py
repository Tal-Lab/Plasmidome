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
from scipy.special import comb
from collections import defaultdict
import matplotlib.patches as mpatches
from station_phys_clustering import Clust_map2, plasmids_by_reads, version
from plasmid_detect import Plasmid_class
from scipy.stats import hypergeom
### Description
# add description

# uncomment relevant path to OS
path = r"../Output"

# working directories
out_dir = f"{path}/data_calculations"
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"

# working files
combined_output = r"../res/filtered_plasmids.fasta"
proteins = r"../res/Filtered_ORFs.fasta"
cog_categories = r"../res/Annotations/cog_cats.csv"
eggnog = r"../res/eggnog_FilteredORFs.csv"
stations = r"../res/stations.txt"
reads_coverage = r"../res/all_cov.csv"
physical = r"../res/station_physical.xlsx"

station_orderCl, station_reorderCl, cluster_st_df, cluster_pl_df = Clust_map2(version,'All_HMannot_', 1150, 1500)

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
    df = pd.read_csv(plasmids_by_reads, index_col=None, header = 0)
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
    sns.set(font_scale = 0.95)
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

def prob_func(x, orfs, genes, df_db, df_orfs):
    ''' Function to calculate COG categories statistics. \
    Function gets COG category (x), number of genes, \
    assigned any COG category (orfs), cog-database, orf_database.'''
    print(genes)
    print(orfs)
    # getting number of genes assigned particular COG category x from cog and orfs databases
    df_db = df_db[df_db['COG cat'] == x]
    cog_db = df_db.iloc[0]['DB representation']
    df_orfs = df_orfs[df_orfs['COG cat'] == x]
    cog_plasmids = df_orfs.iloc[0]['Function count']
    #print(genes, cog_db, orfs, cog_plasmids)
    rv = hypergeom(genes, cog_db, orfs)
    pmf_cog = rv.pmf(cog_plasmids)
    pval = hypergeom.sf(cog_plasmids-1, genes, cog_db, orfs)
    cog_freq = (cog_db/genes)*100
    cog_orf_freq = (cog_plasmids/orfs)*100
    #print (pmf_cog)
    print('The probability of getting %d ORFs assigned COG-%s out of %d ORFs is %s.' % (cog_plasmids, x, orfs, "{:.2e}".format(pmf_cog)))
    print('The probability of getting %d or more ORFs assigned COG-%s out of %d ORFs is %s.' % (cog_plasmids, x, orfs, "{:.2e}".format(pval)))
    print('The frequency of COG-%s in COG-database: (%d/%d)*100=%f' % (x, cog_db, genes, round(cog_freq,2)))
    print('The frequency of COG-%s in our plasmids: (%d/%d)*100=%f' % (x, cog_plasmids, orfs, round(cog_orf_freq, 2)))
    #not working
    #prob = hypergeom_pmf(genes, cog_db, orfs, cog_plasmids)
    #print(prob)
    return pmf_cog

def Frequency_ofCategory():
    """ Statistics for COG categories """
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
    #print(result.head())
    df_pl_count = df_plasmids.groupby('COG cat').size().reset_index(name='Function count').sort_values(by='Function count', ascending = False)
    plasmid_orfs = df_plasmids['Query'].nunique()
    #df_pl_count = df_plasmids['COG cat'].value_counts(normalize = False).to_frame(name = 'Function count').sort_values(by = 'Function count', ascending = False)
    print(df_pl_count.head())
    #print(df_pl_count[df_pl_count['COG cat']=='L'])
    # getting COG elemments in COG database
    df_cogs = csv_reader(cog_categories)
    df_cogs = df_cogs.loc[df_cogs['COG cat'] != 'S']
    # calculating all genes in COG db
    db_genes = df_cogs['DB representation'].astype('int').sum()
    # calculating probability of getting particular ammount of genes assigned a particular COG category by random
    df_pl_count['cog_prob'] = df_pl_count['COG cat'].apply(prob_func, args=(plasmid_orfs,db_genes,df_cogs, df_pl_count))
    df_pl_count['cog_prob'] = df_pl_count['cog_prob'].apply(lambda x: "{:.2e}".format(x))
    print(df_pl_count.head())
    ### getting statistics for COG elements in the sampling points
    df_station['COG cat'] = df_station['COG cat'].replace('-', 'missing').fillna('missing')
    df_station = df_station[(df_station['COG cat'] != 'missing') & (df_station['COG cat'] != 'S')]
    df_grouped = df_station['COG cat'].value_counts(normalize = True).to_frame(name = 'Function count').sort_values(
        by = 'Function count', ascending = False)
    df_grouped['Function count'] = df_grouped['Function count'] * 100
    result = pd.concat([df_grouped, df_categories.set_index('COG cat')['Functional categories']], axis = 1)
    print("******************* COG categories frequency by stations *****************")
    print(result.head())

def PieChart(df, file_name, unknown):
    ''' This function creates Piechart of COG families distribution '''
    df['COG cat'] = df['COG cat'].replace('-', 'missing').fillna('missing')
    df=df.fillna('No category assigned')
    print(df)
    if unknown == 'without':
        df = df[(df['COG cat'] != 'No category assigned') & (df['COG cat'] != 'S')]
    print(df)
    gen_num = df['Query'].nunique()
    print(gen_num)
    #df_grouped = df.groupby[['COG cat', 'Functional categories']].value_counts(normalize = True).to_frame(name = 'Function count').sort_values(by = 'Function count', ascending = False)
    df_grouped = df.groupby(['COG cat', 'Functional categories']).size().reset_index(name="Function count")
    df_grouped['Function frequency'] = (df_grouped["Function count"]/gen_num)*100
    print(df_grouped)
    #df_grouped['Function count'] = df_grouped['Function count'] * 100
    #print(df_grouped)
    df_grouped = df_grouped.sort_values('Function frequency', ascending = False)
    ncolors = df_grouped['Functional categories'].nunique()
    colors = sns.color_palette(cc.glasbey, n_colors = ncolors)
    df_grouped['Func_cat'] = df_grouped[['COG cat', 'Functional categories']].apply(tuple, axis = 1)
    df_grouped['Func_cat'] = df_grouped['Func_cat'].apply(lambda x: ' - '.join(x))
    labels = df_grouped['Func_cat'].unique()
    labels_pie = df_grouped['COG cat'].unique()
    # prepare figure
    sns.set(font_scale = 0.95)
    #plt.pie(df_grouped['Function frequency'], labels = labels_pie, colors = colors, autopct = '%0.0f%%')
    slices_pie = df_grouped['COG cat'].nunique()
    explode = [0] * slices_pie
    print(explode)
    ind = slices_pie-4
    print(ind)
    explode[ind:] = [0.2]*4
    explode = tuple(explode)
    print(explode)
    plt.pie(df_grouped['Function frequency'], colors = colors, explode = explode)
    plt.legend(labels, title = 'Functional categories (COGs)', bbox_to_anchor = (1.01, 1), ncol = 1,title_fontsize = 'medium',fontsize = 'medium', frameon = False, loc = 2, borderaxespad = 0.)
    # fig.tick_params(axis = 'x', rotation = 90, labelsize = 8)
    # save graph in PNG and vector format
    svg_name = file_name + "_" + unknown + str(3) + '.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name =  file_name + "_"  + unknown + str(3) + '.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        plt.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        plt.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()

def table():
    df_cog = MapToFunc()
    df_grouped = df_cog.groupby(['COG cat','Functional categories']).count().sort_values(by = 'Query', ascending = False)
    all_funcs = df_grouped['Query'].sum()
    print(all_funcs)
    no_unknown = all_funcs - df_grouped.loc[("S","Function Unknown"),"Query"]
    print(no_unknown)
    df_grouped['Function frequency'] = df_grouped['Query'].apply(lambda x: round(((x/all_funcs)*100),2))
    df_grouped['Function frequency without unknown function'] = df_grouped['Query'].apply(lambda x: round(((x/no_unknown)*100),2))
    df_grouped.loc[("S","Function Unknown"), 'Function frequency without unknown function'] = '-'
    #df_grouped=df_grouped.round({'Function frequency': 3, 'Function frequency without unknown function': 3})
    print(df_grouped)
    df_grouped.to_excel(f'{tables}/COG_table.xlsx')



#table()
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
#Frequency_ofCategory()
#eggNOGStats()
#BarChart(True, 'with')
#BarChart(True, 'without')
#PieChart(MapToFunc(),'PieChart_COG','with')
#PieChart(MapToFunc(),'PieChart_COG','without')