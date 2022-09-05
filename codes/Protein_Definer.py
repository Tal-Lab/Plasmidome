# -*- coding: utf-8 -*-
"""
Created on 27/10/2021 9:33

Author: Lucy
"""
### Description
# description!!!!

version=12

import pandas as pd
import numpy as np
import re
import os
import sys
from Bio import SeqIO
from pathlib import Path
import timeit
from scipy.cluster.hierarchy import linkage, fcluster
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
from Protein_files_handler import dataset
from station_phys_clustering import Clust_map2, Plasmid_class, Physical, plasmids_by_reads, Station, GetLibSize
from plasmid_detect import Function_ORF

station_orderCl, station_reorderCl, cluster_st_df, cluster_pl_df = Clust_map2(4,Plasmid_class()[2],'PlPutUnc_HMannot_', 1150, 1500)
station_orderPlPut, station_reorderPlPut, cluster_st_dfPlPut, cluster_pl_dfPlPut = Clust_map2(4,Plasmid_class()[1],'PlPut_HMannot_', 800, 900)
station_order7, station_reorder7,cluster_st_df7, cluster_pl_df7 = Clust_map2(4,Plasmid_class()[0],'Pl_HMannot_', 250, 400)

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# uncomment relevant path to OS
# Windows
path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
#path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

# working directories
output_dir = f"{path}/data_calculations"
visuals = f"{path}/visualisations"
Path(output_dir).mkdir(parents=True, exist_ok=True)

# working files
#### blast results from different databases
BacMetExp = f'{dataset}/BacMetExp.csv'
BacMetPred = f'{dataset}/BacMetPred.csv'
AB_res_1 = f'{dataset}/protein_fasta_protein_homolog_model.csv'
AB_res_2 = f'{dataset}/protein_fasta_protein_knockout_model.csv'
AB_res_3 = f'{dataset}/protein_fasta_protein_overexpression_model.csv'
AB_res_4 = f'{dataset}/protein_fasta_protein_variant_model.csv'
vir_fact = f'{dataset}/VFDB_setB_pro.csv'
#### annotations
BMExp_annot = r"../res/Annotations/BacMetExp_annot.txt"
BMPred_annot = r"../res/Annotations/BacMetPred_annot.txt"
CARD_annot = r"../res/Annotations/aro_index.tsv"
AB_classes = r"../res/Annotations/antibiotics.xlsx"
VF_annot = r"../res/Annotations/VFs.xls"
#### stations and reads coverage
proteins_fasta = r"../res/Filtered_ORFs.fasta"
reads_coverage = r"../res/all_cov.csv"

plasmids_byreads=f'{output_dir}/Plasmids_ByReads.csv'

colnames = ['qseqid', 'sseqid', 'stitle', 'evalue', 'length', 'pident', 'mismatch', 'score', 'qcovs', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq']

####### Getting predicted protein's lengths

# reading parsing and creating dataframe from protein fasta file
def ORF_fasta():
    # Create DataFrame
    df_plasmids = pd.DataFrame({'ORF':[],'ORF_length':[]})
    with open(proteins_fasta, 'r') as proteins:
        for rec in SeqIO.parse(proteins, 'fasta'):
            name = rec.id
            seqLen = len(rec)
            df=pd.DataFrame({'ORF':[name],'ORF_length':[seqLen]})
            df_plasmids=df_plasmids.append(df)
        df_plasmids["ORF_length"] = df_plasmids['ORF_length'].astype(float).round().astype(int)
        df_plasmids['Plasmid'] = df_plasmids['ORF'].apply(lambda x: (re.search('^.*\_', x).group(0)[:-1]))
        #print(df_plasmids)
    return df_plasmids

prot_DF=ORF_fasta()
####### Handling BLAST outputs

# retreiving plasmid length from node-name
def Clean_length(node_name):
    if re.search('h\_\d+', node_name):
        pos = int(re.search('h\_\d+', node_name).group(0)[2:])
        return pos
    else:
        return node_name

# retreiving protein name from protein description
def GetProt(prot_name,pattern,start, stop):
    if re.search(pattern, prot_name):
        pos = re.search(pattern, prot_name).group()[start:stop]
        return pos
    else:
        return prot_name

# getting best hits from blast files
def CleanFiltered(df):
    #cleaning data from additional hits
    grouped = df.groupby("qseqid") #grouping of dataframe by query
    result = []
    for name, group in grouped:
        todo = group #reassigning the dataframe to new variable
        while todo.shape[0] >= 1:
            a = max(todo['qstart'])
            b = max(todo['qend'])
            if a > b:
                buffer = todo.loc[todo.loc[:, 'qstart'].between(b, a), :] #defining a subgroup of overlaping proteins
            else:
                buffer = todo.loc[todo.loc[:, 'qend'].between(a, b), :] #defining a subgroup of overlaping proteins
            best_hit = buffer.sort_values('align_perc', ascending = False).groupby(['qseqid']).head(1) #defining the best hit in the subgroup
            result.append(best_hit)
            todo = pd.concat([todo, buffer]).drop_duplicates(keep = False)
    result=pd.concat(result)
    return result

####### Reading BLAST outputs

# handling Metal Resistance outputs
def MetResDF():
    df1 = pd.read_csv(BacMetExp, sep = '\t', index_col = None, header = None)
    df1.columns = colnames
    df1['Protein_ID'] = df1.apply(lambda x: GetProt(x['sseqid'], '^\w+',None,None), axis = 1)
    print(df1)
    #df1=df1.sort_values('length', ascending = False).groupby(['qseqid']).head(1)
    df2 = pd.read_csv(BacMetPred, sep = '\t', index_col = None, header = None)
    df2.columns = colnames
    df2['Protein_ID'] = df2.apply(lambda x: GetProt(x['sseqid'], '\|\d+',1, None), axis = 1)
    #df2 = df2.sort_values('length', ascending = False).groupby(['qseqid']).head(1)
    metDF = pd.concat([df1,df2])
    metDF['query_length'] = metDF['qseqid'].apply(Clean_length)
    hits = 'MetalRes'
    metDF.loc[:, hits] = metDF.groupby(['qseqid'])['qseqid'].transform('count')
    metDF = pd.merge(metDF, prot_DF, left_on = 'qseqid', right_on = "ORF")
    del metDF['ORF']
    metDF['align_perc'] = metDF['length'] / metDF['ORF_length']
    metDF = CleanFiltered(metDF)
    print('############ Print first 10 metDF ###########')
    print(metDF.head())
    return metDF[['qseqid', hits, 'Protein_ID', 'length', 'pident', 'qstart','qend','sseq', 'ORF_length','Plasmid']]

# handling AMR outputs
def ABResDF():
    dataset=[AB_res_1,
             AB_res_2,
             AB_res_3,
             AB_res_4]
    df = pd.DataFrame()
    for file in dataset:
        df_ab = pd.read_csv(file, sep = '\t', index_col = None, header = None)
        df_ab.columns = colnames
        df_ab['Protein_ID'] = df_ab.apply(lambda x: GetProt(x['sseqid'], 'ARO:\d+',None,None), axis = 1)
        #print(df_ab)
        df=df.append(df_ab)
    df['query_length'] = df['qseqid'].apply(Clean_length)
    hits = 'AB_Res'
    df.loc[:, hits] = df.groupby(['qseqid'])['qseqid'].transform('count')
    df = pd.merge(df, prot_DF, left_on = 'qseqid', right_on = "ORF")
    del df['ORF']
    df['align_perc'] = df['length'] / df['ORF_length']
    df= CleanFiltered(df)
    print('############ Print first 10 CARD ###########')
    print(df.head())
    return df[['qseqid', hits, 'Protein_ID', 'length', 'pident', 'qstart', 'qend','sseq', 'ORF_length','Plasmid']]

# handling Virulence factors output
def VF_DF():
    df = pd.read_csv(vir_fact, sep = '\t', index_col = None, header = None)
    df.columns = colnames
    df['VF_name'] = df.apply(lambda x: GetProt(x['stitle'], '(\[\w+.*\()', 1,-2),axis = 1)
    df['Bacteria'] = df.apply(lambda x: GetProt(x['stitle'], '(\[\w+ \w+)', 1,None),axis = 1)
    print(df['VF_name'].unique())
    df['query_length'] = df['qseqid'].apply(Clean_length)
    hits = 'Virulence_Factors'
    df.loc[:, hits] = df.groupby(['qseqid'])['qseqid'].transform('count')

####### Reading annotation files

# reading MR annotation file
def Annotation(file):
    df = pd.read_csv(file, delimiter = "\t", header = 0, index_col = None)
    df = df.rename(columns={df.columns[0]: "Protein_ID" })
    s = df['Compound'].str.split(', ').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Compound'
    del df['Compound']
    df = df.join(s)
    #print(df)
    return df

# reading AMR annotation file
def AnnotCARD():
    df = pd.read_csv(CARD_annot, sep = '\t', header = 0, index_col = None)
    df = df.rename(columns = {df.columns[0]: "Protein_ID", df.columns[9]: "Compound"})
    s = df['Compound'].str.split(';').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Compound'
    del df['Compound']
    df = df.join(s)
    df["Compound"]=df['Compound'].str.replace(' antibiotic','')
    df = df[df['Compound'].notna()]
    df_clean=df[["Protein_ID","Compound"]]
    #print(df_clean)
    df2=pd.read_excel(AB_classes, index_col=None, header = 0)
    #print(df2)
    df_DC= pd.merge(df_clean, df2,left_on = 'Compound', right_on = 'Common')
    return df_DC

####### Handling stations

# defining location for each plasmid
def DefineStation2():
    df = pd.read_csv(plasmids_byreads, index_col=None, header = 0)
    df.dropna(inplace = True)
    #print('######## Printing plasmids by reads df ############')
    #print(df)
    return df

# getting station parameters for each plasmids
def GroupPlaceMapper():
    df1 = DefineStation2()
    df2 = Station()
    df1['St_Depth'] = df1['station_name'].map(df2.set_index('station_name')['St_Depth'])
    df1 = df1.rename(columns = {'NewName': "qseqid"})
    #df1['station_name'] = df1['Station'] + "_coverage"
    df1['Depth'] = df1['St_Depth'].apply(lambda x: (int(re.search("_\d+", x).group(0)[1:])))
    #df1.loc[:, 'hits'] = df1.groupby(['St_Depth'])['qseqid'].transform('count')
    return df1

def FunctionMapper(name):
    if name == "BACMET":
        bmPred = Annotation(BMPred_annot)[['Protein_ID','Compound']]
        bmExp = Annotation(BMExp_annot)[['Protein_ID','Compound']]
        fun_anot = bmPred.append(bmExp)
        df_plasm = MetResDF()
    elif name == "CARD":
        fun_anot = AnnotCARD()[['Protein_ID','Drug_Class']]
        fun_anot = fun_anot.rename(columns = {"Drug_Class": "Compound"})
        print(fun_anot.head())
        df_plasm = ABResDF()
        print(df_plasm.head())
    new_df = pd.merge(df_plasm, fun_anot, on = 'Protein_ID', how = 'inner')
    #print("########### Printing df for %s" % name)
    #print(new_df[['qseqid', 'Compound']])
    return new_df

def Clean_Compound(compound_name):
    if re.search('\w+\)]', compound_name):
        pos = re.search('\w+\)]', compound_name).group(0)[:-2]
        return pos
    else:
        if re.search(':.*', compound_name):
            pos = re.search(':.*', compound_name).group(0)[2:-1]
            return pos
        elif re.search('\w+\)', compound_name):
            pos = re.search('\w+\)', compound_name).group(0)[:-1]
            return pos
        return compound_name

def Metal_Compound():
    df = FunctionMapper("BACMET")
    df = df.rename(columns = {df.columns[10]: "Compound_long"})
    df['Compound'] = df['Compound_long'].apply(Clean_Compound)
    return df

def GroupCARDANnot():
    df1 = FunctionMapper("CARD")
    list_to_compare = df1['Compound'].unique()
    df2 = AnnotCARD()
    df2 = df2[df2['Compound'].isin(list_to_compare)]
    card_name = f'{output_dir}/{"card_annot_byplasmids.csv"}'
    if not os.path.isfile(card_name) or os.stat(card_name).st_size == 0:
        df2.to_csv(card_name, index = False)
    return card_name

def CARDfuncs():
    #card_file = GroupCARDANnot()
    card_file = f'{output_dir}/{"card_annot_byplasmids.csv"}'
    card_df = pd.read_csv(card_file)
    print(card_df)
    card_df = card_df.rename(columns = {card_df.columns[1]: "Protein_name",
                                        card_df.columns[4]: "Compound"})
    #print(card_df)
    s = card_df['Compound'].str.split(', ').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Compound'
    del card_df['Compound']
    card_df = card_df.join(s)
    #print(card_df['Compound'].nunique())
    return card_df

def CARD_to_Plasm():
    df_plasm = ABResDF()
    df_card = CARDfuncs()[['Protein_ID','Compound']]
    new_df = pd.merge(df_plasm, df_card, how = 'inner', on = 'Protein_ID')
    print(new_df)
    return new_df

def MapToPlace(name):
    if name == "BACMET":
        df_plasm = Metal_Compound()
    elif name == "CARD":
        df_plasm = FunctionMapper("CARD")
    df_plasm.drop('qseqid', axis=1,inplace=True)
    df_plasm['Plasmid'] = df_plasm['Plasmid'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    df_stations = GroupPlaceMapper()[['qseqid', 'St_Depth']]
    #new_df = df_stations.merge(df_plasm, left_on = 'qseqid', right_on = 'Plasmid', how='left')
    if name == "BACMET":
        new_df = pd.merge(df_plasm, df_stations, left_on = 'Plasmid', right_on = 'qseqid', how = 'right').reindex(
            columns = ['qseqid', 'St_Depth', 'Compound', 'MetalRes'])
        new_df.fillna('missing', inplace = True)
        new_df['MetalRes'] = new_df['MetalRes'].apply(lambda x: round(x) if x!= 'missing' else x)
    elif name == 'CARD':
        new_df = pd.merge(df_plasm, df_stations, left_on = 'Plasmid', right_on = 'qseqid', how = 'right').reindex(
            columns = ['qseqid', 'St_Depth', 'Compound', 'AB_Res'])
        new_df.fillna('missing', inplace = True)
        new_df['AB_Res'] = new_df['AB_Res'].apply(lambda x: round(x) if x!= 'missing' else x)
    new_df = new_df.rename({'qseqid': 'Plasmid'}, axis = 1)
    new_df = new_df.loc[new_df['Plasmid'] != '94_LNODE_1']
    return new_df

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

def MergeFunct(name):
    func = MapToPlace(name)
    df1 = GroupPlaceMapper()
    new_df = pd.merge(func, df1, left_on=['Plasmid', 'St_Depth'], right_on = ['qseqid','St_Depth'])
    #print(new_df.columns)
    df1_t = df1[['station_name', 'St_Depth']]
    new_df_t = new_df[['station_name', 'St_Depth']]
    df3 = new_df_t.merge(df1_t, how = 'outer' ,indicator=True).loc[lambda x : x['_merge']=='right_only'].drop_duplicates().drop(['_merge'], axis =1)
    print("This is ammount of missing stations with %s: %s" % (name, str(len(df3))))
    prot_list = ['None'] * int(len(df3))
    df3['Compound'] = prot_list
    df3['Function frequency'] = np.zeros(int(len(df3)))
    #print(df3)
    return df3

def FreqFuncStat(name):
    df = MapToPlace(name)[['St_Depth', 'Compound']]
    print('######### Printing mapped resistance genes to sampling points #########')
    print(df.head())
    #df_to_append = MergeFunct(name)
    df2 = GroupPlaceMapper()[['station_name', 'St_Depth']].drop_duplicates()
    #print('############ Printing GroupPlaceMapper ############')
    #print(df2.head())
    df_grouped = df.groupby('St_Depth')['Compound'].value_counts(normalize = False).to_frame(name='Function frequency')
    #print("########### Printing df_grouped from FreqFuncStat ##############")
    #print(df_grouped.head())
    df_grouped = df_grouped.reset_index()
    df_grouped['station_name'] = df_grouped['St_Depth'].map(df2.set_index('St_Depth')['station_name'])
    #df_grouped = df_grouped.append(df_to_append)
    df_grouped = df_grouped.rename(columns = {df.columns[0]: 'Sampling Stations'})
    pivot_df = pd.pivot_table(df_grouped,
                              index='Sampling Stations',
                              columns='Compound',
                              values='Function frequency',
                              dropna = False,
                              fill_value = 0)
    pivot_df.drop('missing', axis=1, inplace = True)
    #print('################### pivoted')
    #print(pivot_df)
    return pivot_df

def Station_Order(order, init_order):
    counter_list = list(enumerate(init_order, 0))
    #print("prinying counter list")
    #print(counter_list)
    station_order = order
    #print("Printing station order")
    #print(station_order)
    init_df = pd.DataFrame(counter_list, columns = ['Number', 'Sampling Station']).set_index('Number').reindex(station_order)
    #print(init_df)
    final_order = init_df['Sampling Station'].to_list()
    #print(final_order)
    return final_order

def Clustermap(name, set_p):
    '''
    This function creates clustermap for AR and MR based on name
    for set of plasmid candidates (plasmids/putative/all based on set_p)
    P.S.to preserve initial clustering order; don't forget uncomment <row_cluster=False> in clustermap
    '''
    ### getting correct set of plasmid candidates
    if set_p == 'Pl':
        station_order = Station_Order(station_reorder7, station_order7)
        candidates = Plasmid_class()[0]['Plasmid'].unique()
        # get row_clusters
        stat_cluster = dict(zip(cluster_st_df7['Sampling points clusters'].unique(),
                                sns.color_palette("colorblind", cluster_st_df7['Sampling points clusters'].nunique())))
        cluster_st = cluster_st_df7['Sampling points clusters'].map(stat_cluster)
    elif set_p == 'PlPut':
        station_order = Station_Order(station_reorderPlPut, station_orderPlPut)
        candidates = Plasmid_class()[1]['Plasmid'].unique()
        # get row_clusters
        stat_cluster = dict(zip(cluster_st_dfPlPut['Sampling points clusters'].unique(),
                                sns.color_palette("colorblind", cluster_st_dfPlPut['Sampling points clusters'].nunique())))
        cluster_st = cluster_st_dfPlPut['Sampling points clusters'].map(stat_cluster)
    else:
        station_order = Station_Order(station_reorderCl, station_orderCl)
        candidates = Plasmid_class()[1]['Plasmid'].unique()
        # get row_clusters
        stat_cluster = dict(zip(cluster_st_df['Sampling points clusters'].unique(),
                                sns.color_palette("colorblind", cluster_st_df['Sampling points clusters'].nunique())))
        cluster_st = cluster_st_df['Sampling points clusters'].map(stat_cluster)
    # df_init = FreqFuncStat(name)
    #print(candidates)
    df_init = MapToPlace(name)[['Plasmid', 'St_Depth', 'Compound']]
    df = df_init.loc[df_init['Plasmid'].isin(candidates)]
    df = df[['St_Depth', 'Compound']]
    #print('######### Printing mapped resistance genes to sampling points #########')
    #print(df.head())
    # df_to_append = MergeFunct(name)
    df2 = GroupPlaceMapper()[['station_name', 'St_Depth']].drop_duplicates()
    #print('############ Printing GroupPlaceMapper ############')
    #print(df2.head())
    df_grouped = df.groupby('St_Depth')['Compound'].value_counts(normalize = False).to_frame(
        name = 'Function frequency')
    #print("########### Printing df_grouped from FreqFuncStat ##############")
    #print(df_grouped.head())
    df_grouped = df_grouped.reset_index()
    df_grouped['station_name'] = df_grouped['St_Depth'].map(df2.set_index('St_Depth')['station_name'])
    # df_grouped = df_grouped.append(df_to_append)
    df_grouped = df_grouped.rename(columns = {df.columns[0]: 'Sampling points'})
    pivot_df = pd.pivot_table(df_grouped,
                              index = 'Sampling points',
                              columns = 'Compound',
                              values = 'Function frequency',
                              dropna = False,
                              fill_value = 0)
    pivot_df.drop('missing', axis = 1, inplace = True)
    pivot_df = pivot_df.reindex(station_order)
    df_norm = pivot_df
    df_norm[:] = np.where(df_norm == 0, 0, 1)

    # get dataframe with library sizes
    library = GetLibSize()[['St_Depth', 'Size']].set_index('St_Depth')
    library = library.rename(columns = {library.columns[0]: 'Library Size'})
    # Sort the values for the bar plot to have the same order as clusters
    library = library.reindex(station_order)
    #print("printing library")
    #print(library)
    sns.set(font_scale = 1.2)
    figure1 = sns.clustermap(data = df_norm,
                             metric = "euclidean",
                             method = 'ward',
                             row_cluster = False,
                             row_colors = cluster_st,
                             linewidths = 0.0,
                             figsize = (14, 14),
                             cmap = sns.color_palette("Blues", as_cmap = True),
                             xticklabels = True,
                             yticklabels = True,
                             rasterized = True,
                             cbar_pos = None
                             )
    figure1.ax_col_dendrogram.remove()
    figure1.ax_row_dendrogram.remove()
    # get heatmap position
    hm = figure1.ax_heatmap.get_position()
    plt.setp(figure1.ax_heatmap.yaxis.get_majorticklabels(), fontsize = 'small')
    plt.setp(figure1.ax_heatmap.xaxis.get_majorticklabels(), fontsize = 'small')
    figure1.ax_heatmap.set_position([hm.x0, hm.y0, hm.width, hm.height])
    figure1.gs.update(right = 0.80)

    # divide existing axes
    divider = make_axes_locatable(figure1.ax_heatmap)

    # create new axes for bar plot
    ax = divider.append_axes("right", size = "20%", pad = 1.1)

    # set parameters for bar plot
    ax.set_facecolor('white')
    ax.set_yticklabels([])
    ax.set_xlabel('Library size, Mbp')
    t = ['0', str(round(((library['Library Size'].max()) / (10 ** 6)), 2))]
    ax.set_xticklabels(t)
    ax.set_ylim(-0.5, len(library.index) - .5)
    ax.invert_yaxis()
    # plot bar plot in ax
    ax.barh(np.arange(len(library.index)), library['Library Size'].values)

    svg_name = name + '_'+ set_p+ str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = name + '_'+ set_p+ str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    if not os.path.isfile(svg_dir) and not os.path.isfile(png_dir):
        plt.savefig(svg_dir, format='svg', dpi=gcf().dpi, bbox_inches='tight')
        plt.savefig(png_dir, format='png', dpi=gcf().dpi, bbox_inches='tight')
    plt.show()

def Function_Frequency(name):
    df = FreqFuncStat(name)
    df_norm = df
    df_norm[:] = np.where(df_norm == 0, 0, 1)
    df.loc['Total']= df.sum()
    number_of_rows = len(df.index)
    df.loc['Frequency'] = (df.loc['Total']/number_of_rows)*100
    df_small = df.loc['Frequency'].sort_values(ascending=False)
    print("******************* %s compounds frequency by stations *****************" % name)
    print(df_small.nlargest())

def ORF_table():
    CARD = FunctionMapper("CARD")[['qseqid', 'Compound']]
    CARD = CARD.rename(columns = {CARD.columns[1]: "Compound-CARD"})
    #print(CARD.head())
    BacMet = FunctionMapper("BACMET")[['qseqid', 'Compound']]
    BacMet = BacMet.rename(columns = {BacMet.columns[1]: "Compound-BacMet"})
    #print(BacMet.head())
    other_funcs = Function_ORF()
    other_funcs2 = other_funcs.merge(BacMet, left_on="ORF_name", right_on = 'qseqid',  how = 'inner')
    del other_funcs2['qseqid']
    print(other_funcs2.head())
    all_functions = other_funcs2.merge(CARD, left_on = "ORF_name", right_on = 'qseqid', how = 'inner')
    del all_functions['qseqid']
    print(all_functions.head())
    orfs_name = f'{output_dir}/{"orf_functions_overall.csv"}'
    if not os.path.isfile(orfs_name) or os.stat(orfs_name).st_size == 0:
        all_functions.to_csv(orfs_name, index = False)

#ORF_table()
#Function_Frequency("BACMET")
#Function_Frequency("CARD")
#Clustermap('BACMET', 'All')
#Clustermap('CARD', 'All')
