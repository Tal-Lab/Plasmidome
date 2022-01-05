# -*- coding: utf-8 -*-
"""
Created on 29/08/2021 11:05

Author: Lucy Androsiuk
"""
### Description
# add description

version=1

import numpy as np
import pandas as pd
import os, re
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import seaborn as sns
from scipy import stats
from plasmid_detect import Plasmid_class

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

# working directories
tables = f"{path}/data_calculations"
visuals = f"{path}/visualisations"
Path(tables).mkdir(parents=True, exist_ok=True)

# working files
reads_coverage = r"../res/all_cov.csv"
proteins = r"../res/Filtered_ORFs.fasta"
library = r"../res/LibrarySize.csv"
stations = r"../res/stations.txt"

plasmids = Plasmid_class()[0]['Plasmid'].unique().tolist()
putative_plasmids = Plasmid_class()[1]['Plasmid'].unique().tolist()
all_candidates = Plasmid_class()[2]['Plasmid'].unique().tolist()

def CoverageDF(work_set):
    ' this function reads coverage file as a dataframe '
    data = pd.read_csv(reads_coverage, sep = ',', index_col = None, header = 0)
    #data = data.reset_index()
    data['rname'] = data['rname'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    data = data.set_index('rname')
    #print([i for i in data.index])
    #print(work_set)
    desired_indices = [i for i in data.index if i in work_set]
    desired_df = data.loc[desired_indices]
    rows = list(desired_df.index)
    arr = desired_df.loc[rows].values
    return (arr, data)
#CoverageDF(plasmids)

def PlasmidsbyReads(work_set):
    ' this function generates plasmid presence matrix where 0 is plasmid absence (plasmids coverage < 99)    \
    and 1 is plasmid presence in the sampling station (plasmids coverage >= 99) '
    arr, df = CoverageDF(work_set)
    df_norm = df
    df_norm[:] = np.where(df_norm < 99, 0, 1)
    #print(df_norm)
    #df_norm = df_norm.append(df_norm.agg(['sum']))
    return df_norm

def DF_plasmids_byReads(work_set, file_name):
    ' this function generates dataframe with plasmids present in stations and writes it to file '
    df = PlasmidsbyReads(work_set)
    df = df.stack().to_frame()
    df = df.reset_index(level=[0,1])
    df.columns = ['NewName', 'station_name', 'value']
    df=df[df['value'] > 0]
    #print(df[['NewName', 'station_name']])
    out_file = f'{tables}/{file_name}'
    #print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        df[['NewName', 'station_name']].to_csv(out_file, index = False)
    return df[['NewName', 'station_name']], out_file

DF_plasmids_byReads(plasmids, 'Plasmids_ByReads.csv')

def Coverage_stat(work_set):
    ' plasmid presence statistics '
    arr, df = CoverageDF(work_set)
    df_norm = df
    df_norm[:] = np.where(df_norm < 99, 0, 1)
    df.loc[:, 'elements'] = df.sum(axis = 1)
    df.sort_values(by = 'elements', ascending = False, inplace = True)
    df2 = df[df['elements'] == 1]
    print("****************** These are plasmids observed at 1 station only *********************")
    #print(df.index.unique())
    print(df2.index.nunique())

#statistics=Coverage_stat()

def ORF_stats():
    ' ORF statistics '
    records = []
    lenghts=[]
    for record in SeqIO.parse(proteins, "fasta"):
        records.append(record.id)
        prot_len=len(record)
        lenghts.append(prot_len)
    plasmids_w_proteins = [(re.search('^.*\_', x).group(0)[:-1]) for x in records]
    plasmids = list(set(plasmids_w_proteins))
    data = {'Plasmids': plasmids_w_proteins, 'Proteins': records, 'Protein Length':lenghts}
    df = pd.DataFrame.from_dict(data)
    min = df['Protein Length'].min()
    max = df['Protein Length'].max()
    mean = df['Protein Length'].mean().round(2)
    print("******************** The ORFs length in plasmids ranges %s - %s amino acids. *****************" % (str(min), str(max)))
    print("******************* The average ORFs length in plasmids is %s amino acids. *************" % str(mean))
    g=sns.histplot(df['Protein Length'])
    plt.axvline(x = df['Protein Length'].median(),
                color = 'blue',
                ls = '--',
                lw = 0.5)
    min_ylim, max_ylim = plt.ylim()
    plt.text(df['Protein Length'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(df['Protein Length'].median()))
    plt.axvline(x = df['Protein Length'].mean(),
                color = 'red',lw = 0.5)
    plt.text(df['Protein Length'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(df['Protein Length'].mean()))
    svg_name = "ProteinSizeHisto" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinSizeHisto" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    df_min=df.loc[df["Protein Length"]<500]
    g_min = sns.histplot(df_min['Protein Length'])
    plt.axvline(x = df_min['Protein Length'].median(),
                color = 'blue',
                ls = '--',
                lw = 2.5)
    min_ylim, max_ylim = plt.ylim()
    plt.text(df_min['Protein Length'].median() * 1.1, max_ylim * 0.9, 'Median: {:.2f}'.format(df_min['Protein Length'].median()))
    plt.axvline(x = df_min['Protein Length'].mean(),
                color = 'red')
    plt.text(df_min['Protein Length'].mean() * 1.1, max_ylim * 0.8, 'Mean: {:.2f}'.format(df_min['Protein Length'].mean()))
    svg_name = "ProteinSizeHistoMin" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinSizeHistoMin" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    df_grouped = df.groupby('Plasmids')['Proteins'].nunique().to_frame(name = 'Number of Proteins').reset_index()
    return df_grouped

def Clean_length(node_name):
    # getting length of each predicted plasmid (node) from its name by pattern
    if re.search('h\_\d+', node_name):
        pos = int(re.search('h\_\d+', node_name).group(0)[2:])
        return pos
    else:
        return node_name

def ORF_byPlasmid_stats():
    ' ORF per plasmid statistics '
    df_grouped=ORF_stats()
    min = df_grouped['Number of Proteins'].min()
    max = df_grouped['Number of Proteins'].max()
    mean = df_grouped['Number of Proteins'].mean().round(2)
    print("Plasmids with minimal number of ORFs:")
    result_min = df_grouped[df_grouped['Number of Proteins'] == 1]
    #print(result_min)
    print(result_min['Plasmids'].nunique())
    print("Plasmids with the highest number of ORFs:")
    result_max=df_grouped.nlargest(5,['Number of Proteins'])
    #print(result_max)
    print("******************** Number of ORFs in plasmids ranges %s - %s *****************" % (str(min), str(max)))
    print("******************* The mean number of ORFs in plasmid is %s *************" % str(mean))
    df_grouped['Plasmid Length']= df_grouped['Plasmids'].apply(Clean_length)
    #print(df_grouped)
    pearson_ORF = stats.pearsonr(df_grouped["Number of Proteins"].to_numpy(), df_grouped["Plasmid Length"].to_numpy())
    print((pearson_ORF))
    size_num = df_grouped.plot.scatter(x = "Plasmid Length", y = "Number of Proteins")
    plt.xticks(rotation = 90)
    plt.xscale('log')
    svg_name = "Len_Proteins_log" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Len_Proteins_log" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    num_hist= sns.histplot(df_grouped, x = "Number of Proteins", bins = 30)
    plt.axvline(x = df_grouped['Number of Proteins'].median(),
                color = 'blue',
                ls = '--',
                lw = 1.0)
    min_ylim, max_ylim = plt.ylim()
    plt.text(df_grouped['Number of Proteins'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(df_grouped['Number of Proteins'].median()))
    plt.axvline(x = df_grouped['Number of Proteins'].mean(),
                color = 'red',
                lw = 1.0)
    plt.text(df_grouped['Number of Proteins'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(df_grouped['Number of Proteins'].mean()))
    svg_name = "ProteinsHisto" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinsHisto" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    df_min=df_grouped.loc[df_grouped["Number of Proteins"]<100]
    num_hist = sns.histplot(df_min, x = "Number of Proteins", bins = 30)
    plt.axvline(x = df_min['Number of Proteins'].median(),
                color = 'blue',
                ls = '--',
                lw = 1.5)
    min_ylim, max_ylim = plt.ylim()
    plt.text(df_min['Number of Proteins'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(df_min['Number of Proteins'].median()))
    plt.axvline(x = df_min['Number of Proteins'].mean(),
                color = 'red',lw = 1.5)
    plt.text(df_min['Number of Proteins'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(df_min['Number of Proteins'].mean()))
    svg_name = "ProteinsHistoMin" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinsHistoMin" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    return df_grouped

def Station():
    ' getting station parameters from station matrix '
    with open(stations) as adressFile:
        matrix = np.loadtxt(adressFile, dtype = "str")
    df = pd.DataFrame(data = matrix, index = None, columns = matrix[0, 0:])
    df.reset_index()
    df.drop(index = 0, inplace = True)
    df['station_name'] = df['Sample'] + "_coverage"
    return df

# this is not the place
def annotate (x,y):
    ' calculation of Pearson correlation '
    r, p = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
            transform = ax.transAxes)

def GetLibSize():
    ' retreiving library sizes from lib-file '
    df = pd.read_csv(library, sep = ',', header = None, index_col = None)
    df.columns = ['Sample', 'Size']
    #print(df)
    df['station_name'] = df['Sample'].astype(str) + "_coverage"
    df2 = Station()
    #print(df2)
    df['St_Depth'] = df['station_name'].map(df2.set_index('station_name')['St_Depth'])
    print(df)
    size_stat = df.plot.scatter(x = "St_Depth", y = "Size")
    plt.xticks(rotation = 90)
    plt.yscale('log')
    plt.show()
    slope, intercept, r_value, p_value, std_err = stats.linregress(df["Sample" ], df["Size"])
    g = sns.regplot(x = "Sample", y = "Size", data = df, line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)})
    g = (g.set(xlim = (256, 304)))
    #g.legend()
    svg_name = "Lib_stat" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Lib_stat" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    return df

def ORF_byStation_stats():
    ' ORF per station statistics '
    df_prot=ORF_stats()
    df_stations=DF_plasmids_byReads()[0]
    df=pd.merge(df_prot, df_stations, left_on='Plasmids', right_on='NewName')
    df = df.drop(['NewName'], axis = 1)
    prot_stat=df.groupby('station_name')['Number of Proteins'].sum()
    lib=GetLibSize()
    prot_stat=pd.merge(prot_stat, lib, on='station_name')
    prot_stat = prot_stat.drop(['Sample'], axis = 1)
    #print(prot_stat)
    min = prot_stat['Number of Proteins'].min()
    max = prot_stat['Number of Proteins'].max()
    mean = prot_stat['Number of Proteins'].mean().round(2)
    result_min = prot_stat.nsmallest(5, ['Number of Proteins'])
    #print(result_min)
    result_max = prot_stat.nlargest(5, ['Number of Proteins'])
    #print(result_max)
    print("******************** Number of ORFs in station ranges %s - %s *****************" % (str(min), str(max)))
    print("******************* The mean number of ORFs in statin is %s *************" % str(mean))
    #print(prot_stat.sort_values(by='Number of Proteins'))
    #lib_num=prot_stat.plot.scatter(x = "Size", y = "Number of Proteins", subplots=True)
    svg_name = "Lib_Proteins" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Lib_Proteins" + str(version) + '.png'
    png_dir=f'{visuals}/{png_name}'
    #plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #stat_num = prot_stat.plot.scatter(x = "St_Depth", y = "Number of Proteins", subplots = True)
    plt.xticks(rotation = 90)
    svg_name = "Stat_Proteins" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Stat_Proteins" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    g = sns.PairGrid(prot_stat, y_vars="Number of Proteins", x_vars=["Size", "St_Depth"],height=4)
    g.map(sns.scatterplot)
    plt.xticks(rotation = 90)
    svg_name = "Stat_Proteins_dub" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Stat_Proteins_dub" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    num_hist = sns.histplot(prot_stat, x = "Number of Proteins", bins = 30)
    plt.axvline(x = prot_stat['Number of Proteins'].median(),
                color = 'blue',
                ls = '--',
                lw = 1.0)
    min_ylim, max_ylim = plt.ylim()
    plt.text(prot_stat['Number of Proteins'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(prot_stat['Number of Proteins'].median()))
    plt.axvline(x = prot_stat['Number of Proteins'].mean(),
                color = 'red',
                lw = 1.0)
    plt.text(prot_stat['Number of Proteins'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(prot_stat['Number of Proteins'].mean()))
    svg_name = "ProteinsStatHisto" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinsStatHisto" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    """df_min = df_grouped.loc[df_grouped["Number of Proteins"] < 100]
    num_hist = sns.histplot(df_min, x = "Number of Proteins", bins = 30)
    plt.axvline(x = df_min['Number of Proteins'].median(),
                color = 'blue',
                ls = '--',
                lw = 1.5)
    min_ylim, max_ylim = plt.ylim()
    plt.text(df_min['Number of Proteins'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(df_min['Number of Proteins'].median()))
    plt.axvline(x = df_min['Number of Proteins'].mean(),
                color = 'red', lw = 1.5)
    plt.text(df_min['Number of Proteins'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(df_min['Number of Proteins'].mean()))
    svg_name = "ProteinsHistoMin" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinsHistoMin" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    # plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    # plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()"""

def Plasmid_Station():
    ' plasmid per station statistics '
    df = ORF_byPlasmid_stats()
    df_stat=DF_plasmids_byReads()[0]
    df = pd.merge(df, df_stat, left_on = 'Plasmids', right_on = 'NewName')
    df = df.drop(['NewName',"Number of Proteins"], axis = 1)
    lib = GetLibSize()
    plasm_stat = pd.merge(df, lib, on = 'station_name')
    plasm_stat = plasm_stat.drop(['Sample'], axis = 1)
    df_grouped = plasm_stat.groupby(['station_name','Size','St_Depth'])['Plasmids'].nunique().to_frame(name = 'Number of Plasmids').reset_index()
    #print(df_grouped)
    PlamidNum_array = df_grouped['Number of Plasmids'].to_numpy()
    LibSize_array = df_grouped['Size'].to_numpy()
    pearson_lib = stats.pearsonr(PlamidNum_array, LibSize_array)
    if pearson_lib[0] > 0:
        print("The larger library size, more plasmids we find there (r=%f, p-value=%f)" % (
        pearson_lib[0].round(3), pearson_lib[1].round(3)))
    elif pearson_lib[0] == 0:
        print("Library size does not correlate with the number of plasmids at the station (r=%f, p-value=%f)" % (
            pearson_lib[0].round(3), pearson_lib[1].round(3)))
    else:
        print("The larger library size, less plasmids we find there (r=%f, p-value=%f)" % (
        pearson_lib[0].round(3), pearson_lib[1].round(3)))
    slope, intercept, r_value, p_value, std_err = stats.linregress(df_grouped["Number of Plasmids"], df_grouped["Size"])
    print("slope=%f; intercept=%f; r_value=%f; p_value=%f, str_err=%f" % (
    slope.round(3), intercept.round(3), r_value.round(3), p_value.round(3), std_err.round(3)))

    g = sns.PairGrid(df_grouped, y_vars = "Number of Plasmids", x_vars = ["Size", "St_Depth"], height = 4)
    g.map(sns.scatterplot)
    plt.xticks(rotation = 90)
    svg_name = "Stat_Plasmids_dub" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Stat_Plasmids_dub" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    plasmid_hist = sns.histplot(df_grouped, x = "Number of Plasmids")
    plt.axvline(x = df_grouped['Number of Plasmids'].median(),
                color = 'blue',
                ls = '--',
                lw = 1.5)
    min_ylim, max_ylim = plt.ylim()
    plt.text(df_grouped['Number of Plasmids'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(df_grouped['Number of Plasmids'].median()))
    plt.axvline(x = df_grouped['Number of Plasmids'].mean(),
                color = 'red', lw = 1.5)
    plt.text(df_grouped['Number of Plasmids'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(df_grouped['Number of Plasmids'].mean()))
    svg_name = "Stat_PlasmidsHisto" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Stat_PlasmidsHisto" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()


#ORF_byStation_stats()
#Plasmid_Station()
#plasmids_byreads=DF_plasmids_byReads()[1]
#print(plasmids_byreads)