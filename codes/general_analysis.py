# -*- coding: utf-8 -*-
"""
Created on 29/08/2021 11:05

Author: Lucy Androsiuk
"""
### Description
# add description

version=8

import numpy as np
import pandas as pd
import os, re, dotenv_setup
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import seaborn as sns
from scipy import stats
import matplotlib.gridspec as gridspec

pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)


path = r"../Output"

# working directories
tables = f"{path}/data_calculations"
visuals = f"{path}/visualisations"
Path(tables).mkdir(parents=True, exist_ok=True)
Path(visuals).mkdir(parents=True, exist_ok=True)

# working files
reads_coverage = r"../res/all_cov.csv"
proteins = r"../res/Filtered_ORFs.fasta"
library = r"../res/LibrarySize.csv"
stations = r"../res/stations.txt"
nt_entries = r"../res/dataset/nt_blast.zip"

colnames = os.getenv('COLS_BLAST')

def CoverageDF():
    ' this function reads coverage file as a dataframe '
    data = pd.read_csv(reads_coverage, sep = ',', index_col = None, header = 0)
    #data = data.reset_index()
    data['rname'] = data['rname'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    data = data.set_index('rname')
    #print([i for i in data.index])
    rows = list(data.index)
    arr = data.loc[rows].values
    return (arr, data)

def PlasmidsbyReads():
    ' this function generates plasmid presence matrix where 0 is plasmid absence (plasmids coverage < 99)    \
    and 1 is plasmid presence in the sampling station (plasmids coverage >= 99) '
    arr, df = CoverageDF()
    df_norm = df
    df_norm[:] = np.where(df_norm < 99, 0, 1)
    #print(df_norm)
    #df_norm = df_norm.append(df_norm.agg(['sum']))
    return df_norm

def DF_plasmids_byReads():
    ' this function generates dataframe with plasmids present in stations and writes it to file '
    df = PlasmidsbyReads()
    df = df.stack().to_frame()
    df = df.reset_index(level=[0,1])
    df.columns = ['NewName', 'station_name', 'value']
    df=df[df['value'] > 0]
    #print(df[['NewName', 'station_name']])
    out_file = f'{tables}/Plasmids_ByReads.csv'
    #print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        df[['NewName', 'station_name']].to_csv(out_file, index = False)
    return df[['NewName', 'station_name']], out_file

plasmids_by_reads=DF_plasmids_byReads()

### getting statistics on plasmids by reads
def Cov_plasmids():
    df = DF_plasmids_byReads()
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

def Coverage_stat():
    ' plasmid presence statistics '
    arr, df = CoverageDF()
    df_norm = df
    df_norm[:] = np.where(df_norm < 99, 0, 1)
    df.loc[:, 'elements'] = df.sum(axis = 1)
    df.sort_values(by = 'elements', ascending = False, inplace = True)
    df2 = df[df['elements'] == 1]
    print("****************** These are plasmids observed at 1 station only *********************")
    print(df.index.unique())
    print(df2.index.nunique())

statistics=Coverage_stat()

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
    #plt.show()
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
    #plt.show()
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
    #print("Plasmids with minimal number of ORFs:")
    result_min = df_grouped[df_grouped['Number of Proteins'] == 1]
    #print(result_min)
    #print(result_min['Plasmids'].nunique())
    #print("Plasmids with the highest number of ORFs:")
    result_max=df_grouped.nlargest(5,['Number of Proteins'])
    #print(result_max)
    #print("******************** Number of ORFs in plasmids ranges %s - %s *****************" % (str(min), str(max)))
    #print("******************* The mean number of ORFs in plasmid is %s *************" % str(mean))
    df_grouped['Plasmid Length']= df_grouped['Plasmids'].apply(Clean_length)
    df_grouped['Plasmids'] = df_grouped['Plasmids'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    #print(df_grouped)

    df_grouped['Length_norm'] = df_grouped['Plasmid Length'].apply(lambda x: round(x / 10 ** 3))
    pearson_ORF = stats.pearsonr(df_grouped["Number of Proteins"].to_numpy(), df_grouped["Plasmid Length"].to_numpy())
    #print((pearson_ORF))
    sns.scatterplot(data=out, x="Length_norm", y="Number of Proteins")
    plt.xlabel('Plasmid length, kb')
    #plt.xticks(rotation = 90)
    #plt.xscale('log')
    svg_name = "Len_Proteins" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Len_Proteins" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()
    print('Total number of plasmids: %s' % str(out['Plasmid'].nunique()))
    out_min = df_grouped.loc[df_grouped['Length_norm'] < 150]
    print(df_grouped.loc[df_grouped['Length_norm'] >= 150])
    print('Number of plasmids < 150 kb: %s' % str(out_min['Plasmid'].nunique()))
    sns.scatterplot(data = out_min, x = "Length_norm", y = "Number of Proteins")
    plt.xlabel('Plasmid length, kb')
    svg_name = "Len150_plasmids" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Len150_plasmids" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()
    sns.histplot(df_grouped, x = "Number of Proteins", bins = 40)
    plt.xlabel('Number of ORFs')
    #sns.histplot(df_grouped_7pl, x = "Number of Proteins", bins = 30, multiple='layer')
    """plt.axvline(x = out['Number of Proteins'].median(),
                color = 'blue',
                ls = '--',
                lw = 1.0)
    min_ylim, max_ylim = plt.ylim()
    plt.text(out['Number of Proteins'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(out['Number of Proteins'].median()))
    plt.axvline(x = out['Number of Proteins'].mean(),
                color = 'red',
                lw = 1.0)
    plt.text(out['Number of Proteins'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(out['Number of Proteins'].mean()))"""
    svg_name = "ProteinsHisto" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinsHisto" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()
    print('Total number of proteins: %s' % str(out['Number of Proteins'].sum()))
    df_min=df_grouped.loc[df_grouped["Number of Proteins"]<=200]
    print('Number of Proteins < 200 per plasmid: %s' % str(df_min['Number of Proteins'].sum()))
    sns.histplot(df_min, x = "Number of Proteins", bins = 40)
    plt.xlabel('Number of ORFs')
    """plt.axvline(x = df_min['Number of Proteins'].median(),
                color = 'blue',
                ls = '--',
                lw = 1.5)
    min_ylim, max_ylim = plt.ylim()
    plt.text(df_min['Number of Proteins'].median() * 1.1, max_ylim * 0.9,
             'Median: {:.2f}'.format(df_min['Number of Proteins'].median()))
    plt.axvline(x = df_min['Number of Proteins'].mean(),
                color = 'red',lw = 1.5)
    plt.text(df_min['Number of Proteins'].mean() * 1.1, max_ylim * 0.8,
             'Mean: {:.2f}'.format(df_min['Number of Proteins'].mean()))"""

    svg_name = "ProteinsHistoMin2" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ProteinsHistoMin2" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    # plt.autoscale()
    plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()
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

def GetLibSize():
    ' retreiving library sizes from lib-file '
    df = pd.read_csv(library, sep = ',', header = None, index_col = None)
    df.columns = ['Sample', 'Size']
    #print(df)
    df['station_name'] = df['Sample'].astype(str) + "_coverage"
    df2 = Station()
    #print(df2)
    df['St_Depth'] = df['station_name'].map(df2.set_index('station_name')['St_Depth'])
    #print(df)
    size_stat = df.plot.scatter(x = "St_Depth", y = "Size")
    plt.xticks(rotation = 90)
    plt.yscale('log')
    #plt.show()
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
    #plt.show()
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
    #plt.show()
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
    #plt.show()
    df_min = prot_stat.loc[prot_stat["Number of Proteins"] < 200]
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
    plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()

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
    #plt.show()
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
    #plt.show()

def Candidates_length():
    """ Statistics for candidates lengths """
    df = ORF_byPlasmid_stats()
    #df['Length_norm'] = df['Plasmid_Length'].apply(lambda x: round(x/10**3))
    print(df.columns)
    # getting the longest plasmid candidate and its length
    df_max = df['Plasmid Length'].max()
    print("The longest candidate length is %d bp" % df_max)
    # calculating number of plasmid candidates with length <100 bp
    df100 = df.loc[df['Length_norm'] <= 150]
    print('Number of plasmids with length <=150 kb is: %d' % df100['Plasmid'].nunique())
    print('Percentage of plasmids with length <=150 kb is: %d' % ((df100['Plasmid'].nunique()/df['Plasmid'].nunique())*100))
    # calculating number of plasmid candidates with length around 140 bp
    df140 = df.loc[df['Length_norm'].between(120,200)]
    print('Number of plasmids around 140 kb is: %d' % df140['Plasmid'].nunique())

    sns.histplot(df, x='Length_norm', bins = 40)
    plt.xlabel('Plasmid length, kb')
    svg_name = "Plasmid_lengths_Histo" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Plasmid_lengths_Histo" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.show()
    # getting candidates <200bp into separate dataframe
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style = 'white', font = 'Helvetica', rc= custom_params)
    gs = gridspec.GridSpec(1, 2)
    fig = plt.figure(figsize = [7.2, 4], dpi = 900)

    # first plot
    ax1 = fig.add_subplot(gs[0])

    df_prot_min = df.loc[df["Number of Proteins"] <= 200]
    print('Number of Proteins < 200 per plasmid: %s' % str(df_prot_min['Number of Proteins'].sum()))
    sns.histplot(df_prot_min, x = "Number of Proteins", bins = 40, ax = ax1)

    ax2 = fig.add_subplot(gs[1])
    df_min = df.loc[df['Length_norm']<=150]
    # plotting histogram for candidate's lengths (<200bp), colored by class
    sns.histplot(df_min, x = 'Length_norm', bins =40, ax = ax2)

    plt.xlabel('Plasmid length, kb')

    fig.tight_layout(pad = 2.0)

    svg_name = "ORFs_per_plasmidBP" + str(version) + '.eps'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ORFs_per_plasmidBP" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    plt.savefig(svg_dir, format = 'eps', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')

    #plt.show()

def nt_counts():
    """Function for nt-database statistics: number of all and significant matches,\
    matches to chromosomal and viral elements; getting number of sampling points,\
    where candidates with no match in this database where found"""
    #creating dataframe from nt-blast file
    df_nt = pd.read_csv(nt_entries, compression='zip', sep = '\t', index_col = None, header = None)
    df_nt.columns = colnames
    #reading dataframe for plasmids stations by reads
    plasmids_by_reads = DF_plasmids_byReads()[0]
    plasmids_by_reads = plasmids_by_reads[plasmids_by_reads.NewName != '94_LNODE_1']
    #extracting plasmid length from the query id
    df_nt['Pl_length'] = df_nt['qseqid'].apply(lambda x: re.search(r'\d+$', x).group(0))
    df_nt['Pl_length'] = df_nt['Pl_length'].astype(int)
    # extracting plasmid name from the query id
    df_nt['Plasmid'] = df_nt['qseqid'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    #counting the number of plasmid candidates with match in nt-database
    n_nt_plasmids = df_nt['Plasmid'].nunique()
    print('Number of plasmid candidates, with %d match in nt-database: %d' % ((round((df_nt['pident'].min()),0)), n_nt_plasmids))
    print('Percentage of plasmid candidates, with %d match in nt-database: %d' % ((round((df_nt['pident'].min()),0)), ((n_nt_plasmids/ plasmids_by_reads['NewName'].nunique()) * 100)))
    # counting the number of plasmid candidates with match to bacterial chromosome in nt-database
    nt_chromosome = df_nt[df_nt['stitle'].str.contains('chromosome')]
    print('Number of plasmid candidates, matched to chromosome: %d' % nt_chromosome['Plasmid'].nunique())
    print('Percentage  of plasmid candidates, matched to chromosome: %d' % ((nt_chromosome['Plasmid'].nunique() / n_nt_plasmids) * 100))
    # counting the number of plasmid candidates with match to viral elements in nt-database
    nt_virus = df_nt[(df_nt['stitle'].str.contains('virus')) | (df_nt['stitle'].str.contains('phage'))]
    print('Number of plasmid candidates, matched to viral elements: %d' % nt_virus['Plasmid'].nunique())
    print('Percentage  of plasmid candidates, matched to viral elements: %d' % ((nt_virus['Plasmid'].nunique() / n_nt_plasmids) * 100))
    # getting candidates matching chromosomal of viral elements
    nt_chromvir = df_nt[df_nt['stitle'].str.contains('chromosome') | (df_nt['stitle'].str.contains('virus')) | (df_nt['stitle'].str.contains('phage'))]
    print('Number of plasmid candidates, matched to chromosomal or viral elements: %d' % nt_chromvir['Plasmid'].nunique())
    print('Percentage  of plasmid candidates, matched to chromosomal or viral elements: %d' % ((nt_chromvir['Plasmid'].nunique() /  plasmids_by_reads['NewName'].nunique()) * 100))
    # getting names of candidates with match in nt-database
    nt_plasmids = df_nt['Plasmid'].unique()
    # print(nt_plasmids)
    # gettiing number of sampling points where each candidate was found
    stat_pl = plasmids_by_reads.groupby('NewName').size().reset_index(name='Number_of_stations')
    # print(stat_pl)
    # getting number of sampling points, where candidates with no match in nt-database where found
    non_nt_pl = stat_pl[~stat_pl['NewName'].isin(nt_plasmids)].reset_index(drop = True)
    # getting candidates, which had no match in nt-database, found in more than 1 sampling point
    non_nt_1st = non_nt_pl[non_nt_pl['Number_of_stations'] > 1]
    print('Number of plasmid candidates, with no significant match in nt-database, appearing more than in 1 sampling point: %d' % non_nt_1st['NewName'].nunique())
    print('Percentage of plasmid candidates, with no significant match in nt-database, appearing more than in 1 sampling point: %d' % ((non_nt_1st['NewName'].nunique() / non_nt_pl['NewName'].nunique()) * 100))
    # counting the number of plasmid candidates with significant (>90% PI and >90% coverage) match in nt-database
    df_nt_high = df_nt[(df_nt['pident']>90) & (df_nt['qcovs']>90)]
    df_nt_high = df_nt_high.loc[abs(df_nt_high['qend']-df_nt_high['qstart']).between((df_nt_high['Pl_length']-round((df_nt_high['Pl_length']*0.1),0)), (df_nt_high['Pl_length']+round((df_nt_high['Pl_length']*0.1),0)))]
    #print(df_nt_high[['Plasmid','stitle','pident','length','qcovs']])
    out_file = f'{tables}/nt_high.csv'
    # print(out_file)
    #writing candidates with significant match to file
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        df_nt_high[['Plasmid','sseqid','stitle','pident','length','mismatch','qcovs','qstart', 'qend', 'sstart', 'send']].to_csv(out_file, index = False)
    print(df_nt_high['Plasmid'].unique())
    print('Number of plasmid candidates, with significant match in nt-database: %d'
          % df_nt_high['Plasmid'].nunique())
    print('Percentage of plasmid candidates, with significant match in nt-database: %d' % ((df_nt_high['Plasmid'].nunique() / plasmids_by_reads['NewName'].nunique()) * 100))

#nt_counts()
#ORF_byPlasmid_stats()
#Candidates_length()
#ORF_byStation_stats()
#Plasmid_Station()
