"""
Created on 01/06/2023

Author: Lucy Androsiuk
"""

import numpy as np
import pandas as pd
import os, re
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import seaborn as sns
from scipy import stats
from plasmid_detect import Plasmid_class, colnames
from general_analysis import ORF_stats
import matplotlib.gridspec as gridspec

path = r"../Output"
candidates_seqs = f'../res/filtered_plasmids.fasta'

# working directories
tables = f"{path}/data_calculations"
visuals = f"{path}/visualisations"
Path(tables).mkdir(parents=True, exist_ok=True)
Path(visuals).mkdir(parents=True, exist_ok=True)

df_class = Plasmid_class()[2]
plasmids = Plasmid_class()[0]['Plasmid'].unique().tolist()
putative_plasmids = Plasmid_class()[1]['Plasmid'].unique().tolist()
all_candidates = Plasmid_class()[2]['Plasmid'].unique().tolist()

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
    out = (df_grouped.merge(df_class, left_on = 'Plasmids', right_on = 'Plasmid')
           .reindex(columns = ['Plasmid', 'Plasmid Length', 'Number of Proteins', 'Class']))
    out.drop_duplicates(subset = None, keep = 'first', inplace = True)
    out.sort_values('Class',inplace = True)
    out.loc[out['Class'] == 'Putative_plasmid', 'Class'] = 'Putative plasmid'
    out.reset_index(inplace = True, drop = True)
    #print(out)
    out['Length_norm'] = out['Plasmid Length'].apply(lambda x: round(x / 10 ** 3))
    pearson_ORF = stats.pearsonr(out["Number of Proteins"].to_numpy(), out["Plasmid Length"].to_numpy())
    #print((pearson_ORF))
    sns.scatterplot(data=out, x="Length_norm", y="Number of Proteins", hue='Class', style='Class')
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
    out_min = out.loc[out['Length_norm'] < 150]
    print(out.loc[out['Length_norm'] >= 150])
    print('Number of plasmids < 150 kb: %s' % str(out_min['Plasmid'].nunique()))
    sns.scatterplot(data = out_min, x = "Length_norm", y = "Number of Proteins", hue = 'Class', style = 'Class')
    plt.xlabel('Plasmid length, kb')
    svg_name = "Len150_plasmids" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Len150_plasmids" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()
    sns.histplot(out, x = "Number of Proteins", bins = 40, hue = 'Class', multiple = 'stack')
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
    df_min=out.loc[out["Number of Proteins"]<=200]
    print('Number of Proteins < 200 per plasmid: %s' % str(df_min['Number of Proteins'].sum()))
    sns.histplot(df_min, x = "Number of Proteins", bins = 40, hue = 'Class', multiple = 'stack')
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
    return df_grouped, out

def Candidates_length():
    """ Statistics for candidates lengths """
    df = ORF_byPlasmid_stats()[1]
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
    # getting plasmid candidates classified as uncertain
    df_uncert = df[df['Class']=='Uncertain']
    pl_uncert = df_uncert['Plasmid'].nunique()
    max_unc = df_uncert['Plasmid Length'].max()
    min_unc = df_uncert['Plasmid Length'].min()
    print("Candidates classified as uncertain are in range from %d to %d bp" % (min_unc, max_unc))
    # getting number of uncertain candidates shorter then 4 kb
    df4kb = df_uncert[df_uncert['Length_norm'] < 4]
    print('Number of plasmids with length < 4 kb is: %d' % df4kb['Plasmid'].nunique())
    print('Percentage of plasmids with length < 4 kb is: %d' % ((df4kb['Plasmid'].nunique() / pl_uncert) * 100))
    # plotting histogram for candidate's lengths, colored by class
    df_plput = df100.loc[df100['Class']!='Uncertain']
    sns.histplot(df_plput, x = 'Length_norm', hue = 'Class', multiple = 'stack', bins = 40)
    #plt.show()


    sns.histplot(df, x='Length_norm', hue = 'Class', multiple = 'stack', bins = 40)
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
    sns.histplot(df_prot_min, x = "Number of Proteins", bins = 40, hue = 'Class', multiple = 'stack', ax = ax1)

    ax2 = fig.add_subplot(gs[1])
    df_min = df.loc[df['Length_norm']<=150]
    # plotting histogram for candidate's lengths (<200bp), colored by class
    sns.histplot(df_min, x = 'Length_norm', hue = 'Class', multiple = 'stack', bins =40, ax = ax2)

    plt.xlabel('Plasmid length, kb')

    fig.tight_layout(pad = 2.0)

    svg_name = "ORFs_per_plasmidBP" + str(version) + '.eps'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "ORFs_per_plasmidBP" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    plt.savefig(svg_dir, format = 'eps', dpi = gcf().dpi, bbox_inches = 'tight')
    plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')

def PieClass():
    df = df_class[['Plasmid', 'Class']].drop_duplicates().reset_index(drop=True)
    df.loc[df['Class']=='Putative_plasmid', 'Class'] = 'Putative plasmid'
    df = df.groupby('Class').count()
    data = df['Plasmid'].to_list()
    label = df.index.to_list()
    colors = sns.color_palette('pastel')[0:3]
    plt.pie(data, labels = label, colors = colors, autopct='%.0f%%')
    svg_name = "Class_piechart" + str(version) + '.svg'
    svg_dir = f'{visuals}/{svg_name}'
    png_name = "Class_piechart" + str(version) + '.png'
    png_dir = f'{visuals}/{png_name}'
    #plt.savefig(svg_dir, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.savefig(png_dir, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')
    #plt.show()

def gc_content():
    # open the FASTA file and iterate over each sequence
    gc_recs = []
    candidates = []
    with open(candidates_seqs) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # calculate the GC content of the sequence
            gc_count = sum(1 for base in record.seq if base in ['G', 'C', 'g', 'c'])
            gc_content = gc_count / len(record.seq) * 100

            # print the ID and GC content of the sequence
            print(f"{record.id}\t{gc_content:.2f}")
            gc_recs.append(round(gc_content,2))
            candidates.append(record.id)
    pl_gc = pd.DataFrame(list(zip(candidates, gc_recs)), columns =['Candidate', 'GC'])
    pl_gc['Candidate'] = pl_gc['Candidate'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    df = df_class[['Plasmid', 'Class']].drop_duplicates().reset_index(drop = True)
    out = pl_gc.merge(df, left_on = 'Candidate', right_on = 'Plasmid')
    out.drop('Plasmid', axis = 1, inplace = True)
    out_put = out.loc[out['Candidate'].isin(putative_plasmids)]
    out_plasm = out.loc[out['Candidate'].isin(plasmids)]
    print(out_plasm)
    print('The GC content values in all candidates ranges between %s and %s' % (str(out['GC'].min()), str(out['GC'].max())))
    print('Average GC content in all candidates is %s' % str(out['GC'].mean()))
    print('The GC content values in putative plasmids ranges between %s and %s' % (
    str(out_put['GC'].min()), str(out_put['GC'].max())))
    print('Average GC content in putative plasmids is %s' % str(out_put['GC'].mean()))
    print('The GC content values in plasmids ranges between %s and %s' % (
    str(out_plasm['GC'].min()), str(out_plasm['GC'].max())))
    print('Average GC content in plasmids is %s' % str(out_plasm['GC'].mean()))
    out_file = f'{tables}/GC_content.csv'
    # print(out_file)
    if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
        out.to_csv(out_file, index = True)

#PieClass()
#gc_content()