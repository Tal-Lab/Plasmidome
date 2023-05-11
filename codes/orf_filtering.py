# -*- coding: utf-8 -*-
"""
Created on 07/08/2022 13:10

Author: Lucy Androsiuk
"""
### Description
# In this code ORFs from concatenated (doubled) plasmid candidates
# are being filtered and rewritten into new fasta file with unique ORFs

import random
import re,os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import multiprocessing as mp, time
from itertools import combinations
from pathlib import Path

pd.set_option('max_columns', None)
pd.set_option('max_rows', None)
pd.set_option('display.max_colwidth', None)

version = 1
# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
path = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome"

# working directories
tables = f"{path}/data_calculations"
visuals = f"{path}/visualisations"
Path(tables).mkdir(parents=True, exist_ok=True)
#proteins_un = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/proteins_single.faa"
proteins_dub = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/proteins_double.faa"
length_dub = r"/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/lengths_double.txt"
out_file2 = r'/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/Unique_ORFs_region4.fasta'
filtered_orf_file = r'/Users/lucyandrosiuk/Documents/bengurion/Plasmidome/Filtered_ORFs4.fasta'

# working files
proteins_dub = r"../res/proteins_double.faa"
length_dub = r"../res/lengths_double.txt"

#output file (will be used further as resource)
filtered_orf_file = r'../res/Filtered_ORFs.fasta'

def LenReader():
    """ Getting region of the doubled plasmid """
    df = pd.read_csv(length_dub, sep='\t')
    # print(df.head(10))
    df['start'] = 1
    # df['start']=df['start'].round(0).astype(int)
    df['end'] = df['#name'].apply(lambda x: re.search(r'\d+$', x).group(0)).astype(int)
    # df['end']=df['end'].round(0).astype(int)
    print(df)
    return df

def ProtReader():
    """ Obtains each protein from doubled protein fasta and writes to dataframe """
    records = []
    for record in SeqIO.parse(proteins_dub, 'fasta'):
        plasmid_name = re.search(r'\w+_', record.description).group(0)[:-1]
        protein_name = re.search(r'\w+', record.description).group(0)
        protein_num = re.search(r'_\d+ ', record.description).group(0)[1:-1]
        regex = r'(?<=\#\s)\d+(?=\s\#)'
        position_match = re.findall(regex, record.description)
        start_pos = position_match[0]
        end_pos = position_match[1]
        strand = re.search(r'.\d # I', record.description).group(0)[:-4]
        partial = re.search(r'\d+;s', record.description).group(0)[:-2]
        seq_ce = record.seq
        sequence = ''.join(seq_ce)
        orf_length = len(record.seq)
        records.append([plasmid_name, protein_name, protein_num, orf_length, strand, sequence, start_pos, end_pos, partial])
    df = pd.DataFrame(records, columns=['Plasmid', 'ORF', 'ORF_num', 'ORF_length', 'Strand', 'Sequence', 'Start', 'End', 'Partial'])
    #print(df.shape)
    return df

def Sequence_trunc(seq_ce, orf_end, orf_start):
    orf_length = round((int(orf_end)-int(orf_start))/3)
    sequence = seq_ce[:orf_length]
    return sequence

def Choose_byRegion():
    """ Filters out proteins which start at the plasmid region identified in LenReader(), gets each ORF length and filters out proteins, which are longer than plasmid """
    region = LenReader()
    orfs = ProtReader()
    # orfs_grouped=orfs.groupby('Plasmid')
    df_merged = orfs.merge(region, left_on='Plasmid', right_on='#name')
    df_merged.drop('#name', axis=1, inplace=True)
    # print(df_merged)
    orfs_filtered = df_merged.loc[df_merged['Start'].astype(int).between((df_merged['start']), (df_merged['end']))]
    orfs_filtered = orfs_filtered[orfs_filtered.Plasmid != '94_LNODE_1_length_975']
    print(orfs_filtered.shape)
    orfs_filtered['ORF_lengthD'] = orfs_filtered['ORF_length'].astype(int)*3
    orfs_filtered1 = orfs_filtered.loc[orfs_filtered['ORF_lengthD'] <= (orfs_filtered['end'])]
    orfs_filtered2 = orfs_filtered.loc[orfs_filtered['ORF_lengthD'] > (orfs_filtered['end'])]
    print(orfs_filtered2)
    orfs_filtered2['End'] = orfs_filtered2['end'].astype(int)
    orfs_filtered2['Sequence'] = orfs_filtered2.apply(lambda x: Sequence_trunc(x['Sequence'], x['End'], x['Start']), axis=1)
    orfs_filtered2['ORF_lengthD'] = orfs_filtered2['ORF_length'].astype(int)*3
    orfs_filtered2['ORF_length'] = orfs_filtered2['Sequence'].apply(lambda x: len(x))
    orfs_filtered_n = pd.concat([orfs_filtered1, orfs_filtered2])
    #print(orfs_filtered_n.shape)
    grouped_df = orfs_filtered_n.groupby('Plasmid')['ORF'].count().sort_values().to_frame(name='Count')
    #print(grouped_df[grouped_df['Count'] == 1].shape)
    return orfs_filtered_n, orfs_filtered_n['ORF'].to_list()

def overlap(r1, r2):
    if r1[0] < r1[1]:
        if r2[0] < r1[1]:
            over = r1[1] - r2[0]
            size1 = r1[2]
            size2 = r2[2]
            if over in range(min(size1, size2), max(size1, size2)):
                return over
            else:
                return

def Pairs(list):
    ''' getting list of overlap of each hit '''
    a = []
    for c1, c2 in combinations(list, 2):
        o = overlap(c1, c2)
        if o is not None:
            a.append((o, c1, c2))
    return a

def Interval(df, x, y, z):
    ''' Function calculates overlaps of specified orf (df)'''
    df_int = df[[x, y, z]]
    df_int['Start'] = df_int['Start'].astype(int)
    df_int['End'] = df_int['End'].astype(int)
    df_int['ORF_lengthD'] = df_int['ORF_lengthD'].astype(int)
    # print(df_int)
    interval_list = df_int.values.tolist()
    # print(interval_list)
    ol_list = Pairs(interval_list)
    return ol_list

def CheckOverlappingORFs():
    ''' Function checks each plasmids ORFs for overlapping ORFs. '''
    # print(" ################### CheckOverlappingORFs  ################### ")
    df = Choose_byRegion()[0]
    df['Plasmid_length'] = df['Plasmid'].apply(lambda x: re.search(r'\d+$', x).group(0))
    #df['ORF_length'] = df['End'].astype(int) - df['Start'].astype(int)
    df.loc[df['Start'].astype(int) > df['Plasmid_length'].astype(int), 'Start'] = df['Start'].astype(int) - df[
        'Plasmid_length'].astype(int)
    df.loc[df['End'].astype(int) > df['Plasmid_length'].astype(int), 'End'] = df['End'].astype(int) - df[
        'Plasmid_length'].astype(int)
    #print(df['Plasmid'].nunique())
    #print(df.shape)
    df_grouped = df.groupby('Plasmid')[
        ['Plasmid_length', 'ORF', 'ORF_num', 'ORF_length', 'Strand', 'Sequence', 'Start', 'End', 'Partial', 'ORF_lengthD']]
    filtered_ORF = pd.DataFrame()
    for name, group in df_grouped:
        print('%%%%%%%%%%%%%%%%%%%%%%')
        print(group)
        group['Start'] = group['Start'].astype(int)
        group['End'] = group['End'].astype(int)
        group_nodupl = group.drop_duplicates(subset=['Start', 'End'], keep='last')
        ref_row = group_nodupl.loc[group_nodupl['Start'].astype(int) > group_nodupl['End'].astype(int)]
        if not ref_row.empty:
            #print("############# Printing reference row #############")
            #print(ref_row)
            ref_Start = ref_row['Start'].values[0]
            ref_End = ref_row['End'].values[0]
            #print("############# Printing group 1: const Start #############")
            group1 = group_nodupl.loc[group_nodupl['Start'].astype(int) == int(ref_Start)]
            #print(group1)
            #print("############# Printing group2: const End #############")
            group2 = group_nodupl.loc[group_nodupl['End'].astype(int) == int(ref_End)]
            #print(group2)
            if (group1.shape[0] != 1 or group2.shape[0] != 1) and (group1['Partial'].any() != 00 or group2['Partial'].any() != 00) and (group1['Strand'].nunique()==1 & group2['Strand'].nunique()==1):
                #print("############# Printing dataframe to drop #############")
                df_to_drop = pd.concat([group1, group2]).drop_duplicates(keep=False)
                #print(df_to_drop)
                #print("############# Printing clean dataframe #############")
                group_nodupl = pd.concat([group_nodupl, df_to_drop]).drop_duplicates(keep=False)
                #print(group_nodupl)
        group_nodupl['Start'] = group_nodupl['Start'].astype(int)
        group_nodupl = group_nodupl.sort_values(by=['Start'])
        #print(group_nodupl)
        group3 = group_nodupl[group_nodupl['Strand'] == ' 1']
        if not group3.empty:
            #print('################ Group 3 #############')
            #print(group3)
            overlap1 = Interval(group3, 'Start', 'End', 'ORF_lengthD')
            if len(overlap1) != 0:
                # print('################ Overlap 1 #############')
                # print(overlap1)
                for ol in overlap1:
                    coordinate1 = ol[1]
                    group_c1 = group3.loc[(group3['Start'] == coordinate1[0]) & (group3['End'] == coordinate1[1])]
                    coordinate2 = ol[2]
                    group_c2 = group3.loc[(group3['Start'] == coordinate2[0]) & (group3['End'] == coordinate2[1])]
                    #print('################ Group 3 to drop #############')
                    group_del = pd.concat([group_c1, group_c2])
                    #print(group_del)
                    group_to_del = group_del.loc[group_del['ORF_lengthD'] == group_del['ORF_lengthD'].min()]
                    group3 = pd.concat([group3, group_to_del]).drop_duplicates(keep=False)
        group4 = group_nodupl[group_nodupl['Strand'] == '-1']
        if not group4.empty:
            #print('################ Group 4 #############')
            #print(group4)
            overlap2 = Interval(group4, 'Start', 'End', 'ORF_lengthD')
            if len(overlap2) != 0:
                # print('################ Overlap 2 #############')
                # print(overlap2)
                for ol in overlap2:
                    coordinate3 = ol[1]
                    group_c3 = group4.loc[(group4['Start'] == coordinate3[0]) & (group4['End'] == coordinate3[1])]
                    coordinate4 = ol[2]
                    group_c4 = group4.loc[(group4['Start'] == coordinate4[0]) & (group4['End'] == coordinate4[1])]
                    #print('################ Group 4 to drop #############')
                    group_del2 = pd.concat([group_c3, group_c4])
                    #print(group_del2)
                    group_to_del2 = group_del2.loc[group_del2['ORF_lengthD'] == group_del2['ORF_lengthD'].min()]
                    group4 = pd.concat([group4, group_to_del2]).drop_duplicates(keep=False)
        group_nodupl2 = pd.concat([group3, group4])
        #print(group_nodupl2)
        orf_num = len(group_nodupl2['ORF'].to_list())
        group_nodupl2['ORF_num'] = range(1, orf_num + 1)
        group_nodupl2['ORF'] = group_nodupl2['Plasmid'] + '_' + group_nodupl2['ORF_num'].astype(str)
        filtered_ORF = filtered_ORF.append(group_nodupl2)
    filtered_ORF = filtered_ORF[filtered_ORF.Plasmid != '94_LNODE_1_length_975']
    print(filtered_ORF.shape)
    print(filtered_ORF['Plasmid'].nunique())
    grouped_df = filtered_ORF.groupby('Plasmid')['ORF'].count().sort_values().to_frame(name='Count')
    orf1 = grouped_df[grouped_df['Count'] == 1]
    print(orf1)
    min = grouped_df['Count'].min()
    max = grouped_df['Count'].max()
    print("Number of ORFs in plasmids ranges %s - %s" % (str(min), str(max)))
    mean = grouped_df['Count'].mean().round(2)
    print("The mean number of ORFs in plasmid is %s" % (str(mean)))
    median = grouped_df['Count'].median().round(2)
    print("The median number of ORFs in plasmid is %s" % (str(median)))
    return filtered_ORF
CheckOverlappingORFs()
def FilteredORF_toFasta():
    df = CheckOverlappingORFs()
    ids = df['ORF'].to_list()
    start = 'Start: ' + df['Start'].astype(str)
    start_list = list(start)
    end = 'End: ' + df['End'].astype(str)
    end_list = list(end)
    partial = list('Partial: ' + df['Partial'].astype(str))
    strand = ('Strand: ' + df['Strand'].astype(str)).to_list()
    sequence = df['Sequence'].to_list()
    zipped = list(zip(ids, strand, start_list, end_list, partial, sequence))
    records = []
    for record in zipped:
        id = record[0]
        description = ' | '.join(record[1:5])
        sequence = record[5]
        the_record = SeqRecord(Seq(sequence), id=id, description=description)
        records.append(the_record)
    if not os.path.isfile(filtered_orf_file) or os.stat(filtered_orf_file).st_size == 0:
        SeqIO.write(records, filtered_orf_file, "fasta")

#FilteredORF_toFasta()