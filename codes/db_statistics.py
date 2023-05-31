# -*- coding: utf-8 -*-
"""
Created on 18/11/2021 17:19

Author: Lucy

Description: analyzing PLSDB alignment result to identify wether the candidate is known or novel.
"""
import pandas as pd
from Bio import SeqIO
import re

# directories
output_dir = r"../Output"
plsdb_plasmids=f"{output_dir}/plsdb_new.csv"
plasmid_candidates=f"../res/filtered_plasmids.fasta"

colnames = ['qseqid', 'sseqid', 'stitle', 'evalue', 'length', 'pident', 'mismatch', 'score', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq']

ap_tresh = 50
pi_tresh = 70

pd.set_option('display.max_colwidth', None)

#get length of read
def Clean_length(node_name):
    if re.search('h\_\d+', node_name):
        pos = int(re.search('h\_\d+', node_name).group(0)[2:])
        return pos
    else:
        return node_name

def PlasmidList():
    with open(plasmid_candidates) as all_plasdmids:
            non_unique_list = [] #list for plasmids
            for record in SeqIO.parse(all_plasdmids, "fasta"):
                non_unique_list.append(record.id)
    return non_unique_list

def FilterOutKnown(db_file):
    data = pd.read_csv(db_file, sep = '\t', index_col=None, header = None)
    data.columns = colnames
    data['query_length'] = data['qseqid'].apply(Clean_length)
    data[['query_length', 'length']] = data[["query_length", "length"]].apply(pd.to_numeric)
    data['align_perc'] = (data['length'] / data['query_length']) * 100
    data_filtered = data[(data['pident'] >= pi_tresh) & (data['align_perc'] >= ap_tresh)]
    #data_filtered = data[(data['pident'] >= pi_tresh)]
    data_filtered = data_filtered[(data_filtered['qseqid'] != data_filtered['qseqid'].shift()) | (data_filtered['stitle'] != data_filtered['stitle'].shift())]
    print(data_filtered[['qseqid', 'stitle', 'length', 'pident', 'align_perc']])
    recs = PlasmidList()
    recsSet = set(recs)
    not_novel = []
    data_filtered['qseqid'].apply(lambda x: not_novel.append(x) if x in recsSet else 0)
    print(not_novel)
    return (data_filtered, (list(dict.fromkeys(not_novel))))

FilterOutKnown(plsdb_plasmids)