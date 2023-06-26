# -*- coding: utf-8 -*-
"""
Created on 29/12/2021 17:05

Author: Lucy Androsiuk
"""
### Description
# add description

import pandas as pd
import re, os, sys, math
import dotenv_setup
import numpy as np
from functools import reduce
import multiprocessing as mp

#pd.set_option('display.max_colwidth', None)
#pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.options.mode.chained_assignment = None

path = r"../Output"

# working directories
visuals = f"{path}/visualisations"
tables = f"{path}/data_calculations"
Path(visuals).mkdir(parents=True, exist_ok=True)
Path(tables).mkdir(parents=True, exist_ok=True)

# working files
aclame_blast = r"../res/dataset/ACLAMEproteins.csv"
aclame_annot = r"../res/Annotations/aclame_family_plasmids_0.4.xls"
ips_prot = r"../res/FilteredORFs.tsv"
egg_nog = r"../res/eggnog_FilteredORFs.csv"
pVerify = r"../res/plasmidVerify_result_table.csv"
vVerify = r"../res/viralVerify_result_table.csv"

colnames = os.getenv('COLS_BLAST')

def Annot_handler():
    df_annot = pd.read_excel(aclame_annot, header = 1)
    df_annot.drop(['Nb proteins'], axis = 1, inplace = True)
    df_annot.set_axis(['family_id', 'function'], axis = 1, inplace = True)
    df_annot = df_annot.fillna('missing')
    s_dict = dict(df_annot.values)
    return df_annot, s_dict

def funct (x):
    '''searching pattern'''
    if re.search(r'\|([^;^|]*)\|', x):
        return re.search(r'\|([^;^|]*)\|', x).group(0)[1:-1]
    else:
        return x

def Aclame_map():
    df_orf = pd.read_csv(aclame_blast, sep = '\t', index_col = None, header = None)
    df_orf.columns = colnames
    df_orf['family'] = df_orf['stitle'].apply(lambda x: re.search(r'family:\w+:\d+', x).group(0))
    df_orf.drop_duplicates(['qseqid', 'family'], inplace=True)
    df_orf2 = df_orf[['qseqid', 'family']]
    df_annot = Annot_handler()[0]
    df_orf3 = df_orf2.merge(df_annot, how='left', left_on='family', right_on='family_id').drop(['family', 'family_id'], axis=1)
    #print(df_orf3['function'])
    df_orf3['func'] = df_orf3['function'].apply(funct)
    df_orf3 = df_orf3.drop('function', axis=1)
    df_orf3 = df_orf3.groupby('qseqid')['func'].apply(','.join).reset_index()
    df_orf3 = df_orf3.rename(columns={"qseqid": "ORF_name", 'func': 'ACLAME_function'})
    #df_orf3['Plasmid'] = df_orf3['qseqid'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    return df_orf3

def IPS_reader():
    df = pd.read_csv(ips_prot, sep='\t', header=None, index_col=None)
    col_names = ['ORF_name', 'Sequence_MD5_digest', 'Sequence_length', 'database', 'Signature_accession', 'Signature_description',
                 'start', 'stop', 'e_value', 'Status', 'Date', 'IP_annot_acs', 'IP_annot_description']
    df.columns = col_names
    df = df.replace({np.nan: ''})
    df = df[['ORF_name', 'database', 'Signature_description', 'IP_annot_description']]
    df['IPS_function'] = df[['Signature_description', 'IP_annot_description']].apply(lambda x: ','.join(x), axis=1)
    df['IPS_function'] = df['IPS_function'].str.split(',').apply(set).str.join(',')
    df = df[['ORF_name', 'database', 'IPS_function']]
    df = df.loc[(df['database'] != 'Coils') & (df['database'] != 'Pfam')]
    df = df.drop_duplicates(["ORF_name", "database"])
    df_pivot = df.pivot(index='ORF_name', columns='database', values='IPS_function')
    df_pivot = df_pivot.fillna('missing')
    df_pivot = df_pivot.add_suffix('_IPS')
    #print(df_pivot)
    #df = df.groupby('ORF_name')['IPS_function'].apply(','.join).reset_index()
    #df['Plasmid'] = df['ORF_name'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    return df_pivot
#IPS_reader()
def EGGNOG_reader():
    df = pd.read_csv(egg_nog, header=0)
    df = df[['#query', 'COG_category', 'Description', 'PFAMs']]
    df = df.rename(columns={"#query": "ORF_name", 'Description': 'EggNOG_function'})
    #df['Plasmid'] = df['ORF_name'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    return df

def PVerify_reader():
    df = pd.read_csv(pVerify)
    df['Plasmid'] = df['Contig name'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    df = df.drop(['Contig name', 'Log-likelihood ratio'], axis=1)
    df = df.rename(columns={"Prediction": "PV_Prediction", 'Predicted HMMs': 'PV_HMMs'})
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return df

def VVerify_reader():
    df = pd.read_csv(vVerify)
    df['Plasmid'] = df['Contig name'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    df = df.drop(['Contig name', 'Length', 'Circular', 'Score'], axis=1)
    df = df.rename(columns={"Prediction": "VV_Prediction", 'Pfam hits': 'VV_Pfam'})
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return df

def Function_ORF():
    df_aclame = Aclame_map()
    df_ips = IPS_reader()
    df_eggnog = EGGNOG_reader()
    dfs = [df_aclame, df_ips, df_eggnog]
    df_final = reduce(lambda left, right: pd.merge(left, right, how='outer', on='ORF_name'), dfs)
    df_final['Plasmid'] = df_final['ORF_name'].apply(lambda x: re.search(r'\w+_l', x).group(0)[:-2])
    cols = df_final.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_final = df_final[cols]
    pl_verify = PVerify_reader()
    vi_verify = VVerify_reader()
    df_final = pd.merge(df_final, pl_verify, how='outer', on='Plasmid')
    df_final = pd.merge(df_final, vi_verify, how='outer', on='Plasmid')
    df_final = df_final[df_final.Plasmid != '94_LNODE_1']
    df_final = df_final.sort_values('ORF_name', ascending=True).reset_index(drop=True)
    #print(df_final['Plasmid'].unique())
    #df_final.to_csv(f'{path}/plasmid_ORF_functions.csv')
    return df_final

def splitter(x):
    if not type(x) == float:
        return ','.join(i for i in x.split(' '))
    else:
        return x

def Plasmid_class():
    'Candidates classification in plasmids, putative plasmids and uncertain'
    df = Function_ORF()
    df = df.fillna('missing')

    ### If viralVerify predicted the candidate as ‘Plasmid’, the candidate is classified as plasmid
    df.loc[(df['VV_Prediction'] == 'Plasmid'), 'Class'] = 'Plasmid'
    df['PV_HMMs'] = df['PV_HMMs'].apply(splitter)
    df['VV_Pfam'] = df['VV_Pfam'].apply(splitter)
    df_def_plasmid = df.loc[df['Class'] == 'Plasmid']

    ### the list of all ORFs functions predictions is collected
    cols = ['ACLAME_function', 'CDD_IPS', 'Gene3D_IPS',
       'Hamap_IPS', 'MobiDBLite_IPS', 'PANTHER_IPS', 'PIRSF_IPS', 'PRINTS_IPS',
       'ProSitePatterns_IPS', 'ProSiteProfiles_IPS', 'SFLD_IPS', 'SMART_IPS',
       'SUPERFAMILY_IPS', 'TIGRFAM_IPS', 'EggNOG_function', 'PFAMs', 'PV_HMMs', 'VV_Pfam']
    params = []
    for col in cols:
        lst = df_def_plasmid[col].to_list()
        for els in lst:
            if not type(els) == float:
                el = els.split(',')
                params.append(el)
            else:
                pass
    params = [item for sublist in params for item in sublist]
    params = list(dict.fromkeys(params))
    # print(params)

    ### to detect the list of plasmid-associated ORFs functions predictions the all ORFs list was manually curated
    plasm_params = ['relaxase activity', 'phage DNA replication', 'plasmid function unknown', 'post-segregational killing', 'regulation of cell division', 'plasmid vegetative DNA replication', 'type IV secretion system coupling protein', 'intercellular transfer by conjugation', 'site-specific DNA excision', 'transpositional DNA recombination', 'site-specific DNA recombination', 'DNA restriction-modification system', 'DNA-methyltransferase activity', 'plasmid mobilization', 'MobA/MobL family', 'MobA/MobL protein', 'Tetracyclin repressor-like', 'Tetracycline Repressor', 'tetR family', 'TetR-type', 'Firmicute plasmid replication protein (RepL)', 'Plasmid replication protein', 'RepL', 'MAPEG family', 'eicosanoid/glutathione metabolism (MAPEG) protein', 'Membrane associated eicosanoid/glutathione metabolism-like domain superfamily', 'MAPEG domain-like', 'PemK-like', 'MazF-like toxin of type II toxin-antitoxin system', 'mRNA interferase PemK-like', 'Plasmid maintenance toxin/Cell growth inhibitor', 'Cell growth inhibitor/plasmid maintenance toxic component', 'Iron dependent repressor', 'DTXR-type HTH domain', 'CobQ/CobB/MinD/ParA nucleotide binding domain', 'ParAB_family', 'Replication initiator protein A', 'Restriction endonuclease-like', 'Restriction endonuclease  type II-like', 'Restriction endonuclease BglII', 'Type-2 restriction enzyme BglII', 'Restriction endonuclease', 'BamHI/BglIII/BstY', 'Type IV secretion-system coupling protein DNA-binding domain', 'Type IV secretion system coupling protein TraD', 'CONJUGATIVE TRANSFER: DNA TRANSPORT', 'TrwC relaxase', 'Origin of replication-binding domain', 'RBD-like', 'SF1_C_RecD', 'DNA2/NAM7 HELICASE FAMILY MEMBER', 'RECBCD ENZYME SUBUNIT RECD', 'relax_trwC: conjugative relaxase domain', 'Conjugative relaxase', 'Resolvase-like', 'DNA-INVERTASE FROM LAMBDOID PROPHAGE', 'Resolvase', 'SR_ResInv', 'resolvase_6', 'Resolvase/invertase-type recombinase catalytic domain profile.', 'SERINE RECOMBINASE PINE-RELATED', 'MT-A70', 'MT-A70-like', 'S-adenosyl-L-methionine-dependent methyltransferases', 'S-adenosyl-L-methionine-dependent methyltransferase', 'MT-A70-like family profile.', 'Vaccinia Virus protein VP39', 'N6-ADENOSINE-METHYLTRANSFERASE', 'COPG FAMILY HELIX-TURN-HELIX PROTEIN-RELATED-RELATED', 'Vibrio phage ICP1', 'Orf50', 'COPG FAMILY HELIX-TURN-HELIX PROTEIN-RELATED', 'Acetyltransferase (GNAT) family', 'GNAT domain', 'GCN5-RELATED N-ACETYLTRANSFERASE', 'L-AMINO ACID N-ACETYLTRANSFERASE', 'Acyl-CoA N-acyltransferases (Nat)', 'Acyl-CoA N-acyltransferase', 'ParB_7', 'ParB/Sulfiredoxin', 'KorB DNA-binding domain-like', 'ParB-like nuclease domain', 'ParB/Sulfiredoxin superfamily', 'ParB_N_like', 'CHROMOSOME 2-PARTITIONING PROTEIN PARB-RELATED', 'CHROMOSOME-PARTITIONING PROTEIN PARB-RELATED', 'Site-specific recombinases signature 2.', 'Recombinase', 'Helix-turn-helix domain of resolvase', 'Site-specific recombinases active site.', 'PIN domain', 'PIN_MtVapC28-VapC30-like', 'PIN domain-like', 'PIN-like domain superfamily', 'Rv0623-like transcription factor', 'Antitoxin VapB-like', 'Tn3 transposase DDE domain', 'PLD-like domain', 'Phospholipase D-like domain', 'Phospholipase D/nuclease', 'Endonuclease Chain A', 'Phospholipase D phosphodiesterase active site profile.', 'Phospholipase D/Transphosphatidylase', 'HsdM N-terminal domain', 'N6 adenine-specific DNA methyltransferase', 'N-6 Adenine-specific DNA methylases signature.', 'DNA methylase', 'N-6 adenine-specific', 'SLR6095 PROTEIN', 'N-6 DNA Methylase', 'adenine-specific', 'N12 class N6 adenine-specific DNA methyltransferase signature', 'TYPE-1 RESTRICTION ENZYME ECOKI SPECIFICITY PROTEIN', 'DNA methylase specificity domains', 'Type I restriction modification DNA specificity domain', 'HsdS', 'DNA methylase specificity domain', 'Bipartite methylase S protein', 'HELICASE SUPERFAMILY 1 AND 2 DOMAIN-CONTAINING PROTEIN', 'Type I restriction enzyme R protein N terminus (HSDR_N)', 'HsdR', 'Superfamilies 1 and 2 helicase ATP-binding type-1 domain profile.', 'Helicase superfamily 1/2', 'ultradead3', 'SWI2/SNF2 ATPase', 'DNA breaking-rejoining enzymes', 'DNA breaking-rejoining enzyme', 'Intergrase catalytic core', 'Integrase-like', 'Tyrosine recombinase domain profile.', 'Integrase', 'DNA_BRE_C', 'Prokaryotic membrane lipoprotein lipid attachment site profile.', 'Bacterial mobilisation protein (MobC)', 'Bacterial mobilisation', 'VirC1 protein', 'TraM recognition site of TraD and TraG', 'TrwC protein', 'Site-specific recombinases DNA invertase Pin homologs', 'ParB domain protein nuclease', 'Transposase', 'Belongs to the N(4) N(6)-methyltransferase family', 'recombinase activity', 'Type I restriction enzyme R Protein', 'MobA_MobL', 'Relaxase', 'Viral_helicase1', 'TetR_C_13', 'TetR_N', 'MAPEG', 'PemK_toxin', 'RPA', 'Endonuc-BglII', 'TrwB_AAD_bind', 'TrwC', 'Acetyltransf_1', 'Acetyltransf_10', 'Acetyltransf_7', 'ParBc', 'PIN', 'PSK_trans_fac', 'DDE_Tnp_Tn3', 'N6_N4_Mtase', 'HsdM_N', 'N6_Mtase', 'Methylase_S', 'HSDR_N', 'ResIII', 'CbiA', 'SWI2_SNF2', 'PLDc_2', 'MobC']

    ### dataframe with candidates, which were not classified yet
    df_ndef_plasmids = df.loc[(df['VV_Prediction'] != 'Plasmid')]
    df_ndef_plasmids['Parameters'] = df_ndef_plasmids[cols].agg(','.join, axis=1)
    df_ndef_plasmids['Parameters'] = df_ndef_plasmids['Parameters'].apply(lambda x: re.sub(r"\B\s+|\s+\B", "", x))
    df_ndef_plasmids['Parameters'] = df_ndef_plasmids['Parameters'].apply(lambda x: ','.join(set(x.split(','))))

    ### the candidate is classified as putative plasmid, if any of the candidate’s ORFs functions predictions includes the word plasmid
    df_ndef_plasmids.loc[df_ndef_plasmids['Parameters'].str.contains('Plasmid', case = False), 'Class'] = 'Putative_plasmid'
    df_plasmid = df[df['Class'] == 'Plasmid'].append(df_ndef_plasmids[df_ndef_plasmids['Class'] == 'Putative_plasmid'])
    # print(df_plasmid['Parameters'].unique())

    ### the list of plasmid-associated ORFs functions predictions from the candidates classified as putative by word plasmid in the function prediction
    new_params = ['VirE_N', 'Thioredoxin', 'PemK_toxin', 'Plasmid maintenance toxin/Cell growth inhibitor',
                  'MazF-like toxin of type II toxin-antitoxin system', 'PemK-like,mRNA interferase PemK-like',
                  'Cell growth inhibitor/plasmid maintenance toxic component', 'ParBc', 'Plasmid replication protein',
                  'RepL', 'Firmicute plasmid replication protein(RepL)', 'Rep_trans', 'Replication initiation factor',
                  'plasmid function unknown', 'relaxase activity', 'Virulence-associated protein E', 'Resolvase',
                  'plasmid vegetative DNA replication']
    plasm_params = plasm_params + new_params
    plasm_params_new = []
    marker2 = set()
    for l in plasm_params:
        ll = l.lower()
        if ll not in marker2:  # test presence
            marker2.add(ll)
            plasm_params_new.append(l)  # preserve order

    ### the candidate is classified as putative plasmid if plasmidVerify predicted the candidate as ‘Plasmid’
    df_ndef_plasmids.loc[df_ndef_plasmids['PV_Prediction'] == 'Plasmid', 'Class'] = 'Putative_plasmid'

    ### the candidate is classified as putative plasmid if plasmidVerify did not classify the candidate as ‘Plasmid’, and any of its ORFs was assigned an ACLAME function
    df_ndef_plasmids.loc[
        (df_ndef_plasmids['PV_Prediction'] != 'Plasmid') & (df_ndef_plasmids['Class'] != 'Plasmid') & df_ndef_plasmids[
            'ACLAME_function'].apply(
            lambda x: all([ele != 'missing' for ele in (x.split(','))])), 'Class'] = 'Putative_plasmid'

    ### the candidate is classified as putative plasmid if plasmidVerify did not classify the candidate as ‘Plasmid’, and any of its ORFs function predictions is plasmid-associated function
    df_ndef_plasmids.loc[
        (df_ndef_plasmids['PV_Prediction'] != 'Plasmid') & (df_ndef_plasmids['Class'] != 'Plasmid') & df_ndef_plasmids[
            'Parameters'].apply(
            lambda x: any([ele in plasm_params_new for ele in (x.split(','))])), 'Class'] = 'Putative_plasmid'

    res_put = df_ndef_plasmids[df_ndef_plasmids['Class'] == 'Putative_plasmid']
    #print(res_put)

    ### list of putative plasmids
    put_plasm_list = res_put['Plasmid'].unique()
    df_ndef_plasmids.loc[df_ndef_plasmids['Plasmid'].isin(put_plasm_list), "Class"] = "Putative_plasmid"

    ### if candidate is not classified as plasmid or putative plasmid, it is classified as uncertain
    df_ndef_plasmids.loc[df_ndef_plasmids['Class'] != 'Putative_plasmid', 'Class'] = 'Uncertain'
    res_unc = df_ndef_plasmids[df_ndef_plasmids['Class'] == 'Uncertain']
    uncert_list = res_unc['Plasmid'].unique()
    df_ndef_plasmids.loc[df_ndef_plasmids['Plasmid'].isin(uncert_list), "Class"] = "Uncertain"
    #df_ndef_plasmids.loc[(df_ndef_plasmids['ORF_name'].str.startswith('3_')) & (df_ndef_plasmids['Class'] == 'Putative_plasmid'), 'Class'] == 'Plasmid'

    ### printing candidates, orf functions and classifications for manual check and reclassification, when appropriate
    #df_ndef_plasmids.to_csv(f'{path}/plasmid_ndef.csv', index = None)
    df_ndef_plasmids = df_ndef_plasmids.drop('Parameters', axis=1)
    result = pd.concat([df_def_plasmid, df_ndef_plasmids])
    # print(result)

    ### if manual curation revealed falsely classified candidates, the true classification is added manually
    result['Class'].mask(result['Plasmid'] == "3_LNODE_1", 'Plasmid', inplace = True)
    #print(result[result['Class']=='Plasmid'])

    ### dataframe with candidates classified as plasmids and their number
    res_plasmid = result[result['Class'] == 'Plasmid']
    print("#### Number of plasmids candidates classified as plasmids: %d" %res_plasmid['Plasmid'].nunique())

    ### dataframe with candidates classified as putative plasmids and their number
    res_put = result[result['Class'] == 'Putative_plasmid']
    print("#### Number of plasmids candidates classified as putative plasmids: %d" %res_put['Plasmid'].nunique())

    ### dataframe with candidates classified as uncertain and their number
    res_unc = result[result['Class'] == 'Uncertain']
    print("#### Number of plasmids candidates classified as uncertain: %d" %res_unc['Plasmid'].nunique())

    ### dataframe with candidates classified as plasmids and putative plasmids and their number
    res_plasmid_put = result[(result['Class'] == 'Plasmid') | (result['Class'] == 'Putative_plasmid')]
    print("#### Number of plasmids candidates classified as plasmids and putative plasmids: %d" % res_plasmid_put['Plasmid'].nunique())

    ### dataframe with all candidates classifications and their number
    re_plasmid_put_unc = result[(result['Class'] == 'Plasmid') | (result['Class'] == 'Putative_plasmid') | (result['Class'] == 'Uncertain')]
    print("#### Number of plasmids candidates: %d" % re_plasmid_put_unc['Plasmid'].nunique())
    result.to_csv(f'{path}/plasmid_classified.csv', index=None)
    return res_plasmid, res_plasmid_put, re_plasmid_put_unc

#Plasmid_class()