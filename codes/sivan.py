# -*- coding: utf-8 -*-
"""
Created on 12/04/2022 11:36

Author: Lucy
"""
import Bio.GenBank
from Bio import SeqIO
import pandas as pd

phylum_df = pd.DataFrame(columns=["id", "phylum"])
print(phylum_df)

with open(r"C:\Users\Lucy\Documents\OneDrive - Israel Oceanograpic & Limnological Research\plasmids\Project_Sivan\CRISPR\bacteria.wgs_mstr.gbff") as handle:
    print("I have read the file")
    for seq_record in Bio.GenBank.parse(handle):
        print(seq_record.accession)
        phylum_df = pd.concat([phylum_df, pd.DataFrame({"id": seq_record.accession[0], "phylum": seq_record.taxonomy[2]}, index = [1])])

print(phylum_df.shape)