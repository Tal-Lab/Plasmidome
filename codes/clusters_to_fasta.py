# -*- coding: utf-8 -*-
"""
Created on 24/04/2023 13:24

Author: Lucy

Description: From similar candidates clusters obtained from Pipeline_Filtering 
generates a fasta file with filtered candidates.
"""

import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

### working files and directories
resource = r"../res"
all_candidates = f'{resource}/CombinedOutput.fasta' #write a full path to the fasta file all plasmid candidates
new_old_names = r'../new_names.csv' #write a full path to the file with old/new names
filtered_candidates = f'{resource}/filtered_plasmids.fasta'  #a full path to the file with filtered plasmids, which will be in the work_dir

def read_csv():
    # Read CSV file into pandas DataFrame and create dictionary of old to new names
    df = pd.read_csv(new_old_names, sep = ',', header = 0, index_col = None)
    name_map = dict(zip(df["old_name"], df["new_name"])) #check if the column names are correct
    return name_map

def FilteredPlasmids_toFasta():
    name_map = read_csv()
    # Iterate over fasta records and write to output with new names
    with open(filtered_candidates, "w") as out:  # Open the output file for writing
        records_filtered= []
        for record in SeqIO.parse(all_candidates, "fasta"):  # Use SeqIO.parse to read the input FASTA file
            if record.id[:record.id.index('_cov_')] in name_map:  # Check if the current record has a new name
                record.id = name_map[record.id[:record.id.index('_cov_')]]  # Update the record's ID with the new name
                record.name = ""  # Set name and description to empty strings to avoid printing them
                record.description = ""
                records_filtered.append(record)
        SeqIO.write(records_filtered, out, "fasta")  # Write the updated record to the output file

FilteredPlasmids_toFasta()