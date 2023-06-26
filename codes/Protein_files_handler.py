# -*- coding: utf-8 -*-
"""
Created on 29/08/2021 11:57

Author: Lucy Androsiuk
"""
### Description
# don't forget description!
import os, glob
import pandas as pd
from pathlib import Path

# uncomment relevant path to OS
# Windows
#path = r"C:\Users\Lucy\iCloudDrive\Documents/bengurion/Plasmidome"
# macOS
path = r"../Output"

# working directories
full_path = f'{path}/blastp'
dataset = r'../res/dataset'
Path(dataset).mkdir(parents=True, exist_ok=True)

# working files

def Blast_file_parser():
    for path in os.listdir(full_path):
        new_path = f'{full_path}/{path}'
        new_file = path + '.csv'
        new_dir=f'{dataset}/{new_file}'
        for file in os.listdir(new_path):
            path_file = f'{new_path}/{file}'
            if not os.stat(path_file).st_size == 0 and os.stat(new_dir).st_size == 0:
                print("I'm parsing files in %s directory. This is file %s." % (path, file))
                df = pd.read_csv(path_file, sep = '\t')
                df.to_csv(new_dir, mode = 'a', header = False, index=False,sep = '\t')
            else:
                print("%s is empty. OR %s is already written" % (file,new_file))
        print("I read all files in %s" % path)
        print("I've merged all files in %s" % path)
        #print(df_concat.head(5))

