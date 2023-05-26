# -*- coding: utf-8 -*-
"""
Created on 19/12/2022 13:56

Author: Lucy

Description: Plasmidome Detection pipeline segments I-III
"""

import time
import re
import numpy as np
import os
import sys
import subprocess
import wget
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio import SearchIO
from Bio import SeqIO
import pandas as pd
import zipfile
import tarfile
import gzip
import logging
import shutil
from shutil import copyfile
from csv import writer
from subprocess import Popen,PIPE

### directions
## working directories
outing = r"../"
db_dir = r"../DBs"
resource = r"../res"
logs = r"../logs"
Path(logs).mkdir(parents=True, exist_ok=True)

### working files
sample_matrix = f'{resource}/samples_matrix.txt'
blast_dir = os.getenv("BLAST")
spades_dir = f'{os.getenv("SPADES")}/spades.py'
recycler_dir = f'{os.getenv("RECYCLER")}/recycle.py'
bash_file = r"./RecyclerBwa.csh"

all_data_log = f'{logs}/all_data.log'
Path(db_dir).mkdir(parents=True, exist_ok=True)
lib_file = f'{resource}/library_size.csv'
start_time = time.time()

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

def setup_logger(name, log_file, level=logging.DEBUG):
    """To setup as many loggers as you want"""
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger

def error_tracker():
    error_file = logs + "/Error.log"
    error_logger = setup_logger("error_logging", error_file)
    exc_type, exc_value, traceback = sys.exc_info()
    error_logger.exception(sys.exc_info())
    all_logger.exception(sys.exc_info())
    print("exception caught: %s" % exc_type, exc_value)

def CreateDirectory (name):
    '''Creating directory with sepcified name in the parent_dir'''
    parent_dir = outing
    path = os.path.join(parent_dir, name)
    if not os.path.exists(path):
        os.mkdir(path)
        #all_logger.info("Directory '%s' created" % path)
    return path

all_out = CreateDirectory('Output')

all_logger = setup_logger('second_logger', all_data_log)

def CreateFolderDirectory (name):
    '''Creating folder in in the directory'''
    parent_dir = CreateDirectory(name)
    FolderName = FileNum + name
    pathFolder = os.path.join(parent_dir, FolderName)
    if not os.path.exists(pathFolder):
        os.mkdir(pathFolder)
        all_logger.info("Directory '%s' created" % pathFolder)
    return pathFolder

def GetFileName (link):
    '''Extracting file name from the link'''
    file_regex = re.search("\/S+\w+.fastq.gz", link)
    file_name = file_regex.group(0)[1:]
    #print(file_name)
    #all_logger.info("Getting link %s" % file_name)
    return file_name

def GetFile (link):
    '''Getting reads file from the link'''
    download_time = time.time()
    parent_dir = input_dir + '/'
    child_dir = GetFileName(link)
    pathFile = os.path.join(parent_dir, child_dir)
    #print(pathFile)
    try:
        while not os.path.isfile(pathFile):
            all_logger.info("Starting download %s" % link)
            wget.download(link, pathFile)
            all_logger.warning("Something may have gone wrong here")
            all_logger.info("Finished download %s" % link)
            all_logger.info("Downloading reads took %s seconds" % (time.time()-download_time))
        else:
            all_logger.warning("File was already downloaded")
        return pathFile
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def GetLibrary(path):
    for file in os.listdir(path):
        all_logger.info("***** Getting information for %s " % file)
        if re.match('.*\.fastq.gz', file):
            fullpath = os.path.join(path, file)
            p1 = Popen(["zcat", fullpath], stdout = PIPE)
            p2 = Popen(["wc", "-l"], stdin = p1.stdout, stdout = PIPE)
            p1.stdout.close()
            result = p2.communicate()[0]
            all_logger.info("***** Number of lines in file is: %s" % str(result))
            numseqs = int(result) / 4.0
            all_logger.info("***** Library size is: %s" % str(numseqs))
    shutil.rmtree(path)
    all_logger.info("***** Library size after dropping duplicates: %s" % numseqs)
    return numseqs

def RunSpades():
    '''Running SPAdes to generate assembly_graphs'''
    spades_time = time.time()
    try:
        if not os.path.isfile(output_dir + "/assembly_graph.fastg"):
            subprocess.call(["python", spades_dir,
                             "-1", file_fw,
                             "-2", file_rv,
                             "-o", output_dir])
            all_logger.info("SPAdes running took %s seconds" % (time.time()-spades_time))
        else:
            all_logger.warning("SPAdes has already created 'assembly_graph.fastg'")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def RecyclerBwa():
    '''Running Bwa to work with BAM files for Recycler'''
    in_dir = output_dir + "/assembly_graph.fastg"
    out_dir = output_dir + "/assembly_graph.nodes.fasta"
    try:
        bwa_time= time.time()
        if not os.path.exists(output_dir+"/reads_pe.bam"):
            all_logger.info("Starting preparations for RECYCLER. Convertion of fastg to bam.")
            subprocess.call(["bash", bash_file, in_dir, out_dir, output_dir, file_fw, file_rv, logs])
            all_logger.warning("Check ERRORs file for records")
            all_logger.info("reads_pe.bam file created")
            all_logger.info("Convertion to BAM took %s seconds" % (time.time() - bwa_time))
            assert os.path.isfile(bash_file), "bash file missing"
        else:
            all_logger.warning("Looks like BWA BAM has already created 'reads_pe.bam'. Check ERRORs file for records")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()
    return in_dir

def Recycler():
    ''' Running Recycler to assemble files'''
    b_dir = output_dir + "/reads_pe_primary.sort.bam"
    try:
        recycler_time = time.time()
        if not os.path.isfile(output_dir + "/assembly_graph.cycs.fasta"):
            subprocess.call(["python", recycler_dir,
                             "-g", RecyclerBwa(),
                             "-k", "55",
                             "-b", b_dir])
            all_logger.info("RECYCLER running took %s seconds" % (time.time()-recycler_time))
            all_logger.warning("Something may have gone wrong here. Check file 'assembly_graph.cycs.fasta'.")
        else:
            all_logger.warning("Looks like RECYCLER has already found nodes.Check file 'assembly_graph.cycs.fasta'.")
        files_in_dir = os.listdir(output_dir + '/')  # get list of files in the directory
        for file in files_in_dir:  # loop to delete each unnecessary file in folder
            important_files = ['scaffolds.fasta', 'assembly_graph.cycs.fasta', 'spades.log', 'assembly_graph.fastg', 'metaplasmidSPAdes']
            path_to_file = f'{output_dir}/{file}'
            if os.path.isfile(path_to_file) and file not in important_files:
                all_logger.info("Will not need: %s. Removing it" % file)
                os.remove(f'{output_dir}/{file}')
            elif os.path.isdir(path_to_file) and file not in important_files:
                all_logger.info("Will not need: %s. Removing it" % file)
                shutil.rmtree(f'{output_dir}/{file}')
            elif os.path.isfile(path_to_file) and file in important_files:
                all_logger.warning("Will need the file %s. Saved it in output direcoty" % file)
            else:
                all_logger.warning("No files to delete anymore")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def MetaPlasmidSPAdes():
    ''' Running metaplasmidSPAdes to assemble plasmids'''
    MPspades_path = os.path.join(output_dir, "metaplasmidSPAdes")
    if not os.path.exists(MPspades_path):
        os.mkdir(MPspades_path)
        all_logger.info("Directory ' %s ' created" % MPspades_path)
    try:
        metaplasmid_time = time.time()
        if not os.path.isfile(output_dir + "/metaplasmidSPAdes/scaffolds.fasta"):
            subprocess.call(["python", spades_dir,
                             "--meta", "--plasmid",
                             "-1", file_fw,
                             "-2", file_rv,
                             "-o", MPspades_path])
            all_logger.info("MetaPlasmidSPAdes running took %s seconds" % (time.time()-metaplasmid_time))
            all_logger.warning("Something may have gone wrong here. Check file 'scaffolds.fasta'.")
        else:
            all_logger.warning("Looks like MetaPlasmidSPAdes has already found nodes. Check file 'scaffolds.fasta' and 'spades.log' for errors.")
        files_in_dir = os.listdir(MPspades_path)
        for file in files_in_dir:  # loop to delete each file in folder
            important_files = ['scaffolds.fasta', 'spades.log']
            path_to_file = f'{MPspades_path}/{file}'
            if os.path.isfile(path_to_file) and file not in important_files:
                all_logger.info("Will not need: %s. Removing it" % file)
                os.remove(f'{MPspades_path}/{file}')
            elif os.path.isdir(path_to_file) and file not in important_files:
                all_logger.info("Will not need: %s. Removing it" % file)
                shutil.rmtree(f'{MPspades_path}/{file}')
            elif os.path.isfile(path_to_file) and file in important_files:
                all_logger.warning("Will need the file %s. Saved it in MetaPlasmidSPAdes output direcoty" % file)
            else:
                all_logger.warning("No files to delete anymore")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def alignNodes():
    '''Comparing candidates obtained from Recycler with ones obtained from metaplasmidSPAdes with BLASTn'''
    query_file = output_dir + "/assembly_graph.cycs.fasta"
    subject_file = output_dir + "/metaplasmidSPAdes/scaffolds.fasta"
    newName = "aligned.tab"
    newFile_path = output_dir + "/aligned.tab"
    newFile = open(newName, 'wt')
    try:
        unique_nodes_time = time.time()
        if not os.stat(query_file).st_size==0 and not os.stat(subject_file).st_size==0:
            blastn = blast_dir + "blastn"
            all_logger.info("Starting the comparison of Recycler and MetaPlasmidSPAdes outputs")
            cline = NcbiblastnCommandline(blastn,
                                          query = query_file,
                                          subject = subject_file,
                                          evalue = 0.001,
                                          out = newFile_path,
                                          outfmt = '6 qseqid sseqid evalue length pident score '
                                                         'qcovs qstart qend sstart send qseq sseq')
            all_logger.info(cline)
            os.system(str(cline))
            stdout, stderr = cline()
            newFile.close()
            all_logger.info("Alignment finished and written to the file aligned.tab in your output directory")
            all_logger.info("Alignment of RECYCLER to MPSPAdes output nodes took %s seconds" % (time.time() - unique_nodes_time))
        else:
            open(newFile_path, "a").close()
            all_logger.warning("Recycler output and/or MetaPlasmidSPAdes output are empty. Could not align.")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

    return newFile_path

def MetaPlasmidOnly():
    ''' Getting list of candidates detected with metaplasmidSPAdes only'''
    toCompare = alignNodes()
    all_logger.info("Analyzing Recycler and MetaPlasmidSPAdes output")
    try:
        if not os.stat(toCompare).st_size == 0:
            df = pd.read_table(toCompare, header = None)
            custom_columns = ['qseqid', 'sseqid', 'evalue', 'length', 'pident',
                              'score', 'qcovs', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq']
            df.columns = custom_columns  # assigning specific columns
            df_filtered = df[(df['pident'] >= 99.0) & (df['qcovs'] >= 99.0)]  # filtering similar
            with open(output_dir + "/metaplasmidSPAdes/scaffolds.fasta") as metaPlasmid_file:
                metaPlasmid_list = []  # list for unique metaplasmids
                for record in SeqIO.parse(metaPlasmid_file, "fasta"):
                    metaPlasmid_list.append(record.id)
                RecordToAppend = df_filtered['sseqid'].values.tolist()
                metaPlasmid_dif = [rec for rec in metaPlasmid_list if
                                   rec not in RecordToAppend]  # writing unique plasmids
            all_logger.info("Unique MetaPlasmidSPAdes nodes defined")
        elif not os.path.isfile(toCompare):
            all_logger.warning("There is no alignment. Recycler or MetaPlasmidSPAdes output is empty.")
            with open(output_dir + "/metaplasmidSPAdes/scaffolds.fasta") as metaPlasmid_file:
                metaPlasmid_list = []  # list for unique metaplasmids
                for record in SeqIO.parse(metaPlasmid_file, "fasta"):
                    metaPlasmid_list.append(record.id)
                metaPlasmid_dif = [rec for rec in metaPlasmid_list]
        else:
            all_logger.warning("There is no alignment. All records in Recycler and MetaPlasmidSPAdes are unique.")
            with open(output_dir + "/metaplasmidSPAdes/scaffolds.fasta") as metaPlasmid_file:
                metaPlasmid_list = []  # list for unique metaplasmids
                for record in SeqIO.parse(metaPlasmid_file, "fasta"):
                    metaPlasmid_list.append(record.id)
                metaPlasmid_dif = [rec for rec in metaPlasmid_list]
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()
    return metaPlasmid_dif

def ReplaceMetaPlasmidID():
    metaplasmidSpades = output_dir + "/metaplasmidSPAdes/scaffolds.fasta"
    new_scaffold = output_dir + "/NewScaffolds.fasta"
    try:
        with open(new_scaffold, "wt") as new_spades:
            records = SeqIO.parse(metaplasmidSpades, 'fasta')
            for record in records:
                all_logger.info("We are replacing MetaPlasmidSPAdes record IDs./n Initial record ID is: %s" % record.id)
                prefix = str(uniqueNum) + "_N"
                record.description = record.description.replace("N", prefix)
                record.id = record.description
                all_logger.info("We are assigning new ID to MetaPlasmidSPAdes record./n New record ID is: %s" % record.id)
                SeqIO.write(record, new_spades, 'fasta')
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()
    return new_scaffold

def AppendToCombined():
    ''' Generating combined fasta for candidates'''
    recycler = output_dir + "/assembly_graph.cycs.fasta"
    try:
        with open(recycler) as recycler_fasta:
            combined_dir = all_out + "/CombinedOutput.fasta"
            combined_fasta = open(combined_dir, 'a')
            for header in recycler_fasta:
                prefix = ">" + str(uniqueNum) + "_"
                combined_fasta.write(header.replace(">", prefix))
        listToCompare = MetaPlasmidOnly()
        spades_records = SeqIO.parse(ReplaceMetaPlasmidID(), 'fasta')
        for record in spades_records:
            name_regex = re.search("([A-Z])\w+.+", record.id)
            list = []
            if name_regex is not None:
                name = name_regex.group(0)
                all_logger.info("Appending ' %s ' to list." % name)
                list.append(name)
            for name in list:
                if name in listToCompare:
                    all_logger.info("Writing ' %s ' to CombinedOutput." % name)
                    SeqIO.write(record, combined_fasta, 'fasta')
                else:
                    all_logger.warning("%s is duplicate" % name)
        combined_fasta.close()
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def MakeDB(name, dbtype, dbpath):
    '''Making database from fasta file for running BLAST'''
    makeblastdb = blast_dir + "makeblastdb"
    db_path = os.path.join(db_dir, name)
    in_file = name + ".fasta"
    if not os.path.exists(dbpath):
        dbpath = db_path
        os.mkdir(dbpath)
        all_logger.info("Directory '% s' created" % dbpath)
    db_file = dbpath + name
    in_path = dbpath + in_file
    try:
        makeDB_time = time.time()
        if os.path.isfile(in_path):
            if dbtype == "nucl":
                cline = NcbimakeblastdbCommandline(makeblastdb,
                                                   dbtype = dbtype,
                                                   input_file = in_path,
                                                   out = db_file)
                all_logger.info(cline)
                os.system(str(cline))
                stdout, stderr = cline()
            elif dbtype == "prot":
                cline = NcbimakeblastdbCommandline(makeblastdb,
                                                   dbtype = dbtype,
                                                   input_file = in_path,
                                                   out = db_file)
                all_logger.info(cline)
                os.system(str(cline))
                stdout, stderr = cline()
            else:
                all_logger.warning("******No database type provided******")
                error_tracker()
            all_logger.info("DataBase %s creation took %s seconds" % (name, (time.time() - makeDB_time)))
            return db_file
        else:
            all_logger.warning("No 'fasta' file provided. The Database was probably ready.")
            return db_file
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def RunBlast(name, dbpath):
    ''' Running BLAST'''
    query_file = all_out + "/CombinedOutput.fasta"
    newFile_path = all_out + "/" + name + ".csv"
    if not os.path.exists(dbpath):
        os.mkdir(dbpath)
    path = dbpath + "/"
    try:
        blast_time = time.time()
        if not os.path.isfile(newFile_path):
            if name == "CombinedOutput":
                copyfile(all_out+"/"+name+".fasta", path+"CombinedOutput.fasta")
                blast_type = "blastn"
                blast = blast_dir + blast_type
                path_to_DB = MakeDB(name, "nucl", path)
                all_vs_all = f'{res}/all_vs_all.csv'
                cline = NcbiblastnCommandline(blast,
                                              query=query_file,
                                              db=path_to_DB,
                                              evalue=0.001,
                                              out=all_vs_all,
                                              outfmt='6 qseqid sseqid evalue length pident score '
                                                     'qcovs qstart qend sstart send qseq sseq'
                                              )
                all_logger.info(cline)
                os.system(str(cline))
                stdout, stderr = cline()
                all_logger.info("BLAST running of %s took %s seconds" % (name, (time.time() - blast_time)))
                if os.path.exists(path):
                    shutil.rmtree(path)
            else:
                all_logger.warning("******Provide with valid blast command*****")
                all_logger.info("BLAST running of %s took %s seconds" % (name, (time.time()- blast_time)))
        else:
            all_logger.warning("Looks like BLAST has already got us '%s.csv'" % name)
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

if not os.path.isfile(lib_file):
    header = ['Sample', 'Library_size']
    with open(lib_file, 'a') as f_object:
        writer_object = writer(f_object)
        writer_object.writerow(header)
        f_object.close()

### Reading sample matrix file to get all read files and info
with open(sample_matrix) as adressFile:
    matrix = np.loadtxt(adressFile, dtype = "str")
    uniqueNum = None
    for i in matrix:
        FileNum = i[0]
        sampleNum = i[1]
        uniqueNum = str(sampleNum[-3:])
        all_logger.info("Assigning unique number: " + uniqueNum)
        forward = i[2]
        #forward = "ftp://" + forward
        reverse = i[3]
        #reverse = "ftp://" + reverse
        input_dir = CreateFolderDirectory('Input')
        fw_name = GetFileName(forward)
        rv_name = GetFileName(reverse)
        file_fw = GetFile(forward)
        file_rv = GetFile(reverse)
        lib_size = GetLibrary(input_dir)
        for_csv = [uniqueNum, lib_size]
        with open(lib_file, 'a') as f_object:
            writer_object = writer(f_object)
            writer_object.writerow(for_csv)
            f_object.close()
        output_dir = CreateFolderDirectory('Output')
        all_logger.info("***************WE START SPAdes for sample #%s****************" % uniqueNum)
        RunSpades()
        all_logger.info("***************SPAdes for sample #%s finished****************" % uniqueNum)
        all_logger.info("***************WE START RECYCLER for sample #%s****************" % uniqueNum)
        Recycler()
        all_logger.info("***************RECYCLER for sample #%s finished****************" % uniqueNum)
        all_logger.info("***************WE START MetaPlasmidSPAdes for sample #%s****************" % uniqueNum)
        MetaPlasmidSPAdes()
        all_logger.info("***************MetaPlasmidSPAdes finished for sample #%s****************" % uniqueNum)
        shutil.rmtree(input_dir)
        all_logger.info("Input directory removed for sample #%s" % uniqueNum )
        all_logger.info("***************APPENDING RECORDS for sample #%s TO COMBINED FILE****************" % uniqueNum)
        AppendToCombined()
        all_logger.info("***************FINISHED APPENDING RECORDS for sample #%s TO COMBINED FILE****************" % uniqueNum)

all_logger.info("Starting BLAST for Combined Output all vs all")

RunBlast("CombinedOutput", db_dir+"/CombinedOutput")

copyfile(f'{all_out}/CombinedOutput.fasta', f'{resource}/CombinedOutput.fasta')
all_logger.info("--- %s seconds ---" % (time.time() - start_time))
