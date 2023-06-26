"""
Code created:
July 13, 2021

Code author:
Lucy Androsiuk
"""

import time,logging,re,os,sys,wget,gzip,tarfile,zipfile,shutil, glob, dotenv_setup
import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline
import dotenv_setup

# Directories
work_dir = r"../"
db_dir = r"../DBs"
resource = r"../res"
logs = r"../logs"
out_dir =f'{work_dir}/Output'
db_links = f"{resource}/dbLinks.txt"
blast_dir = os.getenv("BLAST")
prot_file = f"{resource}/Filtered_ORFs.fasta"

ACLAMEproteins = r"../DBs/ACLAMEproteins.zip"

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
    error_file = logs + "/Error_P5.log"
    error_logger = setup_logger("error_logging", error_file)
    exc_type, exc_value, traceback = sys.exc_info()
    error_logger.exception(sys.exc_info())
    all_logger.exception(sys.exc_info())
    print("exception caught: %s" % exc_type, exc_value)

all_data_log = logs + "/data_P5.log"
all_logger = setup_logger('second_logger', all_data_log)

def UnzipDB(name, dirPath):
    extract_dir1 = db_dir + "/" + name + "/"
    extract_dir2 = db_dir + "/"
    try:
        if dirPath.endswith('.zip'):
            zip_file = zipfile.ZipFile(dirPath, 'r')
            all_logger.info("Extracting file")
            zip_file.extractall(extract_dir1)
            zip_file.close()
        elif dirPath.endswith("tar.gz"):
            tar = tarfile.open(dirPath, "r:gz")
            all_logger.info("Extracting file")
            tar.extractall(extract_dir1)
            tar.close()
        elif dirPath.endswith("tar.bz2"):
            tar = tarfile.open(dirPath, "r:bz2")
            all_logger.info("Extracting file")
            tar.extractall(extract_dir2 + "ABresDB/")
            tar.close()
        elif dirPath.endswith("tar"):
            tar = tarfile.open(dirPath, "r:")
            all_logger.info("Extracting file")
            tar.extractall(extract_dir2 + "ABresDB/")
            tar.close()
        elif dirPath.endswith("fasta.gz"):
            file_name = extract_dir2 + "BacMetDB/" + name + ".fasta"
            file_to_write = open(file_name, "wb")
            with gzip.open(dirPath, "rb") as f:
                bindata = f.read()
            all_logger.info("Extracting file")
            file_to_write.write(bindata)
            file_to_write.close()
        elif dirPath.endswith("fas.gz"):
            file_name = extract_dir2 + "Virulence/" + name + ".fasta"
            file_to_write = open(file_name, "wb")
            with gzip.open(dirPath, "rb") as f:
                bindata = f.read()
            all_logger.info("Extracting file")
            file_to_write.write(bindata)
            file_to_write.close()
        else:
            all_logger.warning("*****Could not extract file, as now extractor found******")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def MakeDB(name, dbtype, dbpath):
    makeblastdb = blast_dir + "makeblastdb"
    db_path = os.path.join(db_dir, name)
    in_file = name + ".fasta"
    if not os.path.exists(dbpath):
        dbpath = db_path
        os.mkdir(db_path)
        all_logger.info("Directory '% s' created" % dbpath)
    db_file = f'{dbpath}/{name}'
    in_path = f'{dbpath}/{in_file}'
    try:
        makeDB_time = time.time()
        if not any(fname.endswith('.pdb') for fname in os.listdir(dbpath)):
            if os.path.isfile(in_path):
                if dbtype == "prot":
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
        else:
            all_logger.warning("The database was ready")
            return db_file

    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def RunBlast(in_file,name, dbtype, DBpath, out_file):
    try:
        if not os.path.isfile(out_file) or os.stat(out_file).st_size == 0:
            blast_time = time.time()
            blast = blast_dir + 'blastp'
            path_to_DB = MakeDB(name, dbtype, DBpath)
            cline = NcbiblastpCommandline(blast, query = in_file,
                                      db = path_to_DB,
                                      evalue = 0.001,
                                      out = out_file,
                                      outfmt = '6 qseqid sseqid stitle evalue length pident mismatch score '
                                               'qcovs qstart qend sstart send qseq sseq'
                                      )
            all_logger.info(cline)
            os.system(str(cline))
            stdout, stderr = cline()
            all_logger.info("BLAST running of %s took %s seconds" % (name, (time.time() - blast_time)))
        else:
            all_logger.warning("------- Blast of this protein file:\n %s \n in %s database was already performed --------" % (in_file,name))
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def GetDB():
    try:
        with open(db_links) as db_adress:
            matrix = np.loadtxt(db_adress, dtype="str")
            matrix_del = np.delete(matrix, 0, 0)
            for line in matrix_del:
                all_logger.info("Getting DataBase from %s" % line)
                dirDB_name = line[0]
                DB_link = line[1]
                name = line[2]
                dbtype = line[3]
                if dbtype == 'nucl':
                    continue
                DBpath = os.path.join(db_dir, dirDB_name)
                if not os.path.exists(DBpath):
                    os.mkdir(DBpath)
                    print("Directory '%s' created" % DBpath)
                all_logger.info("DataBase name is %s" % name)
                #DB_name = "/" + name
                DBFile = os.path.join(DBpath, name)

                file_name = os.path.basename(DBFile)
                index_of_dot = file_name.index('.')
                new_name = file_name[:index_of_dot]
                if not os.path.exists(DBFile):
                    all_logger.info("Downloading Database %s by link %s", DBFile, DB_link)
                    wget.download(DB_link, DBFile)
                if zipfile.is_zipfile(DBFile) or tarfile.is_tarfile(DBFile) or DBFile.endswith(".gz"):
                    all_logger.info("Extracting DataBase")
                    UnzipDB(new_name, DBFile)
                if dirDB_name == "ABresDB":
                    all_logger.info("Removing all non-protein databases from AntiBioticResistanceDB")
                    for file in os.listdir(DBpath):
                        filename_re = re.search("^protein.*.fasta$", file)
                        if filename_re:
                            ab_file = os.path.basename(file)
                            index_dot = ab_file.index('.')
                            ab_name = ab_file[:index_dot]
                            all_logger.info("Starting BLAST for Predicted genes vs %s" % ab_name)
                            RunBlast(prot_file, ab_name, dbtype, DBpath, f'{out_dir}/{ab_name}.csv')
                        else:
                            all_logger.info(filename_re)
                            all_logger.warning("File %s is not protein. Removing" % file)
                            try:
                                all_logger.info("TRY delete %s" % file)
                                os.remove(DBpath + "/" + file)
                            except:
                                all_logger.error("Cant delete file %s" % file)
                else:
                    all_logger.info("Starting BLAST for Predicted genes vs %s" % new_name)
                    RunBlast(prot_file, new_name, dbtype, DBpath, f'{out_dir}/{new_name}.csv')
        all_logger.info("Extracting ACLAME DataBase")
        UnzipDB("ACLAMEproteins", ACLAMEproteins)
        all_logger.info("Starting BLAST for Corrected Output vs ACLAMEproteins")
        RunBlast(prot_file, "ACLAMEproteins", 'prot', db_dir + "/ACLAMEproteins", f'{out_dir}/ACLAMEproteins.csv')

    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

all_logger.info("************************ Starting BLASTp searches for Predicted genes *************************")
RunBlast(prot_file, "all_proteins", "prot", db_dir + "/all_proteins", f'{out_dir}/all_proteins.csv')
GetDB()
all_logger.info("--- %s seconds ---" % (time.time() - start_time))
