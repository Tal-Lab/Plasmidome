import time, re, os, sys, subprocess, wget, logging, shutil
import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio import SearchIO
from Bio import SeqIO
from zipfile import ZipFile
from Bio import bgzf
import tarfile
import gzip
from shutil import copyfile

start_time = time.time()

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

parent_dir = r"../"
combined = sys.argv[1]
bash_file_Ref = r"./BowtieReference.csh"
bash_file_Map = r"./BowtieMapper.csh"
cov_bash = r"./length_cov.csh"
sample_matrix = r"../res/samples_matrix.txt"
logs = r'../logs'
mycwd = os.getcwd()

def setup_logger(name, log_file, level=logging.DEBUG):
    """To setup as many loggers as you want"""
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger

def error_tracker():
    error_file = logs + "/Error_Bowtie2.log"
    error_logger = setup_logger("error_logging", error_file)
    exc_type, exc_value, traceback = sys.exc_info()
    error_logger.exception(sys.exc_info())
    all_logger.exception(sys.exc_info())
    print("exception caught: %s" % exc_type, exc_value)

def CreateDirectory (name):
    path = os.path.join(parent_dir, name)
    if not os.path.exists(path):
        os.mkdir(path)
        #all_logger.info("Directory '%s' created" % path)
    return path

all_out = CreateDirectory('Output')

all_data_log = logs + "/Bowtie2.log"

all_logger = setup_logger('second_logger', all_data_log)

def CreateFolderDirectory (name):
    parent_dir = CreateDirectory(name)
    FolderName = FileNum + name
    pathFolder = os.path.join(parent_dir, FolderName)
    if not os.path.exists(pathFolder):
        os.mkdir(pathFolder)
        all_logger.info("Directory '%s' created" % pathFolder)
    return pathFolder

def Zip_bio(name):
    os.chdir(output_dir)
    # create a ZipFile object
    compress_time = time.time()
    zip_file_name = name + uniqueNum + '.bgzf'
    zipObj = bgzf.BgzfWriter(zip_file_name, 'a')
    # Add multiple files to the zip
    all_logger.info("Compressing old BAM/SAM files.")
    [zipObj.write(f) for f in os.listdir(output_dir) if re.match(r'^\d+', f)]
    all_logger.info("Compressing old BAM/SAM files took %s seconds" % (time.time() - compress_time))
    zipObj.close()
    # Remove uncompressed files
    [os.remove(f) for f in os.listdir(output_dir) if re.match(r'^\d+', f)]
    all_logger.warning("Files before compression removed.")
    os.chdir(mycwd)

def GetFileName (link):
    file_regex = re.search("\/S+\w+.fastq.gz", link)
    file_name = file_regex.group(0)[1:]
    #print(file_name)
    #all_logger.info("Getting link %s" % file_name)
    return file_name

def GetFile (input_dir,link):
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

def Build_Reference():
    ref_folder = os.path.join(all_out, "reference")
    if not os.path.exists(ref_folder):
        os.mkdir(ref_folder)
    ref_name = "plasmids_double"
    path_to_file = f'{ref_folder}/{ref_name}'
    try:
        ref_time = time.time()
        if not os.path.exists(path_to_file):
            all_logger.info("Starting preparations for Bowtie2. Convertion of Doubled OUTPUT fasta to reference genome.")
            subprocess.call(["bash", bash_file_Ref, combined, ref_name, ref_folder])
            all_logger.warning("Check ERRORs file for records")
            all_logger.info("reference files created")
            all_logger.info("Convertion to Genome Reference took %s seconds" % (time.time() - ref_time))
            assert os.path.isfile(bash_file_Ref), "bash file missing"
        else:
            all_logger.warning("Looks like Bowtie2-Build has already created reference genome. Check ERRORs file for records")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()
    return path_to_file

def Bowtie2():
    outname = uniqueNum + "_double.sam"
    out_dir = f'{output_dir}/{outname}'
    try:
        bowtie_time = time.time()
        if not os.path.isfile(out_dir):
            input_dir = CreateFolderDirectory('Input')
            file_fw = GetFile(input_dir, forward)
            file_rv = GetFile(input_dir, reverse)
            subprocess.call(["bash", bash_file_Map, file_fw, file_rv, output_dir, Build_Reference(), out_dir, uniqueNum, combined])
            all_logger.info("Bowtie2 mapping took %s seconds" % (time.time()-bowtie_time))
            all_logger.warning("Something may have gone wrong here. Check file output directory for '.sam' file.")
            # all_logger.info("***************Cleaning out fastq-reads input folder******************")
            # shutil.rmtree(input_dir)
        else:
            all_logger.warning("Looks like Bowtie2 has already mapped reads. Check file 'assembly_graph.cycs.fasta'.")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

mapping_analysis = os.path.join(all_out, "Mapping")
if not os.path.exists(mapping_analysis):
    os.mkdir(mapping_analysis)
    all_logger.info("Directory '% s' created" % mapping_analysis)

def Samtools_Run():
    in_name = uniqueNum + "_double.sorted.bam"
    in_dir = f'{output_dir}/{in_name}'
    out_covn = uniqueNum + "_cov.csv"
    out_covd = f'{mapping_analysis}/{out_covn}'
    out_depn = uniqueNum + "_depth.csv"
    out_depd = f'{mapping_analysis}/{out_depn}'
    try:
        bowtie_time = time.time()
        if not os.path.isfile(out_covd) and not os.path.isfile(out_depd):
            subprocess.call(
                ["bash", cov_bash, combined, in_dir, out_covd])
            all_logger.info("Coverage and depth calculation took %s seconds" % (time.time() - bowtie_time))
            all_logger.warning("Something may have gone wrong here. Check file output directory for '.sam' file.")
        else:
            all_logger.warning("Looks like coverage was already obtained.")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def MakeDF(uniqueNum):
    name = uniqueNum+"_cov.csv"
    print(name)
    cov_file = f'{mapping_analysis}/{name}'
    df = pd.read_csv(cov_file, sep='\t')
    print(df)
    df.sort_values(by=['#rname'])
    df_cov = df[["#rname", "coverage"]]
    new_name = uniqueNum + "_coverage"
    df_cov = df_cov.rename(columns={'#rname': 'rname',
        'coverage': new_name})
    print(df_cov)
    outfile = f'{parent_dir}Output/{"all_cov.csv"}'
    output = open(outfile, 'w')
    if uniqueNum == "258" and os.stat(outfile).st_size == 0:
        df_cov.to_csv(output, index=False)
        output.close()
    return df_cov

def Merger(num):
    if num == "258" and not os.path.isfile(f'{output_dir}/{"all_cov.csv"}'):
        MakeDF(num)
    else:
        df_prim = pd.read_csv(f'{output_dir}/{"all_cov.csv"}')
        df_to_join = MakeDF(num)
        df_joined = pd.merge(df_prim, df_to_join, on='rname')
        outfile = f'{output_dir}/{"all_cov.csv"}'
        os.remove(outfile)
        output = open(outfile, 'w')
        df_joined.to_csv(output, index=False)
        output.close()

with open(sample_matrix) as adressFile:
    matrix = np.loadtxt(adressFile, dtype = "str")
    uniqueNum = None
    for i in matrix:
        FileNum = i[0]
        sampleNum = i[1]
        uniqueNum = str(sampleNum[-3:])
        all_logger.info("Assigning unique number: " + uniqueNum)
        forward = i[2]
        reverse = i[3]
        output_dir = CreateFolderDirectory('Output')
        Zip_bio('single')
        all_logger.info("*************** WE START ALIGNMENT OF FASTQ READS TO COMBINED OUTPUT #%s ****************" % uniqueNum)
        Bowtie2()
        all_logger.info("*************** ALIGNMENT for sample #%s finished ****************" % uniqueNum)
        all_logger.info("***************WE START EXTRACTING ALIGNMENT COVERAGE AND DEPTH FOR SAMPLE #%s****************" % uniqueNum)
        Samtools_Run()
        all_logger.info("*************** EXTRACTING INFO for sample #%s finished ****************" % uniqueNum)
        Merger(uniqueNum)
        all_logger.info("*************** Writing INFO for sample #%s into general all.cov file finished ****************" % uniqueNum)
        Zip_bio('double')

all_logger.info("--- %s seconds ---" % (time.time() - start_time))
