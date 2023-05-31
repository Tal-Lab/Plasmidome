# Plasmidome detection and analysis

## Table of Contents

- [Project Description](#project-description)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Project Description

To extend the repertoire of environmental marine plasmids, we established a pipeline for the de novo assembly of plasmids in the marine environment by analyzing available microbiome metagenomic sequencing data. By applying the pipeline to data from the Red Sea, we identified 362 plasmid candidates. We showed that the distribution of plasmids corresponds to environmental conditions, particularly, depth, temperature and physical location. At least seven of the 362 candidates are most probably real plasmids, based on a functional analysis of their open reading frames (ORFs). Only one of the seven has been described previously. Three plasmids were identified in other public marine metagenomic data from different locations all over the world; these plasmids contained different cassettes of functional genes at each location. Analysis of antibiotic- and metal-resistance genes revealed that the same positions that were enriched with genes encoding resistance to antibiotics were also enriched with resistance to metals, suggesting that plasmids contribute site-dependent phenotypic modules to their ecological niches. Finally, half of the ORFs (50.8%) could not be assigned to a function, emphasizing the untapped potential of the unique marine plasmids to provide proteins with multiple novel functions. Paper under review.

## Requirements

- BLAST+ 
- <a href="https://github.com/ablab/spades">SPAdes 3.14+</a>
- bwa, samtools, bowtie2, seqkit, bcftools
- <a href="https://github.com/Shamir-Lab/Recycler">Recycler</a>
- <a href="https://github.com/hyattpd/Prodigal">Prodigal</a>
- <a href="http://eggnog-mapper.embl.de/">eggNOG-mapper</a>
- <a href="https://github.com/ebi-pf-team/interproscan">interproscan</a>
- <a href="https://github.com/santirdnd/COPLA">COPLA</a> 
- plasmid dataset in FASTA format (in this study the plasmid database <a href="https://ccb-microbe.cs.uni-saarland.de/plsdb">PLSDB</a>)
- File with sampling stations description and data, called stations.txt: 
sample_id, Station number, station location, temperature, depth, any other physical parameters 
- File with samples raw reads links, called samples_matrix.txt: 
station_id (stationNumber_stationDepth), sample_id, forward link, reverse link
- File with databases dbLinks.txt
- ACLAMEproteins.zip upload to DBs folder manually

## Installation

Clone the source code from GitHub:

```
git clone https://github.com/Tal-Lab/Plasmidome && ./Plasmidome/setup.sh
```

Make sure you set correct directories in the <b>.env</b> file to BLAST, Recycler, SPAdes. 

Make sure you have bwa, samtools, bowtie2, seqkit, bcftools in you PATH. 

## Usage

1. To perform steps I-III of the Plasmidome detection Pipeline run python script **PlasmidomeAssembly.py**.
It should generate among all files **all_vs_all.csv** and **CombinedOutput.fasta** in the resource (/res) folder.

2. To perfrom steps IV – V run python script **Pipeline_Filtering.py**. It should generate file **Plasmids_byClusters.csv** in data_calculations folder. Manually go over it. Choose 1 representative of each cluster, where possible. Make a **new_names.csv** file in your working directory with two columns: "old_name" and "new_name". Where old_name is a representative of the cluster, and new_name is a new distinct name for the cluster. In cases, where there’s only one plasmid in the cluster, you may leave the old name. 
Run python script **clusters_to_fasta.py** to generate **filtered_plasmids.fasta** in resource (/res) folder. 

3. To perform step VI of the pipeline, run **BowtieRunner.csh**. You will have file **plasmids_double.fasta** in Output folder.


4. Run **blast_plsdb.csh** to compare detected candidates with plasmids in PLSDB. This script will generate file **plsdb.csv** in Output folder. Run **db_statistics.py** to detect known plasmids - the result will be prented out in the console.
5. 
6. To get predicted ORFs first run **prodigal_runner.csh**, which will generate **plasmids_double_proteins.faa** in Output folder. Next run <b>orf_filtering.py</b> to obtain ORFs only single plasmid length. The filtered ORFs will be written into **Filtered_ORFs.fasta**.
7. 

## Feedback

If you have any feedback, please reach out to us at lucyandrosyuk@gmail.com or leave a comment in Discussions.

Please, include “PLASMIDOME github” in the subject and specify your issue.

## Acknowledgements

This study was supported (in part) by grant no. 3-17700 from the Office of the Chief Scientist, Israel Ministry of Health. L.A. is the recipient of a Hi-Tech, Bio-Tech, and Chemo-tech fellowship of Ben-Gurion University of the Negev.