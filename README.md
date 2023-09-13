# ![indrop](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) Indrop-Flow 

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Nextflow version](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/indrops.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/indrops/builds)


Indrops analysis pipeline at BioCore@CRG

The pipeline is based on the DropEST tool:
https://github.com/hms-dbmi/dropEst

Rewrite the .nf to make the pipeline compatible in DSL2 without needing to install BioNextflow. 

## Installing
1. install docker or singularity.
   - Pull docker images used in the workflow.
   - The final list of docker images used in this workflow:

   ```
   REPOSITORY                       TAG                  IMAGE ID       
   dropest_centos                   latest               590beb8096e7   
   quay.io/biocontainers/multiqc    1.15--pyhdfd78af_0   585f4076bcc7    
   quay.io/biocontainers/kallisto   0.46.2--h4f7b962_1   dec25607ea35  
   biocontainers/samtools           v1.9-4-deb_cv1       f210eb625ba6  
   biocontainers/fastqc             v0.11.5              fb9766783947   
   ```
   
2. git clone this repo; create dropest docker images (rebuilt because of compiling errors 20230830):

   ```
   cd ./docker
   docker build -t dropest_centos .
   ```

3. run nextflow using indrop_new.nf


## Running the pipeline
The parameters are listed when using ```nextflow run indrop_new.nf --help``` command.

Run with `nextflow run indrop_new.nf -resume --read_len 38`

```
N E X T F L O W  ~  version 23.04.3
Launching `indrop_new.nf` [determined_ampere] DSL2 - revision: fab418aaa4


╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╦┌┐┌┌┬┐┬─┐┌─┐┌─┐╔═╗╦  ╔═╗╦ ╦
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ║│││ ││├┬┘│ │├─┘╠╣ ║  ║ ║║║║
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ╩┘└┘─┴┘┴└─└─┘┴  ╚  ╩═╝╚═╝╚╩╝

====================================================
BIOCORE@CRG indropSEQ - N F  ~  version 2.0
====================================================
pairs                         : /your/path/indrop/data/test_{1,2}.fastq.gz
read_len                      : 38
genome                        : /your/path/indrop/anno/Homo_sapiens.GRCh38.cdna.all.fa.gz
annotation                    : /your/path/indrop/anno/gencode.v21.annotation.gtf
config                        : /your/path/indrop/conf/indrop_v1_2.xml
barcode_list                  : /your/path/indrop/conf/indrop_v1_2_barcodes.txt
email                         : yourmail@yourdomain
mtgenes                       :
dbdir                         : /your/path/indrop/db
version                       : 1_2
keepmulti                     : YES
library_tag                   : AGATATAA
output (output folder)        : output_v1_2
```

You can change them either by using the command line:
```
nextflow run indrop_new.nf --pairs "data/{1,2}.fastq.gz" --version 1_2 > log
```
or changing the params.file
You can use the nextflow options for sending the execution in with resuming a failed one (-resume).

```
nextflow run indrop_new.nf --pairs "data/{1,2}.fastq.gz" --version 1_2 -resume > log
```


## Indrop versions v1, v2 and v3 are supported
### Version 1 and 2
**Parameter version: "V1-2"**
* File 1: barcode reads. Structure:
  * Cell barcode, part 1
  * Spacer
  * Cell barcode, part 2
  * UMI
* File 2: gene reads

### Version 3
**Parameter version: "V3_3"**
* File 1: cell barcode
* File 2: 
  * cell barcode 
  * UMI 
* File 3: gene read

**Parameter version: "V3_4"**
* File 1: cell barcode
* File 2: 
  * cell barcode 
  * UMI 
* File 3: gene read
* File 4: library_tag

The parameter **library_tag** is only needed with version **V3_4**


# Parameters
1. Parameters are specified within the **params.config** file

## The pipeline
1. QC: Run FastQC on raw reads. It stores the results within **QC** folder.
1. Indexing: It makes the index of the genome by using Kallisto.
1. dropTag: It creates a "tagged" fastq file with information about the single cell that originated that read in the header. 
1. Alignment: It aligns tagged reads to the indexed genome by using Kallisto. Reasults are stored in **Alignments** folder.
1. dropEst: It provides the estimation of read counts per gene per single cell. The results are in **Estimated_counts** folder and consists of an R data object, a file with a list of cells (aka barcode combinations), another with a list of genes and a matrix in Matrix Market format (https://en.wikipedia.org/wiki/Matrix_Market_exchange_formats).
1. dropReport: It reads the R data oject produced by the dropEst step to produce a quality report. It needs a list of mitochondrial genes. 
1. multiQC: It wraps the QC from fastQC and Kallisto mapping in a single output. 

## Other notes
1. This pipeline is tested with the test dataset under ./data with indrop v1_2 setting on a personal laptop with 8GB RAM.
2. Included a jupyter notebook in `./notebooks` to explore, process, and visualize the gene-cell matrix generated by dropest.
   - Explore: Aggregate per cell, per gene; PCA/t-SNE/UMAP for dimention reduction.
   - Process: Using scanpy to QC for mitochondrial genes, standardize, and transform the matrix; Run tests to find out marker genes.
   - Visualize: Plot QC, cluster, marker gene heatmap figures.
   - Additional libraries required: seaborn, matplotlib, gget, sklearn, umap, scanpy

