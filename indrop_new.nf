#!/usr/bin/env nextflow

// Rewrite workflow in DSL2 syntax
// Test fastq files from https://www.refine.bio/experiments/SRP053052/droplet-barcoding-for-single-cell-transcriptomics-applied-to-embryonic-stem-cells

/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '2.0'

params.help            = false
params.resume          = false

log.info """

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╦┌┐┌┌┬┐┬─┐┌─┐┌─┐╔═╗╦  ╔═╗╦ ╦
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ║│││ ││├┬┘│ │├─┘╠╣ ║  ║ ║║║║
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ╩┘└┘─┴┘┴└─└─┘┴  ╚  ╩═╝╚═╝╚╩╝
                                                                                       
====================================================
BIOCORE@CRG indropSEQ - N F  ~  version ${version}
====================================================
pairs                         : ${params.pairs}
genome                        : ${params.genome}
annotation                    : ${params.annotation}
config                        : ${params.config}
barcode_list                  : ${params.barcode_list}
email                         : ${params.email}
mtgenes                       : ${params.mtgenes}
dbdir                         : ${params.dbdir}
version                       : ${params.version}
keepmulti                     : ${params.keepmulti}
library_tag                   : ${params.library_tag}
output (output folder)        : ${params.output}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// Input parameters
// genomeFile          = params.genome
// annotationFile      = params.annotation
configFile          = params.config
// barcodeFile         = params.barcode_list
// if (params.mtgenes != "") mitocgenesFile  = params.mtgenes
db_folder		    = params.dbdir
params.dropestScript       = "$projectDir/docker/dropestr/dropReport.Rsc"

// Output folders
params.outputfolder    = "${params.output}"
params.outputQC        = "${params.outputfolder}/QC"
params.outputMultiQC   = "${params.outputfolder}/multiQC"
params.outputMapping   = "${params.outputfolder}/Alignments"
params.filt_folder     = "${params.outputfolder}/Tagged_reads"
params.est_folder      = "${params.outputfolder}/Estimated_counts"
params.rep_folder      = "${params.outputfolder}/Reports"

/*
* Check essential files and emit them as value channel
*/
genome_ch = Channel.fromPath(params.genome, checkIfExists: true).collect()
annotation_ch = Channel.fromPath(params.annotation, checkIfExists: true).collect()
barcode_ch = Channel.fromPath(params.barcode_list, checkIfExists: true).collect()
if (params.mtgenes != "") {
    mitocgenes_ch  = Channel.fromPath(params.mtgenes, checkIfExists: true).collect()
}

if (params.version != "1_2" && params.version != "3_3" && params.version != "3_4") 
	exit 1, "Please define a valid version! It can be 1_2, 3_3, 3_4.\nRespectively version 1 or 2, version 3 with 3 files and version 4 with 4 files."

if (params.keepmulti != "NO" && params.keepmulti != "YES") 
	exit 1, "Please define a valid keepmulti value! It can YES or NO"


/*
*   Process
*/

process FASTQC {
    publishDir params.outputQC, mode : 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "FASTQC_$sample_id"

    script:
    """
    mkdir FASTQC_$sample_id
    fastqc -o FASTQC_$sample_id -f fastq -q ${reads}
    """
}

process DROPTAG {
    
}

/*
*   Workflow
*/

workflow {  
    /*
    * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
    * three elements: the pair ID, the first read-pair file and the second read-pair file
    */
    read_pairs = Channel
        .fromFilePairs( params.pairs, checkIfExists: true, size: (params.version == "1_2") ? 2 : (params.version == "3_3") ? 3 : 4)
        .ifEmpty { error "Cannot find any read pairs matching: ${params.pairs}" }  
    read_pairs.view()  

    // FASTQC
    FASTQC(read_pairs)
}