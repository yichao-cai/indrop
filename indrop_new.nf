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
read_len                      : ${params.read_len}
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
// configFile          = params.config
// barcodeFile         = params.barcode_list
// if (params.mtgenes != "") mitocgenesFile  = params.mtgenes
// db_folder		    = params.dbdir
params.dropestScript       = "$projectDir/bin/dropestr/dropReport.Rsc"

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
    tag "$sample_id"
    publishDir params.outputQC, mode : 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val("$sample_id"), path("FASTQC_$sample_id")

    script:
    """
    mkdir FASTQC_${sample_id}
    fastqc -o FASTQC_${sample_id} -f fastq -q ${reads}
    """
}

process GET_READ_LEN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastqcDir)

    output:
    tuple val("$sample_id"), stdout

    script: // Get length from read2, as read1 is barcode.
    """
    unzip -o -q ${fastqcDir}/${sample_id}_2_fastqc.zip -d ./${fastqcDir}
    get_read_len.py --workdir ${fastqcDir} --prefix ${sample_id}_2
    """
}


process DROPTAG {
    // label 'indrop_multi_cpus'
    tag "$pair_id"
    publishDir params.filt_folder, mode:'copy'
    
    input:
    tuple val(pair_id), path(reads)
    path configFile
    
    output:
    tuple val(pair_id), path("*.tagged.fastq.gz") 
    tuple val(pair_id), path("*.tagged.params.gz")
    tuple val(pair_id), path("*.tagged.rds")

    script:
    def v3_params = ""
    if (params.library_tag != "") {
    	v3_params = "-t ${params.library_tag}"
    }
    """
    droptag -r 0 -S -s ${v3_params} -p ${task.cpus} -c ${configFile} ${reads}
    """
}

process KALLISTO_INDEX {
    label 'kallisto'
    storeDir params.dbdir

    input:
    path genome

    output:
    path "Kallisto_hg38_transcripts.idx"
    
    script:
    """
    kallisto index -i Kallisto_hg38_transcripts.idx ${genome}
    """
}

process KALLISTO_QUANT {
    label 'kallisto'
    tag "$pair_id"
    publishDir params.outputMapping, mode:'copy'
    // errorStrategy { task.exitStatus=1 ? 'ignore' : 'terminate' }    // Ignore exit error code 1 in kallisto quant, which means no reads mapped. But publishDir not triggered.

    input:
    tuple val(pair_id), val(reads)
    path genome_index

    output:
    tuple val(pair_id), path("kallisto_${pair_id}")
    tuple val(pair_id), path("kallisto_${pair_id}/pseudoalignments.bam")
    
    shell:
    """
    #!/bin/bash 
    kallisto quant --pseudobam --single -i ${genome_index} -o kallisto_${pair_id} -l 100 -s 10 ${reads}
    if [ \$? -eq 1 ]; then
        echo "Exit error code 1 catched. Cause: no reads mapped in Kallisto quant."
        echo "Ignoring..."
        exit 0
    fi
    """
}

process REMOVE_MULTI {
    // label 'big_mem_cpus'
    tag "$pair_id"
    publishDir params.outputMapping, mode:'copy'

    input:
    tuple val(pair_id), path(aln)

    output:
    tuple val(pair_id), path ("${pair_id}_univoc_s.bam")

    script:
    """
    samtools view -H ${aln} > ${pair_id}_univoc_s.sam
    samtools view -@ ${task.cpus} ${aln} | grep \"\\<NH:i:1\\>\" >> ${pair_id}_univoc_s.sam
    samtools view -@ ${task.cpus} -Sb ${pair_id}_univoc_s.sam > ${pair_id}_univoc_s.bam 
    rm ${pair_id}_univoc_s.sam
    """
}

process DROPEST {
    // label 'indrop_one_cpu'
    tag "$pair_id"
    publishDir params.est_folder, mode:'copy'

    input:
    tuple val(pair_id), path(tags), path(params_est)
    path barcodeFile
    path annotationFile
    path configFile

    output:
    tuple val(pair_id), path ("*.rds")
    tuple val(pair_id), path ("*.tsv")
    tuple val(pair_id), path ("*.mtx")

    script:     
    """
    mv ${barcodeFile} barcode_file.txt
    dropest -m -w -P -r ${params_est} -c ${configFile} -o ${pair_id}.est ${tags} 
    """
}

process DROP_REPORT {
    // label 'dropreport'
    tag "$pair_id"
    errorStrategy = 'ignore'    // Ignore error. Possible cases: 1) empty alignment generated a empty .rds.
    publishDir params.rep_folder, mode: 'copy'

    input:
    tuple val(pair_id), path(estimate), path(droptag) 
    path dropestScript
    
    output:
    tuple val(pair_id), path("${pair_id}_report.html")

    script:
    def mitopar = ""
    def mitocmd = ""
    if (params.mtgenes != "") {
        mitopar = " -m mitoc.rds" 
        mitocmd = "gene_to_rds.r ${mitocgenesFile} mitoc.rds"
    }
    """
    ${mitocmd}
    echo Rscript --vanilla ${dropestScript} -t ${droptag} -o ${pair_id}_report.html ${mitopar} ${estimate} 
    Rscript --vanilla ${dropestScript} -t ${droptag} -o ${pair_id}_report.html ${mitopar} ${estimate} 
    """
}

process MULTIQC {
    tag "$pair_id"
    publishDir params.outputMultiQC, mode:'copy'
    

    input:
    tuple val(pair_id), path(raw_fastqc_files), path(kallisto_out_files) 

    output:
    path "multiqc_report_${pair_id}.html"

    script:
    """
    multiqc .
    mv multiqc_report.html multiqc_report_${pair_id}.html 
    """
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
    fastqc_ch = FASTQC(read_pairs)
    // Get read length from fastq if not speficied in params. Needed for STAR.
    if ( params.read_len == "" ) {
        readLen_ch = GET_READ_LEN(fastqc_ch)
    } else {
        readLen_ch = read_pairs
            .map{ tuple( it[0], params.read_len) }
    }
    // DropTag
    (tagged_reads_ch, tagged_params_ch, tagged_rds_ch) = DROPTAG(read_pairs, params.config)
    // Kallisto
    ref_ch = KALLISTO_INDEX(params.genome)
    (kallisto_out, bam_ch) = KALLISTO_QUANT(tagged_reads_ch, ref_ch)
    // Remove multi-mapping alignment
    if (params.keepmulti == "NO" ){
        bam_rmdup_ch = REMOVE_MULTI(bam_ch)
    } else {
        bam_rmdup_ch = bam_ch
    }
    // DropEst
    est_in_ch = bam_rmdup_ch.join(tagged_params_ch)
    (est_rds_ch, est_tsv_ch, est_mtx_ch) = DROPEST(est_in_ch, params.barcode_list, params.annotation, params.config)
    // Report
    report_in_ch = est_rds_ch.join(tagged_rds_ch)
    report_ch = DROP_REPORT(report_in_ch, params.dropestScript)
    // MultiQC
    fastqc_ch
        .join(kallisto_out)
        | MULTIQC
}

// Actions after completion
workflow.onComplete {
    println "Pipeline BIOCORE_indrop_nf_new completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
* send mail
*/
workflow.onComplete {
    def subject = 'indrop_new execution'
    def recipient = "${params.email}"
    // def attachment = "${params.outputMultiQC}/multiqc_report.html"
    def attachment = Channel
        .fromPath( "${params.outputMultiQC}/*.html" )
        .collect()

    ['mail', '-s', subject, '-a', attachment, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}