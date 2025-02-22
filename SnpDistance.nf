#!/usr/bin/env nextflow

params.outdir = "SNP_Distance"
params.project_name = "Project"

log.info """
Pipeline_11
SNP Distance Analysis
RC3ID & CentraBioRes
Universitas Padjadjaran
================================
project    : $params.project_name
outdir     : $params.outdir
================================
"""

include { SnpDistanceAnalysis } from '${projectDir}/modules/SnpDistanceAnalysis.nf'

process RemoveDuplicates {

    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${projectName}_SNP_Distance.csv"}

    input:
        val projectName
        path snpDistance

    output:
        path "${snpDistance}.csv"

    script:
    """
    Rscript ${projectDir}/Scripts/Remove_Duplicates.R ${snpDistance} > ${snpDistance}.csv
    """

}

workflow {

    projectName_ch = Channel.of(params.project_name)
    inputFasta_ch = Channel.fromPath(params.input_fasta)

    SnpDistanceAnalysis(projectName_ch, inputFasta_ch)
    RemoveDuplicates(projectName_ch, SnpDistanceAnalysis.out.dupedSnpDistance)

}
