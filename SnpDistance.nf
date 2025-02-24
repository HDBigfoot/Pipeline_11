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

include { SnpDistanceAnalysis } from './modules/SnpDistanceAnalysis.nf'
include { RemoveDuplicates } from './modules/RemoveDuplicates.nf'

workflow {

    projectName_ch = Channel.of(params.project_name)
    inputFasta_ch = Channel.fromPath(params.input_fasta)

    SnpDistanceAnalysis(projectName_ch, inputFasta_ch)
    RemoveDuplicates(projectName_ch, SnpDistanceAnalysis.out.dupedSnpDistance)

}
