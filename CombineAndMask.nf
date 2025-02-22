#!/usr/bin/env nextflow

params.outdir = "Masking"
params.project_name = "Project"

log.info """
Pipeline_11
Combine and Mask Fastas
RC3ID & CentraBioRes
Universitas Padjadjaran
================================
project    : $params.project_name
outdir     : $params.outdir
================================
"""

include { CombineFasta } from '${projectDir}/modules/CombineFasta.nf'
include { Masking } from '${projectDir}/modules/Masking.nf'


workflow {

    inputdir_ch = Channel.fromPath(params.inputdir)
    fastas_ch = Channel.fromPath(params.sample_list)
    projectName_ch = Channel.of(params.project_name)

    CombineFasta(projectName_ch, inputdir_ch, fastas_ch)
    Masking(projectName_ch, CombineFasta.out.unmasked_fasta, CombineFasta.out.sampleList)

}
