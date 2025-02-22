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

process Masking {

    publishDir params.outdir, mode: 'copy', saveAs: { filename -> if (filename.endsWith(".bed")) {"${projectName}_to_be_masked.bed"}
                                                    else if (filename.endsWith(".masked.fasta")) {"${projectName}.masked.fasta"}}

    conda 'bedtools'

    input:
        val projectName
        path unmaskedFasta
        path sampleList

    output:
        path "${unmaskedFasta}_to_be_masked.bed"
        path "${unmaskedFasta}.masked.fasta"

    script:
    """
    python ${projectDir}/Scripts/prepare_multifasta_bedfile.py ${sampleList} > ${unmaskedFasta}_to_be_masked.bed
    bedtools maskfasta -fi ${unmaskedFasta} -bed ${unmaskedFasta}_to_be_masked.bed -fo ${unmaskedFasta}.masked.fasta
    """

}


workflow {

    inputdir_ch = Channel.fromPath(params.inputdir)
    fastas_ch = Channel.fromPath(params.sample_list)
    projectName_ch = Channel.of(params.project_name)

    CombineFasta(projectName_ch, inputdir_ch, fastas_ch)
    Masking(projectName_ch, CombineFasta.out.unmasked_fasta, CombineFasta.out.sampleList)

}
