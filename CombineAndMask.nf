#!/usr/bin/env nextflow

params.outdir = "Masking"
params.project_name = "Project"

process CombineFasta {

    conda 'seqkit'

    publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${projectName}_unmasked.fasta" }

    input:
        val projectName
        path inputdir
        path fastas

    output:
        path "final.fasta", emit: unmasked_fasta

    script:
    """
    cat ${inputdir}/*.fasta > all.fasta
    seqkit grep -n -f ${fastas} all.fasta -o final.fasta
    """

}

workflow {

    inputdir_ch = Channel.fromPath(params.inputdir)
    fastas_ch = Channel.fromPath(params.sample_list)
    sampleList_ch = Channel.fromPath(params.sample_list).splitText().toList()
    projectName_ch = Channel.of(params.project_name)

    CombineFasta(projectName_ch, inputdir_ch, fastas_ch)

}
