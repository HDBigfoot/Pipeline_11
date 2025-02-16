#!/usr/bin/env nextflow

params.outdir = "Masking"

process CombineFasta {

    conda 'seqkit'

    publishDir params.outdir, mode: 'copy'

    input:
        path inputdir
        path fastas

    output:
        path "final.fasta"

    script:
    """
    cat ${inputdir}/*.fasta > all.fasta
    seqkit grep -n -f ${fastas} all.fasta -o final.fasta
    """

}

workflow {

    inputdir_ch = Channel.fromPath(params.inputdir)
    fastas_ch = Channel.fromPath(params.sample_list)

    CombineFasta(inputdir_ch, fastas_ch)

}
