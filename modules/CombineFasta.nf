#!/usr/bin/env nextflow

process CombineFasta {

    conda 'seqkit'

    publishDir params.outdir, mode: 'copy', saveAs: { filename -> if (filename.endsWith(".fasta")) {"${projectName}_unmasked.fasta"}
                                                    else if (filename.endsWith(".csv")) {"${projectName}_sample_list.csv"} }

    input:
        val projectName
        path inputdir
        path fastas

    output:
        path "final.fasta", emit: unmasked_fasta
        path "sampleList.csv", emit: sampleList

    script:
    """
    cat ${inputdir}/*.fasta > all.fasta
    seqkit grep -n -f ${fastas} all.fasta -o final.fasta
    seqkit seq -n final.fasta | tr '\n' ',' > sampleList.csv
    """

}
