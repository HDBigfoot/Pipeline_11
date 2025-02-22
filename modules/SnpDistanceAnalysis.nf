#!/usr/bin/env nextflow

process SnpDistanceAnalysis {

    conda 'snp-dists'

    //publishDir params.outdir, mode: 'copy', saveAs: { filename -> "${projectName}_SNP_Distance.csv" }

    input:
        val projectName
        path maskedFasta

    output:
        path "${maskedFasta}_snpDistance.csv", emit: dupedSnpDistance

    script:
    """
    snp-dists -mc ${maskedFasta} > ${maskedFasta}_snpDistance.csv
    """

}
