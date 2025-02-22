#!/usr/bin/env nextflow

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
