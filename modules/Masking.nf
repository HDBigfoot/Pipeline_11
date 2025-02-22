#!/usr/bin/env nextflow

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
