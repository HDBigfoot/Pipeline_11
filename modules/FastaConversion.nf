#!/usr/bin/env nextflow

process FastaConversion {

    conda 'gatk4'

    publishDir params.outdir + "/FASTA", mode: 'copy', saveAs: { filename -> "${sampleName}.fasta"}

    input:
        val sampleName
        path clean_vcf
        path clean_idx
        path ref
        path ref_index
        path ref_dict

    output:
        path "${clean_vcf}_clean.fasta"

    script:
    """
    gatk FastaAlternateReferenceMaker --R ${ref} --V ${clean_vcf} --O ${clean_vcf}_raw.fasta
    sed 's/1 NC_000962.3:1-4411532/'${sampleName}'/' ${clean_vcf}_raw.fasta > ${clean_vcf}_clean.fasta
    """

}
