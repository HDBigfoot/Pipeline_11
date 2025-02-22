#!/usr/bin/env nextflow

params.sample_name = "Sample"
params.outdir = "Results"
params.ref = "${projectDir}/Reference/NC_000962.3.fasta"
params.ref_index = "${projectDir}/Reference/NC_000962.3.fasta.fai"
params.ref_dict = "${projectDir}/Reference/NC_000962.3.dict"

log.info """
Pipeline_11
RC3ID & CentraBioRes
Universitas Padjadjaran
================================
sample     : $params.sample_name
reads      : $params.raw_read1 & $params.raw_read2
outdir     : $params.outdir
================================
"""

include { Trimming } from '${projectDir}/modules/Trimming.nf'
include { Mapping } from '${projectDir}/modules/Mapping.nf'
include { Dedup } from '${projectDir}/modules/Dedup.nf'
include { Calling } from '${projectDir}/modules/Calling.nf'
include { Filtering } from '${projectDir}/modules/Filtering.nf'

process fastaConversion {

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

workflow {

    ref_file = file(params.ref)
    ref_index_file = file(params.ref_index)
    ref_dict_file = file(params.ref_dict)

    sampleName_ch = Channel.of(params.sample_name)
    rawRead1_ch = Channel.fromPath(params.raw_read1)
    rawRead2_ch = Channel.fromPath(params.raw_read2)

    trimming(sampleName_ch, rawRead1_ch, rawRead2_ch)
    mapping(sampleName_ch, trimming.out.fastp_R1, trimming.out.fastp_R2, ref_file, ref_index_file, ref_dict_file)
    dedup(sampleName_ch, mapping.out.bwa_aligned, ref_file, ref_index_file, ref_dict_file)
    calling(sampleName_ch, dedup.out.bam_processed, ref_file, ref_index_file, ref_dict_file)
    filtering(sampleName_ch, calling.out.called_vcf, calling.out.called_idx, ref_file, ref_index_file, ref_dict_file)
    fastaConversion(sampleName_ch, filtering.out.clean_vcf, filtering.out.clean_idx, ref_file, ref_index_file, ref_dict_file)

}
