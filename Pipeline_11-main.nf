#!/usr/bin/env nextflow

params.outdir = "Results"
params.ref = "${ProjectDir}/Reference/NC_000962.3.fasta"

process trimming {

    conda '${projectDir}/main_conda.yml'

    publishDir params.ourdir + "/Trimming", mode: 'copy'

    input:
        val sampleName
        path rawRead1
        path rawRead2

    output:
        path "${sampleName}_R1_fastp.fastq.gz", emit: fastp_R1
        path "${sampleName}_R2_fastp.fastq.gz", emit: fastp_R2
        path "${sampleName}.fastp.html"

    script:
    """
    fastp -i ${rawRead1} -I ${rawRead2} -o ${sampleName}_R1_fastp.fastq.gz -O ${sampleName}_R2_fastp.fastq.gz --length_required 50 --html ${sampleName}.fastp.html
    """

}

process mapping {

    conda '${projectDir}/main_conda.yml'

    publishDir params.outdir + "/Aligned", mode: 'copy'

    input:
        val sampleName
        path fastp_R1
        path fastp_R2
        path ref

    output:
        path "${sampleName}_Aligned.sam", emit: bwa_Aligned

    script:
    """
    bwa mem -M ${ref} ${fastp_R1} ${fastp_R2} > ${sampleName}_Aligned.sam
    """

}

workflow {

    ref_file = file(params.ref)

    sampleName_ch = Channel.fromPath(params.sample-name)
    rawRead1_ch = Channel.fromPath(params.raw-read1)
    rawRead2_ch = Channel.fromPath(params.raw-read2)

    trimming(sampleName_ch, rawRead1_ch, rawRead2_ch)
    mapping(sampleName_ch, trimming.out.fastp_R1, trimming.out.fastp_R2, ref_file)

}
