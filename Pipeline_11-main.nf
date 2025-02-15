#!/usr/bin/env nextflow

params.sample_name = "Sample"
params.outdir = "Results"
params.ref = "${projectDir}/Reference/NC_000962.3.fasta"
params.ref_index = "${projectDir}/Reference/NC_000962.3.fasta.fai"
params.ref_dict = "${projectDir}/Reference/NC_000962.3.fasta.dict"

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

process trimming {

    conda 'fastp'

    publishDir params.outdir + "/Trimming", mode: 'copy', saveAs: { filename -> if (filename.endsWith("1_fastp.fastq.gz")) {"${sampleName}_1_fastp.fastq.gz"}
                                                                  else if (filename.endsWith("2_fastp.fastq.gz")) {"${sampleName}_2_fastp.fastq.gz"}
                                                                  else if (filename.endsWith("html")) {"${sampleName}.fastp.html"}}

    input:
        val sampleName
        path rawRead1
        path rawRead2

    output:
        path "${rawRead1}_fastp.fastq.gz", emit: fastp_R1
        path "${rawRead2}_fastp.fastq.gz", emit: fastp_R2
        path "${rawRead1}.fastp.html"

    script:
    """
    fastp -i ${rawRead1} -I ${rawRead2} -o ${rawRead1}_fastp.fastq.gz -O ${rawRead2}_fastp.fastq.gz --length_required 50 --html ${rawRead1}.fastp.html
    """

}

process mapping {

    conda 'bwa'

    publishDir params.outdir + "/Aligned", mode: 'copy', saveAs: { filename -> "${sampleName}_Aligned.sam"}

    input:
        val sampleName
        path fastp_R1
        path fastp_R2
        path ref
        path ref_index
        path ref_dict

    output:
        path "${fastp_R1}_Aligned.sam", emit: bwa_Aligned

    script:
    """
    bwa index ${ref}
    bwa mem -M ${ref} ${fastp_R1} ${fastp_R2} > ${fastp_R1}_Aligned.sam
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

}
