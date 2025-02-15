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

    publishDir params.outdir + "/Aligned", mode: 'copy', saveAs: { filename -> "${sampleName}_aligned.sam"}

    input:
        val sampleName
        path fastp_R1
        path fastp_R2
        path ref
        path ref_index
        path ref_dict

    output:
        path "${fastp_R1}_aligned.sam", emit: bwa_aligned

    script:
    """
    bwa index ${ref}
    bwa mem -M ${ref} ${fastp_R1} ${fastp_R2} > ${fastp_R1}_aligned.sam
    """

}

process dedup {

    conda 'gatk4'

    publishDir params.outdir + "/Dedup", mode: 'copy', saveAs: { filename -> if (filename.endsWith("_metrics.txt")) {"${sampleName}_metrics.txt"}
                                                               else if (filename.endsWith("_dedup.bam")) {"${sampleName}_dedup.bam"}}

    input:
        val sampleName
        path bwa_aligned
        path ref
        path ref_index
        path ref_dict

    output:
        path "${bwa_aligned}_metrics.txt"
        path "${bwa_aligned}_dedup.bam", emit: bam_processed

    script:
    """
    gatk FixMateInformation --ASSUME_SORTED false --I ${bwa_aligned} --O ${bwa_aligned}_fixed.bam
    gatk SortSam --I ${bwa_aligned}_fixed.bam --O ${bwa_aligned}_sorted.bam --SO coordinate
    gatk AddOrReplaceReadGroups --I ${bwa_aligned}_sorted.bam --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM ${sampleName} --O ${bwa_aligned}_rg.bam
    gatk MarkDuplicates --REMOVE_DUPLICATES true --CREATE_INDEX true --ASSUME_SORTED true --I ${bwa_aligned}_rg.bam --M ${bwa_aligned}_metrics.txt --O ${bwa_aligned}_dedup.bam
    """

}

process calling {

    conda 'gatk4'

    publisdhDir params.outdir + "/Calling", mode: 'copy', saveAs: {filename -> if (filename.endsWith(".vcf")) {"${sample_name}_raw.snps.indels.vcf"}
                                                                  else if (filename.endsWith(".idx")) {"${sampleName}_raw.snps.indels.vcf.idx"}}

    input:
        val sampleName
        path bam_processed
        path ref
        path ref_index
        path ref_dict

    output:
        path "${bam_processed}_raw.snps.indels.vcf", emit: vcf
        path "${bam_processed}_raw.snps.indels.vcf.idx", emit: idx

    script:
    """
    gatk HaplotypeCaller --sample-ploidy 1 --R ${ref} --I ${bam_processed} --O ${bam_processed}_raw.snps.indels.vcf --max-assembly-region-size 600 --standerd-min-confidence-threshold-for-calling 30.0 --min-assembly-region-size 300
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
    calling(sampleName_ch, dedup.out.bam_processed, ref_file, ref_index_file, ref_dict_file

}
