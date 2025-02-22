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

    /*
     *publishDir params.outdir + "/Calling", mode: 'copy', saveAs: {filename -> if (filename.endsWith(".vcf")) {"${sampleName}_raw.snps.indels.vcf"}
     *                                                             else if (filename.endsWith(".idx")) {"${sampleName}_raw.snps.indels.vcf.idx"}}
     *
     */

    input:
        val sampleName
        path bam_processed
        path ref
        path ref_index
        path ref_dict

    output:
        path "${bam_processed}_raw.snps.indels.vcf", emit: called_vcf
        path "${bam_processed}_raw.snps.indels.vcf.idx", emit: called_idx

    script:
    """
    gatk HaplotypeCaller --sample-ploidy 1 --R ${ref} --I ${bam_processed} --O ${bam_processed}_raw.snps.indels.vcf --max-assembly-region-size 600 --standard-min-confidence-threshold-for-calling 30.0 --min-assembly-region-size 300
    """

}

process filtering {

    conda 'gatk4'

    publishDir params.outdir + "/VCF", mode: 'copy', saveAs: {filename -> if (filename.endsWith("_clean.snps.vcf")) {"${sampleName}.vcf"}
                                                             else if (filename.endsWith("_clean.snps.vcf.idx")) {"${sampleName}.vcf.idx"}}

    input:
        val sampleName
        path called_vcf
        path called_idx
        path ref
        path ref_index
        path ref_dict

    output:
        path "${called_vcf}_clean.snps.vcf", emit: clean_vcf
        path "${called_vcf}_clean.snps.vcf.idx", emit: clean_idx

    script:
    """
    gatk SelectVariants --R ${ref} --V ${called_vcf} --select-type-to-include SNP --O ${called_vcf}_raw.snps.vcf
    gatk VariantFiltration --R ${ref} --V ${called_vcf}_raw.snps.vcf --filter-expression "QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 12" --filter-name "FAILED" --O ${called_vcf}_flagged.snps.vcf
    gatk SelectVariants --R ${ref} --V ${called_vcf}_flagged.snps.vcf --exclude-filtered --O ${called_vcf}_clean.snps.vcf
    """

}

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
