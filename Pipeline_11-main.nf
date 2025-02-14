#!/usr/bin/env nextflow

params.outdir = "Results"

process trimming {

    conda '${projectDir}/main_conda.yml'

    publishDir params.ourdir + "/Trimming", mode: 'copy'

    input:
        val sampleName
        path inputdir

    output:
        path "${sampleName}_R1_fastp.fastq.gz"
        path "${sampleName}_R2_fastp.fastq.gz"
        path "${sampleName}.fastp.html"

    script:
    """
    fastp -i ${inputdir}/${sampleName}_1.fastq.gz -I ${inputdir}/${sampleName}_2.fastq.gz -o ${sampleName}_R1_fastp.fastq.gz -O ${sampleName}_R2_fastp.fastq.gz --length_required 50 --html ${sampleName}.fastp.html
    """

}

workflow {

    sampleName_ch = Channel.fromPath(params.sample-names).splitCsv().view{ "After splitCsv: $it" }.flatten().view{ "After flatten: $it"}
    inputdir_ch = Channel.fromPath(params.inputdir)

    trimming(sampleName_ch, inputdir_ch)

}
