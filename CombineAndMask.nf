#!/usr/bin/env nextflow

params.outdir = "Masking"
params.project_name = "Project"

process CombineFasta {

    conda 'seqkit'

    publishDir params.outdir, mode: 'copy', saveAs: { filename -> if (filename.endsWith(".fasta")) {"${projectName}_unmasked.fasta"}
                                                    else if (filename.endsWith(".csv")) {"${projectName}_sample_list.csv"} }

    input:
        val projectName
        path inputdir
        path fastas

    output:
        path "final.fasta", emit: unmasked_fasta
        path "sampleList.csv", emit: sampleList

    script:
    """
    cat ${inputdir}/*.fasta > all.fasta
    seqkit grep -n -f ${fastas} all.fasta -o final.fasta
    seqkit seq -n final.fasta | tr '\n' ',' > sampleList.csv
    """

}

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


workflow {

    inputdir_ch = Channel.fromPath(params.inputdir)
    fastas_ch = Channel.fromPath(params.sample_list)
    projectName_ch = Channel.of(params.project_name)

    CombineFasta(projectName_ch, inputdir_ch, fastas_ch)
    Masking(projectName_ch, CombineFasta.out.unmasked_fasta, CombineFasta.out.sampleList)

}
