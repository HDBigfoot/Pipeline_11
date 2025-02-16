This is a bioinformatics pipeline I wrote for the analysis of Mycobacterium tuberculosis genomes for the Research Center for Care and Control of Infectious Diseases, Universitas Padjadjaran.

This pipeline is modified from Mycodentifier (https://github.com/JordyCoolen/MyCodentifier.git)

Schildkraut JA, Coolen JPM, Severin H, Koenraad E, Aalders N, Melchers WJG, Hoefsloot W, Wertheim HFL, van Ingen J. MGIT Enriched Shotgun Metagenomics for Routine Identification of Nontuberculous Mycobacteria: a Route to Personalized Health Care. J Clin Microbiol. 2023 Mar 23;61(3):e0131822. doi: 10.1128/jcm.01318-22. Epub 2023 Feb 22. PMID: 36840602; PMCID: PMC10035320.

## Installation

Clone this repository:

```bash
git clone https://github.com/HDBigfoot/Pipeline_11.git
```

## Usage

Running main pipeline:

```bash
nextflow run /PATH/TO/PROJECT/Pipeline_11/Pipeline_11-main.nf --raw_read1 /PATH/TO/RAW/READS/<sample_name>_1.fastq.gz --raw_read2 /PATH/TO/RAW/READS/<sample_name>.fastq.gz --sample_name <sample_name>
```

Combining and masking FASTAs:

```bash
nextflow run /PATH/TO/PROJECT/Pipeline_11/CombineAndMask.nf --sample_list <list-of-samples>.txt --inputdir /PATH/TO/FASTA/FILES/ --project_name <project_name>
```

Running SNP-Distance Analysis:

```bash
nextflow run /PATH/TO/PROJECT/Pipeline_11/SnpDistance.nf --input_fasta /PATH/TO/COMBINED/FASTA/<combined_fasta>.masked.fasta --project_name <project_name>
```
