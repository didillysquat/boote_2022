#!/usr/bin/env nextflow

/*
Here we will run preprocessing of the transcript sequencing for the Boote project
There are two species of coral: Acropora intermedia and Montipora digitata

For A. intermedia there is a reference genome here and an associated gff that has all of the august predictions
We will first extract the transcript CDS and map to these. We will then mapp the transciptome mappings back to parent gene in DESEQ
To extract the transcript isoform CDS we used gffread from the dockerhub image zavolab/gffread:0.11.7
Command was: gffread aint.gff -x acropora.transcript.isoforms.CDS.fasta -g aint.fasta
It gave errors because of the transcription_end_site features not having IDs.
We then made kallisto indices using the jennylsmith/kallistov45.0:latest dockerhub
command: kallisto index -i acropora.transcript.isoforms.CDS.fasta.kallisto.index acropora.transcript.isoforms.CDS.fasta

For M. digitata we will use the CDS seqeunces of the transcriptome produced by Gonzalez-Pech
Paper: https://www.frontiersin.org/articles/10.3389/fmars.2017.00403/full
Transcriptome CDS fasta: https://github.com/PalMuc/Montipora_digitata_resources/tree/master/Gonzalez-Pech_et_al/Annotated_transcriptome
This fasta contains isoform CDS. To get to Trinity component (which is the closest analogue for gene) you have to compine the first part of the name e.g. TRXXXX with the first
two components of the second part of the name cX_gX. The isoform is given by the ix. We will map to the isoform CDS
and then consolidate to the Trinity component (gene) level for the DeSeq analysis.
Kallisto index command: kallisto index -i Mdigitata.annotated.fa.kallisto.index Mdigitata.annotated.fa

For preprocessing we will use fastp to adapter and low quality trim before pseudo mapping using kallisto
*/

acro_samples_ch = Channel.fromFilePairs("/home/humebc/projects/boote/raw_seq_data/Acropora/*_{1,2}.fq.gz").map{["Acropora", it[0], it[1]]}
monti_samples_ch = Channel.fromFilePairs("/home/humebc/projects/boote/raw_seq_data/Montipora/*_{1,2}.fq.gz").map{["Montipora", it[0], it[1]]}

acropora_kallisto_index = file("/home/humebc/projects/boote/references/acropora.transcript.isoforms.CDS.fasta.kallisto.index")
montipora_kallisto_index = file("/home/humebc/projects/boote/references/Mdigitata.annotated.fa.kallisto.index")

// Run fastp on all samples
process fastp{
    tag "${sample_name}"
    container "biocontainers/fastp:v0.20.1_cv1"
    publishDir "/home/humebc/projects/boote/fastp/${species}", pattern: "*.html", mode: "copy"

    input:
    tuple val(species), val(sample_name), path(reads) from acro_samples_ch.mix(monti_samples_ch)

    output:
    file "${sample_name}.fastp.html" into fastp_out_ch
    tuple val(species), val(sample_name), file("${sample_name}_1.clean.fq.gz"), file("${sample_name}_2.clean.fq.gz") into kallisto_in_ch

    script:
    """
    fastp -q 20 -i ${reads[0]} -I ${reads[1]} -o ${sample_name}_1.clean.fq.gz -O "${sample_name}_2.clean.fq.gz"
    mv fastp.html ${sample_name}.fastp.html
    """
}

process kallisto{
    tag "${sample}"
    container "jennylsmith/kallistov45.0:latest"
    publishDir "/home/humebc/projects/boote/kallisto_quant_results/${species}/${sample}", mode: "copy"
    cpus 10

    input:
    tuple val(species), val(sample), path(read_1_clean), path(read_2_clean) from kallisto_in_ch
    file acropora_kallisto_index
    file montipora_kallisto_index

    output:
    file "*" into kallisto_out_ch

    script:
    if (species == "Acropora")
    """
    kallisto quant -i $acropora_kallisto_index -o . -t ${task.cpus} $read_1_clean $read_2_clean
    """
    else if (species == "Montipora")
    """
    kallisto quant -i $montipora_kallisto_index -o . -t ${task.cpus} $read_1_clean $read_2_clean
    """
}
