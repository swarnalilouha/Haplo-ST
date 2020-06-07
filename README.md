# Haploid Sequence-Typer (Haplo-ST)
Haplo-ST is a pipeline for generating standardized whole-genome multi-locus sequence typing (wgMLST) of *Listeria monocytogenes (Lm)* isolates from whole-genome sequencing (WGS) data. Along with wgMLST profiles, this pipeline also generates assembled allele sequences and identifies paralogous genes for each *Lm* isolate.

Haplo-ST works in several steps. It takes in (i) raw WGS reads as input, (ii) cleans the raw data according to user-specified parameters, (iii) assembles genes across loci by mapping to genes from reference strains, (iv) assigns allelic profiles to assembled genes and provides a wgMLST profile for each isolate. Data standardization and exchangeability relies on the tool assigning allelic profiles based on the centralized nomenclature defined by the BIGSdb-*Lm* database. More broadly, Haplo-ST is flexible and can be adapted to provide allele-based subtyping for any haploid organism, with the installation of an organism-specific gene database with associated allelic nomenclature within the Virtual Machine (VM) on which this pipeline runs.

## Requirements:

Haplo-ST should work on the Linux environment within the VM available at https://www.dropbox.com/s/zh10mpk0lukhadn/Haplo-ST.ova?dl=0. This VM can be imported within any local installation of Oracle VM VirtualBox.

Haplo-ST uses FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) for cleaning sequencing reads, YASRA (https://github.com/aakrosh/YASRA) for assembling alleles by mapping to reference genes, and BIGSdb-Lm (https://bigsdb.pasteur.fr/listeria/) for assigning allelic profiles to assembled alleles.

## Starting Inputs:

1. raw_data : This folder contains the forward and reverse reads for each sample in fastq format.

Read names can contain a combination of alphabets and/or integers and should follow the format:

Forward-reads_R1_001.fastq : Forward reads in fastq format

Reverse-reads_R2_001.fastq : Reverse reads in fastq format

2.	genes : This folder contains the reference genes used for subtyping. Each reference gene should be included as a separate file in fasta format. The choice and number of reference genes used for subtyping depends on the user; fewer reference genes can be used for low resolution typing, whereas higher resolution can be achieved by increasing the number of genes used in the analysis.

3.	header : tab-delimited file containing the names of reference genes used for subtyping.

Please keep in mind that the first line of this file starts with two tabs followed by reference gene names separated by tabs

4.	param.txt : file containing parameters specified by user for cleaning raw fastq files.

Raw reads are cleaned with the following tools from the FASTX-Toolkit in a sequential manner:

  •	fastx_trimmer :  Trims low quality bases from the beginning or end of reads.

  •	fastq_quality_trimmer  :  Trims nucleotides with lower quality from the end of the sequence and purges reads lower than     a specified length after trimming.

  •	fastq_quality_filter  :  Filters and retains reads which contain a minimum percent of bases with a desired quality s         score.

5.	Makefile : file used by YASRA for assembling alleles.

## Usage:

Simply run:


