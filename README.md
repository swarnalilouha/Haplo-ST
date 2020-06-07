# Haploid Sequence-Typer (Haplo-ST)
===================================
Haplo-ST is a pipeline for generating standardized whole-genome multi-locus sequence typing (wgMLST) of *Listeria monocytogenes (Lm)* isolates from whole-genome sequencing (WGS) data. Along with wgMLST profiles, this pipeline also generates assembled allele sequences and identifies paralogous genes for each *Lm* isolate.

Haplo-ST works in several steps. It takes in (i) raw WGS reads as input, (ii) cleans the raw data according to user-specified parameters, (iii) assembles genes across loci by mapping to genes from reference strains, (iv) assigns allelic profiles to assembled genes and provides a wgMLST profile for each isolate. Data standardization and exchangeability relies on the tool assigning allelic profiles based on the centralized nomenclature defined by the BIGSdb-*Lm* database. More broadly, Haplo-ST is flexible and can be adapted to provide allele-based subtyping for any haploid organism, with the installation of an organism-specific gene database with associated allelic nomenclature within the Virtual Machine (VM) on which this pipeline runs.

## Requirements:
Haplo-ST should work on the Linux environment within the VM available at https://www.dropbox.com/s/zh10mpk0lukhadn/Haplo-ST.ova?dl=0. This VM can be imported within any local installation of Oracle VM VirtualBox.

Haplo-ST uses FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) for cleaning sequencing reads, YASRA (https://github.com/aakrosh/YASRA) for assembling alleles by mapping to reference genes, and BIGSdb-Lm (https://bigsdb.pasteur.fr/listeria/) for assigning allelic profiles to assembled alleles.

## Starting Inputs:
1. raw_data : This folder contains the forward and reverse reads for each sample in fastq format.

   Read names can contain a combination of alphabets and/or integers and should be in the format:

     Forward-reads_R1_001.fastq : Forward reads in fastq format

     Reverse-reads_R2_001.fastq : Reverse reads in fastq format

2.	genes : This folder contains the reference genes used for subtyping. Each reference gene should be included as a separate file in fasta format. The choice and number of reference genes used for subtyping depends on the user; fewer reference genes can be used for low resolution typing, whereas higher resolution can be achieved by increasing the number of genes used in the analysis.

3.	header : tab-delimited file containing the names of reference genes used for subtyping.

      Please keep in mind that the first line of this file starts with two tabs followed by reference gene names separated 
      by tabs

4.	param.txt : file containing parameters specified by user for cleaning raw fastq files.

    Raw reads are cleaned with the following tools from the FASTX-Toolkit in a sequential manner:
   
      • fastx_trimmer : Trims low quality bases from the beginning or end of reads.
      
      • fastq_quality_trimmer : Trims nucleotides with lower quality from the end of the sequence and purges reads lower
        than a specified length after trimming.
          
      • fastq_quality_filter : Filters and retains reads which contain a minimum percent of bases with a desired quality
        score.

5.	Makefile : file used by YASRA for assembling alleles.

## Usage:
Place all the required inputs and the script 'pipeline.pl' in the working directory and run:

    perl pipeline.pl param.txt

The program asks for the value of the ‘length criteria’: 

    Enter the value for the ‘length criteria’

Length criteria = length of assembled allele / length of reference gene used for assembly
This is a measure of the minimum length of the assembled allele that will be retained by the pipeline and assigned an allelic designation. All alleles having a length of less than the ‘length criteria’ will be filtered out by Haplo-ST. Insert the value of this parameter and press Enter.

If everything runs smoothly, the file ‘header’ will be populated with wgMLST profiles of Lm isolates, for which WGS data had been provided in “raw_data”. This file can be opened up in excel for better readability. Each row in this file corresponds to the allele-based profile of a specific Lm isolate, whose name is written in the first column of the row. Integers in every cell correspond to the allele ID of the associated gene present in the column header. New alleles not having a defined allele ID in BIGSdb-Lm are designated as close matches to another allele of the corresponding gene. For example, ‘closest match: 114’ refers to a new allele which closely matches allele 114 of the corresponding gene. Haplo-ST also generates the following folders:

   • consensus : Each file within this folder contains alleles assembled for a sample in fasta format. Files are named 
     according to names of samples provided to the pipeline.
     
   • paralog_warnings : Each file within this folder contains a list of paralogs identified within a sample.
   
## Test Dataset
A sample toy dataset can be found in the 'test-dataset' sub-directory.

•	The ‘raw_data’ folder contains paired reads for two isolates of Lm. These raw files should be unzipped before use.

•	The value of the parameters in the ‘param.txt’ can be changed to suit your needs for a particular project.

•	The ‘header’ file can be modified to include less or more reference genes used for subtyping.

•	The ‘genes’ folder contains fasta files for 2554 reference genes listed in the BIGSdb-Lm dataset.

•  The ‘Makefile’ remains unchanged for all analysis. 

    


