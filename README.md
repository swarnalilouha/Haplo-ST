# Haploid Sequence-Typer (Haplo-ST)
Haplo-ST is a pipeline for generating standardized whole-genome multi-locus sequence typing (wgMLST) of Listeria monocytogenes (Lm) isolates from whole-genome sequencing (WGS) data. Along with wgMLST profiles, this pipeline also generates assembled allele sequences and identifies paralogous genes for each Lm isolate.

Haplo-ST works in several steps. It takes in (i) raw WGS reads as input, (ii) cleans the raw data according to user-specified parameters, (iii) assembles genes across loci by mapping to genes from reference strains, (iv) assigns allelic profiles to assembled genes and provides a wgMLST profile for each isolate. Data standardization and exchangeability relies on the tool assigning allelic profiles based on the centralized nomenclature defined by the BIGSdb-Lm database. More broadly, Haplo-ST is flexible and can be adapted to provide allele-based subtyping for any haploid organism, with the installation of an organism-specific gene database with associated allelic nomenclature within the Virtual Machine (VM) on which this pipeline runs.

#Requirements:

Haplo-ST should work on the Linux environment within the VM available at https://www.dropbox.com/s/zh10mpk0lukhadn/Haplo-ST.ova?dl=0. This VM can be imported within any local installation of Oracle VM VirtualBox.

Haplo-ST uses FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) for cleaning sequencing reads, YASRA (https://github.com/aakrosh/YASRA) for assembling alleles by mapping to reference genes, and BIGSdb-Lm (https://bigsdb.pasteur.fr/listeria/) for assigning allelic profiles to assembled alleles.
