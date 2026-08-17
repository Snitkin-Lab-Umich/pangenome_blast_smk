# pangenome_blast_smk
A snakemake pipeline for comparing Panaroo pangenome results against a BLAST search.

The goal of these scripts is to run a custom BLAST search for a provided list of accessory genes against a dataset of assemblies.
This pipeline is best used on the output files produced by Panaroo.
To start, modify the fields in config/config.yaml to match your dataset. You will need a list of accessory genes, directories of .fasta and .gff files, and a gene family presence/absence matrix in the same format used by Panaroo.
More detailed instructions are provided directly in the above config file and workflow/pangenome_blast.smk.
