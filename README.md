# Strand-Seq pooling library identification from SNP genotyping (1kgp samples)

This repository aims to develop a pipeline to automatically reidentify cells sequenced in a pooling Strand-Seq experiments (limited to 1kgp samples).

## Rationale

1kgp is the single large-scale genomic project where genotype information is available at the sample level (GT FORMAT columns in the VCF file). Data is available [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 

## Aim

As in the HGSVC large-scale inversion project, Strand-Seq libraries are pooled on a 96-well plate, it is not possible to know directly which sample is corresponding to a specific cell. Using the 1kgp genotype information will allow us to identify cells origin.

## Code 

Pipeline was created using [snakemake](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) and relies essentially on bcftools, freebayes, samtools and custom python scripts.

## Troubleshooting

Main issue at the moment is the generation of the Directed Acyclic Graph (DAG) of execution, as each cell needs to be compared to 3202 samples. For a 96-well plate, the number of jobs to be done will be around 308k including one-to-one comparison and final analysis.

## TODO

- [ ] Split in 2 snakefiles: one complete & one to be bind to mosaicatcher pooling version
- [ ] Simplify to reduce DAG time creation (lower AF, isec on a single file?)