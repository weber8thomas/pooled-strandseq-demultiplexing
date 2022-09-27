import subprocess
subprocess.Popen("module load Java/1.8.0_221", shell=True)


import hail as hl
hl.init(min_block_size=256)
import os
os.environ['PYSPARK_SUBMIT_ARGS'] = "--driver-memory 32G pyspark-shell"

import pandas as pd

# gs = "/g/korbel/weber/MosaiCatcher_files/EXTERNAL_DATA/snv_sites_to_genotype/1000G_SNV_with_GT/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz"
# gs = "/data/scratch/gnomAD/v3/gnomad.genomes.v3.1.sites.ht/"

# mt = hl.import_vcf(gs, force_bgz=True, reference_genome="GRCh38", call_fields=["GT"])
mt = hl.import_vcf(snakemake.input.onekgp, force_bgz=True, reference_genome="GRCh38", call_fields=["GT"])
mt = mt.select_cols()


# first_gt = mt.filter_cols(mt.s == "HG00096")
first_gt = mt.filter_cols(mt.s == snakemake.wildcards.sample)
first_gt = first_gt.filter_rows(hl.is_snp(first_gt.alleles[0], first_gt.alleles[1]))

first_gt = first_gt.filter_rows(hl.agg.any(first_gt.GT.is_het()))


first_gt = first_gt.annotate_rows(
    info= first_gt.info.annotate(
        ID= hl.str(first_gt.locus) + "_" + hl.str(first_gt.alleles[0]) + "_" + hl.str(first_gt.alleles[1]),
        AC= first_gt.info.AC[0],
        AF= first_gt.info.AF[0],
    )

)

first_gt = first_gt.select_rows(first_gt.info.ID, first_gt.info.AC, first_gt.info.AF)

first_gt = first_gt.filter_rows(hl.agg.any(first_gt.AF <= 0.01))

first_gt.rows().to_pandas().to_csv(snakemake.output[0], compression='gzip', index=False, sep="\t")
# first_gt.rows().to_pandas().to_csv("test.tsv.gz", compression='gzip', index=False, sep="\t")

