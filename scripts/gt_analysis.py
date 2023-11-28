import multiprocessing as mp
import parmap
import sys, os
import pandas as pd
import gzip
import numpy as np


## List of inputs
# ref = pd.read_csv(snakemake.input.ref, compression="gzip", sep="\t")
ref = pd.read_csv(
    "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CONCAT_TAB/LanexHGSVCpool2NEW/merge.txt.gz",
    compression="gzip",
    sep="\t",
)


# Prepare multiprocessing instances
m = mp.Manager()
# l_df_gb = m.list()
# l_df = m.list()
l_df_gb = list()
l_df = list()


# def mp_vcf(vcf, l_df_gb, l_df):
def mp_vcf(vcf):
    # Nb of rows to skip at the beginning of the VCF
    skip = len(
        [
            e.decode("ISO-8859–1").split("\n")
            for e in gzip.open("{}".format(vcf, "rb"))
            if e.decode("ISO-8859–1").split("\n")[0].startswith("##")
        ]
    )
    # Sample name
    sample = os.path.basename(vcf).replace(".vcf.gz", "")

    # Read dataframe
    df = pd.read_csv(vcf, compression="gzip", skiprows=skip, sep="\t")

    # Process CHROM column & add ID column
    df["#CHROM"] = df["#CHROM"].astype(str)
    df["#CHROM"] = df["#CHROM"].str.replace("chr", "")
    df["ID"] = (
        df["#CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"]
        + ":"
        + df["ALT"]
    )

    # Merge ref & df
    merge_df = pd.merge(ref, df, on="ID", how="right")

    merge_df["SAMPLE"] = merge_df["SAMPLE"].fillna("NAN")
    merge_df["QUERY_SAMPLE_CELL"] = sample

    # Groupby to count SNPs number / sample, then sort values
    merge_df_gb = (
        merge_df.groupby("SAMPLE")["ID"]
        .nunique()
        .sort_values(ascending=False)
        .reset_index()
    )
    merge_df_gb["QUERY_SAMPLE_CELL"] = sample
    # Compute the rank
    merge_df_gb["Rank"] = list(range(1, 1 + merge_df_gb.shape[0]))
    # Compute total nb of SNPs
    merge_df_gb["Total_SNP"] = df.shape[0]
    print(merge_df_gb)

    l_df_gb.append(merge_df_gb)
    l_df.append(merge_df)


# Process list of cells aplhabetically
# f_ldir = sorted(list(snakemake.input.vcf))
f_ldir = sorted(
    list(
        [
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A68.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A67.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G27.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C55.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E82.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A69.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G56.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E45.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C09.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A51.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A66.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E51.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C48.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E37.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E40.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G66.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G40.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E35.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E52.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A87.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A43.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C69.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A46.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A45.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C96.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C61.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C60.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E23.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E71.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C45.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C89.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G30.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G11.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G65.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E19.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A50.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G85.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G89.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A86.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G82.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C59.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C11.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A24.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A30.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E68.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E91.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A05.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E86.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E20.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C64.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A75.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C50.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C49.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C80.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E70.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C06.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G54.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G08.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A04.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C39.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C93.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E15.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G53.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G01.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E74.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G05.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A47.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E67.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C29.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E26.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G69.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C37.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C23.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A09.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E88.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A01.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G42.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A91.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A73.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A35.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G18.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G31.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A36.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G13.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G16.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A07.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A55.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G84.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A14.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G25.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A74.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E32.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G70.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E12.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A96.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E18.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C58.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G09.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C32.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E43.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G91.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A53.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C35.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C67.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E80.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A29.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A23.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C73.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C07.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C38.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C87.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A03.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G41.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E24.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C85.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A88.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C44.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A80.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A15.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E08.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C03.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E93.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C18.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G72.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G44.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A62.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C86.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C91.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A31.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C53.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E57.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E85.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C16.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E22.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A83.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G33.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G14.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C77.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E47.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A63.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C88.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C81.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A93.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C01.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G28.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G77.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E63.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G37.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C28.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A40.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C46.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G75.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C62.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C71.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G06.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E69.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G49.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A61.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G57.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G02.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A25.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G76.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E10.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G52.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G83.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G32.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E02.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E77.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G74.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C10.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A16.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G78.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E42.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C57.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E04.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A08.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G34.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A82.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A18.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E16.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A42.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G73.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G15.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G90.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G59.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E64.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G63.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E96.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E84.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A70.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C08.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G92.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C72.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G64.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C13.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G46.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A72.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A81.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G17.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E31.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C19.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C83.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C31.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C36.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C02.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E33.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A06.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E76.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G50.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A64.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E27.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G45.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E11.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E49.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C26.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E79.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A85.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E44.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E65.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C78.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G20.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G94.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G71.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A59.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C22.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A48.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A84.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A65.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G86.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A39.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C95.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C66.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C74.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E92.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E21.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E87.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A22.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C90.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E48.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G04.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C84.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G93.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C82.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G36.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A54.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G03.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E36.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G21.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E17.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G39.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A52.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G19.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E72.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A32.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A44.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C47.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E03.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E62.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E56.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E60.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G24.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C27.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A38.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G88.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E34.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A26.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E58.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C42.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E06.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G67.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C20.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A37.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A02.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G55.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E13.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C21.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A41.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A76.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A57.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E90.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E30.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A89.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G10.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E50.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E66.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C52.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A71.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C94.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E07.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E54.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A11.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G35.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E95.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G96.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C40.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C43.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A21.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E28.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G58.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A56.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G62.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G68.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E38.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C05.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E61.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E01.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A77.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C65.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G60.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A17.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A60.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E55.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E81.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A34.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A58.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A13.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E46.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G51.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C14.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E94.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G47.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A95.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G87.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A10.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A78.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G23.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C51.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A94.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C76.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E59.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E78.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A20.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E53.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C34.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E41.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E83.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C12.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E29.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G48.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G22.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E73.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G81.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E09.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C79.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A90.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E89.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C68.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A19.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G79.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C41.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G26.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C70.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E75.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C04.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C25.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A79.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A28.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G95.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A92.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G61.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A49.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C56.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A27.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G29.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C15.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E39.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C54.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E14.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G43.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A33.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G38.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E05.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G80.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G12.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C33.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C24.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A12.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE5E25.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C75.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C30.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C63.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C17.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU3C92.vcf.gz",
            "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRUE7G07.vcf.gz",
        ]
    )
)

# Core statement to fire parallelised above-mentionned function on f_ldir
# parmap.starmap(mp_vcf, list(zip(f_ldir)), l_df_gb, l_df, pm_processes=64)

l_df = list()
for e in list(f_ldir)[:10]:
    mp_vcf(e)

# parmap.starmap(mp_vcf, list(zip(f_ldir)), l_df_gb, l_df, pm_processes=snakemake.threads)

# Concat output
concat_df = pd.concat(list(l_df))
print(concat_df)


######

# Read coverage

coverage_df = pd.read_csv(
    "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/COVERAGE/LanexHGSVCpool2NEW/MERGE/merge_coverage.txt",
    sep="\t",
)
# coverage_df = pd.read_csv(snakemake.input.coverage, sep="\t")

# Read ashleys predictions, rename & sort
ashleys_predictions = (
    pd.read_csv(
        "/scratch/tweber/DATA/MC_DATA/STOCKS/2023-11-09-HW3YVAFX5/LanexHGSVCpool2NEW/cell_selection/labels.tsv",
        sep="\t",
    )
    # pd.read_csv(list(snakemake.input.ashleys_predictions)[0], sep="\t")
    .rename(
        {
            "cell": "QUERY_SAMPLE_CELL",
            "prediction": "ashleys_prediction",
            "probability": "ashleys_probability",
        },
        axis=1,
    ).sort_values(by="QUERY_SAMPLE_CELL")
)
ashleys_predictions["QUERY_SAMPLE_CELL"] = ashleys_predictions[
    "QUERY_SAMPLE_CELL"
].str.replace(".sort.mdup.bam", "")

# Read MC predictions & rename columns
mc_predictions = (
    pd.read_csv(
        "/scratch/tweber/DATA/MC_DATA/STOCKS/2023-11-09-HW3YVAFX5/LanexHGSVCpool2NEW/counts/LanexHGSVCpool2NEW.info_raw",
        skiprows=13,
        sep="\t",
    )
    # pd.read_csv(list(snakemake.input.mc_predictions)[0], skiprows=13, sep="\t")
    .rename({"cell": "QUERY_SAMPLE_CELL"}, axis=1)[
        ["QUERY_SAMPLE_CELL", "mapped", "good", "pass1"]
    ].rename(
        {
            "mapped": "reads_mapped",
            "good": "reads_used",
            "pass1": "mc_coverage_compliant",
        },
        axis=1,
    )
)

# Merge ashleys predictions & MC bool column
ashleys_mc_predictions = pd.merge(
    ashleys_predictions, mc_predictions, on="QUERY_SAMPLE_CELL"
)
ashleys_mc_predictions.loc[
    (ashleys_mc_predictions["ashleys_prediction"] == 1)
    & (ashleys_mc_predictions["mc_coverage_compliant"] == 1),
    "Used/Not used in MC",
] = 1

# Compute final column to know if used or not in the pipeline
ashleys_mc_predictions["Used/Not used in MC"] = ashleys_mc_predictions[
    "Used/Not used in MC"
].fillna(0)
ashleys_mc_predictions["Used/Not used in MC"] = ashleys_mc_predictions[
    "Used/Not used in MC"
].astype(int)


########

# Compute to check if working with SNP displaying QUAL above 10

qual_cutoff = 10
merge_df = pd.merge(
    concat_df.groupby("QUERY_SAMPLE_CELL")["ID"]
    .nunique()
    .rename("QUAL >= 0")
    .reset_index(),
    concat_df.loc[concat_df["QUAL"] >= qual_cutoff]
    .groupby("QUERY_SAMPLE_CELL")["ID"]
    .nunique()
    .rename("QUAL >= {}".format(str(qual_cutoff)))
    .reset_index(),
    on="QUERY_SAMPLE_CELL",
    how="outer",
).sort_values(by="QUERY_SAMPLE_CELL")


#######

# THIS PART IS TO FILL THE TABLE WITH MISSING CELLS WHEN THIS IS THE CASE
# SHOULD BE REMOVED IN THE FUTUR


# CHECKME
def findstem(arr):
    """Retrieve most common substring among a list of strings

    Args:
        arr (_type_): _description_

    Returns:
        _type_: _description_
    """

    # Determine size of the array
    n = len(arr)

    # Take first word from array
    # as reference
    s = arr[0]
    l = len(s)

    res = ""

    for i in range(l):
        for j in range(i + 1, l + 1):
            # generating all possible substrings
            # of our reference string arr[0] i.e s
            stem = s[i:j]
            k = 1
            for k in range(1, n):
                # Check if the generated stem is
                # common to all words
                if stem not in arr[k]:
                    break

            # If current substring is present in
            # all strings and its length is greater
            # than current result
            if k + 1 == n and len(res) < len(stem):
                res = stem

    return res


# Retrieve most common string
common = findstem(
    [os.path.basename(e).replace(".vcf.gz", "") for e in f_ldir]
    # [os.path.basename(e).replace(".vcf.gz", "") for e in list(snakemake.input.vcf)]
)
index = common[-1]
common = common[:-1]

# vcf_list = list(concat_df.QUERY_SAMPLE_CELL.unique())
# empty_samples_list = [
#     {
#         "QUERY_SAMPLE_CELL": "{}{}".format(common, str(e)),
#         "QUAL >= 0": np.nan,
#         "QUAL >= {}".format(str(qual_cutoff)): np.nan,
#     }
#     for e in list(range((int(index) * 100) + 1, (int(index) * 100) + 97))
#     if "{}{}".format(common, str(e)) not in vcf_list
# ]

# merge_df = pd.concat([merge_df, pd.DataFrame(empty_samples_list)]).sort_values(
#     by="QUERY_SAMPLE_CELL"
# )

########


# Counts SNPs / sample & pivot table

concat_counts_df = (
    concat_df.groupby("QUERY_SAMPLE_CELL")["SAMPLE"]
    .value_counts()
    .rename("Counts")
    .reset_index()
)
concat_counts_df["Rank"] = concat_counts_df.groupby("QUERY_SAMPLE_CELL")[
    "SAMPLE"
].transform(lambda r: range(1, len(r) + 1))

merge_df = pd.merge(
    merge_df,
    pd.pivot_table(
        concat_counts_df.loc[concat_counts_df["Rank"] <= 3],
        index="QUERY_SAMPLE_CELL",
        columns=["Rank"],
        values=["SAMPLE"],
        aggfunc=lambda x: " ".join(x),
    ),
    on="QUERY_SAMPLE_CELL",
    how="left",
)
merge_df = pd.merge(
    merge_df,
    pd.pivot_table(
        concat_counts_df.loc[concat_counts_df["Rank"] <= 3],
        index="QUERY_SAMPLE_CELL",
        columns=["Rank"],
    ),
    on="QUERY_SAMPLE_CELL",
    how="left",
)

concat_counts_df = (
    concat_df.loc[concat_df["QUAL"] >= qual_cutoff]
    .groupby("QUERY_SAMPLE_CELL")["SAMPLE"]
    .value_counts()
    .rename("Counts")
    .reset_index()
)
concat_counts_df["Rank"] = concat_counts_df.groupby("QUERY_SAMPLE_CELL")[
    "SAMPLE"
].transform(lambda r: range(1, len(r) + 1))

merge_df = pd.merge(
    merge_df,
    pd.pivot_table(
        concat_counts_df.loc[concat_counts_df["Rank"] <= 3],
        index="QUERY_SAMPLE_CELL",
        columns=["Rank"],
        values=["SAMPLE"],
        aggfunc=lambda x: " ".join(x),
    ),
    on="QUERY_SAMPLE_CELL",
    how="left",
)
merge_df = pd.merge(
    merge_df,
    pd.pivot_table(
        concat_counts_df.loc[concat_counts_df["Rank"] <= 3],
        index="QUERY_SAMPLE_CELL",
        columns=["Rank"],
    ),
    on="QUERY_SAMPLE_CELL",
    how="left",
)
merge_df = pd.merge(merge_df, coverage_df, on="QUERY_SAMPLE_CELL", how="left")
merge_df = pd.merge(
    merge_df, ashleys_mc_predictions, on="QUERY_SAMPLE_CELL", how="left"
)

pd.options.display.max_columns = 50
# merge_df.columns = [
#     "QUERY_SAMPLE_CELL",
#     "QUAL>=0",
#     "QUAL>={}".format(str(qual_cutoff)),
#     "QUAL>=0_SAMPLE_1",
#     "QUAL>=0_SAMPLE_2",
#     "QUAL>=0_SAMPLE_3",
#     "QUAL>=0_Counts_1",
#     "QUAL>=0_Counts_2",
#     "QUAL>=0_Counts_3",
#     "QUAL>={}_SAMPLE_1".format(str(qual_cutoff)),
#     "QUAL>={}_SAMPLE_2".format(str(qual_cutoff)),
#     "QUAL>={}_SAMPLE_3".format(str(qual_cutoff)),
#     "QUAL>={}_Counts_1".format(str(qual_cutoff)),
#     "QUAL>={}_Counts_2".format(str(qual_cutoff)),
#     "QUAL>={}_Counts_3".format(str(qual_cutoff)),
#     "meandepth",
#     "ashleys_prediction",
#     "ashleys_probability",
#     "reads_mapped",
#     "reads_used",
#     "mc_coverage_compliant",
#     "Used/Not used in MC",
# ]

# ##########

# # Compare with DAVID results

# merge_df["meandepth"] = merge_df["meandepth"] * 100

# hgsvc_batch = merge_df["QUERY_SAMPLE_CELL"].values.tolist()[0].split("xpool")[1][0]
# ground_truth_path = "/g/korbel2/weber/MosaiCatcher_files/POOLING/POOLING_POOL{index}/HGSVCxpool{index}/all"

# david_df = pd.DataFrame(
#     [
#         {
#             "QUERY_SAMPLE_CELL": os.path.basename(f).split("_")[0],
#             "DAVID_Sample": os.path.basename(f).split("_sorted_")[1].replace(".sort.mdup.bam", ""),
#         }
#         for f in os.listdir(ground_truth_path.format(index=hgsvc_batch))
#         if f.endswith(".bam")
#     ]
# )
# merge_df = pd.merge(merge_df, david_df, on="QUERY_SAMPLE_CELL", how="left")

# ########

# Output file

merge_df.to_excel(
    "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/ANALYSIS/LanexHGSVCpool2NEW.xlsx",
    index=False,
)
# merge_df.to_excel(snakemake.output.stats, index=False)
print(merge_df)
