import multiprocessing as mp
import parmap
import sys, os
import pandas as pd
import gzip

## List of inputs
# filelist_isec = ["/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_NA19443/0002.vcf", "/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_HG00109/0002.vcf"]
ref = pd.read_csv("/scratch/tweber/DATA/1000G_SNV_with_GT/BCFTOOLS_CLEAN/merge.vcf.gz", compression="gzip", sep="\t", nrows=10000)
ref = pd.read_parquet("/scratch/tweber/DATA/1000G_SNV_with_GT/BCFTOOLS_CLEAN/merge.parquet")
# ref = pd.read_csv(snakemake.input.ref, compression="gzip", sep="\t", nrows=10000)
print(ref)

# vcf = snakemake.input.vcf[0]
vcf = "/scratch/tweber/DATA/1000G_SNV_with_GT/GENOTYPING/HG00268_401.vcf.gz"

# Nb of rows to skip at the beginning of the VCF
skip = len([e.decode("ISO-8859–1").split('\n') for e in gzip.open("{}".format(vcf, "rb")) if e.decode("ISO-8859–1").split('\n')[0].startswith("##")])
# chrom = vcf.split("/")[-3]
sample = vcf.split("/")[-2].split("_")[-1]
df = pd.read_csv(vcf, compression="gzip", skiprows=skip, sep="\t")
df["#CHROM"] = df["#CHROM"].astype(str)
df["#CHROM"] = df["#CHROM"].str.replace("chr", "")
try:
    df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + ":" + df["REF"] + ":" + df["ALT"]
except:
    print(vcf)
print(df["ID"])

merge_df = pd.merge(ref, df, on="ID")

merge_df.groupby("SAMPLE")["ID"].nunique().sort_values(ascending=False)

# # Create dataframe from list of dict
# final_df = pd.DataFrame(list(m_l)).sort_values(by="Overlap_size", ascending=False)
# # Output the complete DF
# ## COMMENT: useful to have complete info at the chrom level
# final_df.to_csv(snakemake.output.detailed, index=False, sep="\t")

# # Simplify and sum at the sample level + output 
# final_df.groupby("Sample")["Overlap_size"].sum().sort_values(ascending=False).to_csv(snakemake.output.summary, index=True, sep="\t")

