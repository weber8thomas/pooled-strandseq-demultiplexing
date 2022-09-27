import pandas as pd
from ast import literal_eval
import multiprocessing as mp
import parmap
import sys, os

m = mp.Manager()
l_out_vcf = m.list()
l_out_df = m.list()

l = list(snakemake.input)

def process_df(f, l_out_vcf, l_out_df):
    tmp_df = pd.read_csv(f, compression='gzip', sep="\t")
    tmp_df["alleles"] = tmp_df["alleles"].apply(literal_eval)
    tmp_df["alleles"] = tmp_df["alleles"].apply(lambda r: [e[0] for e in r if e][:2])
    tmp_df["#CHROM"] = tmp_df["locus"].apply(lambda r: r.split(':')[0])
    tmp_df["POS"] = tmp_df["locus"].apply(lambda r: r.split(':')[1])
    tmp_df["REF"] = tmp_df["alleles"].apply(lambda r: r[0])
    tmp_df["ALT"] = tmp_df["alleles"].apply(lambda r: r[1])
    tmp_df["ID"] = tmp_df["ID"].apply(lambda r: "_".join([e[0] if j > 0 else e for j, e in enumerate(r.split('_'))]))
    # tmp_df["INFO"] = "AC={AC};AF={AF}".format(AC=tmp_df["AC"].astype(str), AF=tmp_df["AF"].astype(str))
    tmp_df["SAMPLE"] = os.path.basename(f).strip(".tsv.gz")
    tmp_df["INFO"] = "AC=" + tmp_df["AC"].astype(str) + ";AF=" + tmp_df["AF"].astype(str) + ";SAMPLE=" + tmp_df["SAMPLE"]
    tmp_df["QUAL"] = "."
    tmp_df["FILTER"] = "PASS"

    tmp_df_vcf = tmp_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]
    l_out_vcf.append(tmp_df_vcf)

    tmp_df_df = tmp_df[["ID", "SAMPLE", "AC", "AF"]]
    l_out_df.append(tmp_df_df)

parmap.starmap(process_df, list(zip(l)), l_out_vcf, l_out_df, pm_pbar=True, pm_processes=snakemake.threads)

df_vcf_final = pd.concat(list(l_out_vcf))

# GETTING CHR LIST & ORDER CATEGORICALLY
chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX"]
df_vcf_final["#CHROM"] = pd.Categorical(df_vcf_final["#CHROM"], categories=chroms, ordered=True)
df_vcf_final["POS"] = df_vcf_final["POS"].astype(int)
df_vcf_final = df_vcf_final.sort_values(by=["#CHROM", "POS"])

df_final = pd.concat(list(l_out_df))


vcfheader = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
##bcftools_viewCommand=view /gpfs/commons/home/usevani/compbio/CCDG/Project_CCDG_14151_B01_GRM_WGS/final_annotated_vcfs/tmp_dir_annt/CCDG_14151_B01_GRM_WGS_2020-08-05_chr10.recalibrated_variants.annotated.repeats.vcf.gz; Date=Mon Sep 28 17:31:05 2020
"""

with open(snakemake.output.vcf, "w") as w:
    w.write(vcfheader)
df_vcf_final.to_csv(snakemake.output.vcf, mode="a", index=False, sep="\t")

df_final.to_csv(snakemake.output.tsv, compression='gzip', sep="\t", index=False)


# print(df_final.groupby("SAMPLE")["ID"].apply(list).to_dict())