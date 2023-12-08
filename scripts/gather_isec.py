import multiprocessing as mp
import parmap
import sys, os
import pandas as pd
import gzip

## List of inputs
# filelist_isec = ["/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_NA19443/0002.vcf", "/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_HG00109/0002.vcf"]
# ref = pd.read_csv("/scratch/tweber/DATA/1000G_SNV_with_GT/BCFTOOLS_CLEAN/merge.vcf.gz", compression="gzip", sep="\t", nrows=10000)
# ref = pd.read_parquet("/scratch/tweber/DATA/1000G_SNV_with_GT/BCFTOOLS_CLEAN/merge.parquet")


ref = pd.read_csv(snakemake.input.ref, compression="gzip", sep="\t", nrows=10000)

# # vcf = snakemake.input.vcf[0]
# vcf = "/scratch/tweber/DATA/1000G_SNV_with_GT/GENOTYPING/HG00268_401.vcf.gz"

# # Nb of rows to skip at the beginning of the VCF
# skip = len([e.decode("ISO-8859–1").split('\n') for e in gzip.open("{}".format(vcf, "rb")) if e.decode("ISO-8859–1").split('\n')[0].startswith("##")])
# # chrom = vcf.split("/")[-3]
# sample = vcf.split("/")[-2].split("_")[-1]
# df = pd.read_csv(vcf, compression="gzip", skiprows=skip, sep="\t")
# df["#CHROM"] = df["#CHROM"].astype(str)
# df["#CHROM"] = df["#CHROM"].str.replace("chr", "")
# try:
#     df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + ":" + df["REF"] + ":" + df["ALT"]
# except:
#     print(vcf)
# print(df["ID"])

# merge_df = pd.merge(ref, df, on="ID")

# merge_df.groupby("SAMPLE")["ID"].nunique().sort_values(ascending=False)

# # Create dataframe from list of dict
# final_df = pd.DataFrame(list(m_l)).sort_values(by="Overlap_size", ascending=False)
# # Output the complete DF
# ## COMMENT: useful to have complete info at the chrom level
# final_df.to_csv(snakemake.output.detailed, index=False, sep="\t")

# # Simplify and sum at the sample level + output
# final_df.groupby("Sample")["Overlap_size"].sum().sort_values(ascending=False).to_csv(snakemake.output.summary, index=True, sep="\t")


l_df_gb = []
l_df = []


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
    # chrom = vcf.split("/")[-3]
    sample = vcf.split("/")[-1].replace(".vcf.gz", "")
    df = pd.read_csv(vcf, compression="gzip", skiprows=skip, sep="\t")
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

    # df = df.loc[df["QUAL"] >= 10]

    merge_df = pd.merge(ref, df, on="ID", how="inner")

    # merge_df["SAMPLE"] = merge_df["SAMPLE"].fillna("NAN")
    merge_df["QUERY_SAMPLE_CELL"] = sample

    merge_df_gb = (
        merge_df.groupby("SAMPLE")["ID"]
        .nunique()
        .sort_values(ascending=False)
        .reset_index()
    )
    merge_df_gb["QUERY_SAMPLE_CELL"] = sample
    merge_df_gb["Rank"] = list(range(1, 1 + merge_df_gb.shape[0]))
    merge_df_gb["Total_SNP"] = df.shape[0]

    l_df_gb.append(merge_df_gb)
    l_df.append(merge_df)


# f = snakemake.input.gt_sample_folder
# f_ldir = sorted([f"{f}/{e}" for e in os.listdir(f) if e.endswith(".vcf.gz")])
f_ldir = sorted(list(snakemake.input.vcf))
# f_ldir = ["/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/GENOTYPING_OTF/LanexHGSVCpool2NEW/LanexHGSVCpool2NEWiTRU1A06.vcf.gz"]
print(f_ldir)
# f_ldir = sorted(list(snakemake.input.vcf))
# parmap.starmap(mp_vcf, list(zip(f_ldir)), l_df_gb, l_df)

skip = len(
    [
        e.decode("ISO-8859–1").split("\n")
        for e in gzip.open("{}".format(f_ldir[0], "rb"))
        if e.decode("ISO-8859–1").split("\n")[0].startswith("##")
    ]
)

# df = pd.read_csv(f_ldir[0], compression="gzip", skiprows=skip, sep="\t")
# df["#CHROM"] = df["#CHROM"].astype(str)
# df["#CHROM"] = df["#CHROM"].str.replace("chr", "")
# df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + ":" + df["REF"] + ":" + df["ALT"]


for e in f_ldir:
    print(e)
    mp_vcf(e)


pd.options.display.max_rows = 250
concat_df = pd.concat(list(l_df))
concat_df["GT"] = concat_df["LanexHGSVCpool2NEW"].str.split(":", expand=True)[0]
concat_df = concat_df.sort_values(by=["ID"])
concat_df["QUERY_SAMPLE_CELL"] = concat_df["QUERY_SAMPLE_CELL"].apply(
    lambda x: x.split("/")[-1].replace(".vcf.gz", "")
)
# concat_df.loc[(concat_df["GT"] != "1/1") & (concat_df["GT"] != "0/0")]


# concat_counts_df = concat_df.loc[concat_df["QUAL"] >= 10].groupby("QUERY_SAMPLE_CELL")["SAMPLE"].value_counts().rename("Counts").reset_index()

concat_counts_df = (
    concat_df.groupby("QUERY_SAMPLE_CELL")["SAMPLE"]
    .value_counts()
    .rename("Counts")
    .reset_index()
)
concat_counts_df["Rank"] = concat_counts_df.groupby("QUERY_SAMPLE_CELL")[
    "SAMPLE"
].transform(lambda r: range(1, len(r) + 1))


coverage_df = pd.read_csv(
    "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/COVERAGE/LanexHGSVCpool2NEW/MERGE/merge_coverage.txt",
    sep="\t",
)


ashleys_predictions = (
    pd.read_csv(
        "/scratch/tweber/DATA/MC_DATA/STOCKS/2023-11-09-HW3YVAFX5/LanexHGSVCpool2NEW/cell_selection/labels.tsv",
        sep="\t",
    )
    .rename(
        {
            "cell": "QUERY_SAMPLE_CELL",
            "prediction": "ashleys_prediction",
            "probability": "ashleys_probability",
        },
        axis=1,
    )
    .sort_values(by="QUERY_SAMPLE_CELL")
)
# ashleys_predictions = pd.read_csv(snakemake.input.ashleys_predictions, sep="\t").rename({"cell": "QUERY_SAMPLE_CELL", "prediction" : "ashleys_prediction", "probability": "ashleys_probability"}, axis=1).sort_values(by="QUERY_SAMPLE_CELL")
# ashleys_predictions = pd.read_csv("/g/korbel2/weber/MosaiCatcher_files/POOLING/POOL1_RESEQ/HGSVCxPool1/cell_selection/labels.tsv", sep="\t").rename({"cell": "QUERY_SAMPLE_CELL", "prediction" : "ashleys_prediction", "probability": "ashleys_probability"}, axis=1).sort_values(by="QUERY_SAMPLE_CELL")
ashleys_predictions["QUERY_SAMPLE_CELL"] = ashleys_predictions[
    "QUERY_SAMPLE_CELL"
].str.replace(".sort.mdup.bam", "")
# mc_predictions = pd.read_csv("/g/korbel2/weber/MosaiCatcher_files/POOLING/POOL1_RESEQ/HGSVCxPool1/counts/HGSVCxPool1.info_raw", skiprows=13, sep="\t").rename({"cell": "QUERY_SAMPLE_CELL"}, axis=1)[["QUERY_SAMPLE_CELL", "mapped", "good", "pass1"]].rename({"mapped": "reads_mapped", "good" : "reads_used", "pass1": "mc_coverage_compliant"}, axis=1)
mc_predictions = (
    pd.read_csv(
        "/scratch/tweber/DATA/MC_DATA/STOCKS/2023-11-09-HW3YVAFX5/LanexHGSVCpool2NEW/counts/LanexHGSVCpool2NEW.info_raw",
        skiprows=13,
        sep="\t",
    )
    .rename({"cell": "QUERY_SAMPLE_CELL"}, axis=1)[
        ["QUERY_SAMPLE_CELL", "mapped", "good", "pass1"]
    ]
    .rename(
        {
            "mapped": "reads_mapped",
            "good": "reads_used",
            "pass1": "mc_coverage_compliant",
        },
        axis=1,
    )
)
# mc_predictions = pd.read_csv(snakemake.input.mc_predictions, skiprows=13, sep="\t").rename({"cell": "QUERY_SAMPLE_CELL"}, axis=1)[["QUERY_SAMPLE_CELL", "mapped", "good", "pass1"]].rename({"mapped": "reads_mapped", "good" : "reads_used", "pass1": "mc_coverage_compliant"}, axis=1)

ashleys_mc_predictions = pd.merge(
    ashleys_predictions, mc_predictions, on="QUERY_SAMPLE_CELL"
)
ashleys_mc_predictions.loc[
    (ashleys_mc_predictions["ashleys_prediction"] == 1)
    & (ashleys_mc_predictions["mc_coverage_compliant"] == 1),
    "Used/Not used in MC",
] = 1
ashleys_mc_predictions["Used/Not used in MC"] = ashleys_mc_predictions[
    "Used/Not used in MC"
].fillna(0)
ashleys_mc_predictions["Used/Not used in MC"] = ashleys_mc_predictions[
    "Used/Not used in MC"
].astype(int)


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


concat_counts_df = (
    concat_df.groupby("QUERY_SAMPLE_CELL")["SAMPLE"]
    .value_counts()
    .rename("Counts")
    .reset_index()
)
concat_counts_df["Rank"] = concat_counts_df.groupby("QUERY_SAMPLE_CELL")[
    "SAMPLE"
].transform(lambda r: range(1, len(r) + 1))


# First pivot table
pivot_table1 = pd.pivot_table(
    concat_counts_df.loc[concat_counts_df["Rank"] <= 3],
    index="QUERY_SAMPLE_CELL",
    columns=["Rank"],
    values=["SAMPLE"],
    aggfunc=lambda x: " ".join(x),
).reset_index()
pivot_table1.columns = [
    "_".join(map(str, col)).strip() for col in pivot_table1.columns.values
]
pivot_table1.columns = ["QUERY_SAMPLE_CELL", "SAMPLE_1", "SAMPLE_2", "SAMPLE_3"]


merge_df = pd.merge(merge_df, pivot_table1, on="QUERY_SAMPLE_CELL", how="left")


pivot_table2 = pd.pivot_table(
    concat_counts_df.loc[
        concat_counts_df["Rank"] <= 3, ["QUERY_SAMPLE_CELL", "SAMPLE", "Counts", "Rank"]
    ],
    index="QUERY_SAMPLE_CELL",
    columns=["Rank"],
    values=["Counts"],
).reset_index()
pivot_table2.columns = [
    "_".join(map(str, col)).strip() for col in pivot_table2.columns.values
]
pivot_table2.columns = [
    "QUERY_SAMPLE_CELL",
    "SAMPLE_1_SNP_nb",
    "SAMPLE_2_SNP_nb",
    "SAMPLE_3_SNP_nb",
]


merge_df = pd.merge(merge_df, pivot_table2, on="QUERY_SAMPLE_CELL", how="left")


merge_df = pd.merge(merge_df, coverage_df, on="QUERY_SAMPLE_CELL", how="left")
merge_df = pd.merge(
    merge_df, ashleys_mc_predictions, on="QUERY_SAMPLE_CELL", how="left"
)

# pd.options.display.max_columns = 50
# merge_df.columns = ['QUERY_SAMPLE_CELL', 'QUAL>=0', 'QUAL>={}'.format(str(qual_cutoff)), 'QUAL>=0_SAMPLE_1', 'QUAL>=0_SAMPLE_2', 'QUAL>=0_SAMPLE_3', 'QUAL>=0_Counts_1', 'QUAL>=0_Counts_2', 'QUAL>=0_Counts_3',
# 'QUAL>={}_SAMPLE_1'.format(str(qual_cutoff)), 'QUAL>={}_SAMPLE_2'.format(str(qual_cutoff)), 'QUAL>={}_SAMPLE_3'.format(str(qual_cutoff)), 'QUAL>={}_Counts_1'.format(str(qual_cutoff)), 'QUAL>={}_Counts_2'.format(str(qual_cutoff)), 'QUAL>={}_Counts_3'.format(str(qual_cutoff)),  'Average_Coverage', 'ashleys_prediction', 'ashleys_probability', 'reads_mapped', 'reads_used', 'mc_coverage_compliant', 'Used/Not used in MC']


# # merge_df = merge_df[
# #         [
# #              'QUERY_SAMPLE_CELL',   'QUAL>={}'.format(str(qual_cutoff)), 'QUAL>={}_SAMPLE_1'.format(str(qual_cutoff)), 'QUAL>={}_SAMPLE_2'.format(str(qual_cutoff)), 'QUAL>={}_SAMPLE_3'.format(str(qual_cutoff)), 'QUAL>={}_Counts_1'.format(str(qual_cutoff)), 'QUAL>={}_Counts_2'.format(str(qual_cutoff)), 'QUAL>={}_Counts_3'.format(str(qual_cutoff)),  'Average_Coverage', 'ashleys_prediction', 'ashleys_probability', 'reads_mapped', 'reads_used', 'mc_coverage_compliant', 'Used/Not used in MC'
# #         ]
# # ]
# merge_df["Average_Coverage"] = merge_df["Average_Coverage"]*100

# # merge_df = merge_df.loc[merge_df["mc_coverage_compliant"] == 1].sort_values(by="Average_Coverage")

merge_df.to_excel(snakemake.output.stats, index=False)
