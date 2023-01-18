import multiprocessing as mp
import parmap
import sys, os
import pandas as pd
import gzip
import numpy as np


## List of inputs
ref = pd.read_csv(snakemake.input.ref, compression="gzip", sep="\t")

# Prepare multiprocessing instances
m = mp.Manager()
l_df_gb = m.list()
l_df = m.list()


def mp_vcf(vcf, l_df_gb, l_df):

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
    df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + ":" + df["REF"] + ":" + df["ALT"]

    # Merge ref & df
    merge_df = pd.merge(ref, df, on="ID", how="right")

    merge_df["SAMPLE"] = merge_df["SAMPLE"].fillna("NAN")
    merge_df["QUERY_SAMPLE_CELL"] = sample

    # Groupby to count SNPs number / sample, then sort values
    merge_df_gb = merge_df.groupby("SAMPLE")["ID"].nunique().sort_values(ascending=False).reset_index()
    merge_df_gb["QUERY_SAMPLE_CELL"] = sample
    # Compute the rank
    merge_df_gb["Rank"] = list(range(1, 1 + merge_df_gb.shape[0]))
    # Compute total nb of SNPs
    merge_df_gb["Total_SNP"] = df.shape[0]

    l_df_gb.append(merge_df_gb)
    l_df.append(merge_df)


# Process list of cells aplhabetically
f_ldir = sorted(list(snakemake.input.vcf))

# Core statement to fire parallelised above-mentionned function on f_ldir
parmap.starmap(mp_vcf, list(zip(f_ldir)), l_df_gb, l_df, pm_processes=snakemake.threads)

# Concat output
concat_df = pd.concat(list(l_df))


######

# Read coverage

coverage_df = pd.read_csv(snakemake.input.coverage, sep="\t")

# Read ashleys predictions, rename & sort
ashleys_predictions = (
    pd.read_csv(list(snakemake.input.ashleys_predictions)[0], sep="\t")
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
ashleys_predictions["QUERY_SAMPLE_CELL"] = ashleys_predictions["QUERY_SAMPLE_CELL"].str.replace(".sort.mdup.bam", "")

# Read MC predictions & rename columns
mc_predictions = (
    pd.read_csv(list(snakemake.input.mc_predictions)[0], skiprows=13, sep="\t")
    .rename({"cell": "QUERY_SAMPLE_CELL"}, axis=1)[["QUERY_SAMPLE_CELL", "mapped", "good", "pass1"]]
    .rename(
        {
            "mapped": "reads_mapped",
            "good": "reads_used",
            "pass1": "mc_coverage_compliant",
        },
        axis=1,
    )
)

# Merge ashleys predictions & MC bool column
ashleys_mc_predictions = pd.merge(ashleys_predictions, mc_predictions, on="QUERY_SAMPLE_CELL")
ashleys_mc_predictions.loc[
    (ashleys_mc_predictions["ashleys_prediction"] == 1) & (ashleys_mc_predictions["mc_coverage_compliant"] == 1),
    "Used/Not used in MC",
] = 1

# Compute final column to know if used or not in the pipeline
ashleys_mc_predictions["Used/Not used in MC"] = ashleys_mc_predictions["Used/Not used in MC"].fillna(0)
ashleys_mc_predictions["Used/Not used in MC"] = ashleys_mc_predictions["Used/Not used in MC"].astype(int)


########

# Compute to check if working with SNP displaying QUAL above 10

qual_cutoff = 10
merge_df = pd.merge(
    concat_df.groupby("QUERY_SAMPLE_CELL")["ID"].nunique().rename("QUAL >= 0").reset_index(),
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
common = findstem([os.path.basename(e).replace(".vcf.gz", "") for e in list(snakemake.input.vcf)])
index = common[-1]
common = common[:-1]

vcf_list = list(concat_df.QUERY_SAMPLE_CELL.unique())
empty_samples_list = [
    {
        "QUERY_SAMPLE_CELL": "{}{}".format(common, str(e)),
        "QUAL >= 0": np.nan,
        "QUAL >= {}".format(str(qual_cutoff)): np.nan,
    }
    for e in list(range((int(index) * 100) + 1, (int(index) * 100) + 97))
    if "{}{}".format(common, str(e)) not in vcf_list
]

merge_df = pd.concat([merge_df, pd.DataFrame(empty_samples_list)]).sort_values(by="QUERY_SAMPLE_CELL")

########


# Counts SNPs / sample & pivot table

concat_counts_df = concat_df.groupby("QUERY_SAMPLE_CELL")["SAMPLE"].value_counts().rename("Counts").reset_index()
concat_counts_df["Rank"] = concat_counts_df.groupby("QUERY_SAMPLE_CELL")["SAMPLE"].transform(lambda r: range(1, len(r) + 1))

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
    concat_df.loc[concat_df["QUAL"] >= qual_cutoff].groupby("QUERY_SAMPLE_CELL")["SAMPLE"].value_counts().rename("Counts").reset_index()
)
concat_counts_df["Rank"] = concat_counts_df.groupby("QUERY_SAMPLE_CELL")["SAMPLE"].transform(lambda r: range(1, len(r) + 1))

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
merge_df = pd.merge(merge_df, ashleys_mc_predictions, on="QUERY_SAMPLE_CELL", how="left")

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

##########

# Compare with DAVID results

merge_df["meandepth"] = merge_df["meandepth"] * 100

hgsvc_batch = merge_df["QUERY_SAMPLE_CELL"].values.tolist()[0].split("xpool")[1][0]
ground_truth_path = "/g/korbel2/weber/MosaiCatcher_files/POOLING/POOLING_POOL{index}/HGSVCxpool{index}/all"

david_df = pd.DataFrame(
    [
        {
            "QUERY_SAMPLE_CELL": os.path.basename(f).split("_")[0],
            "DAVID_Sample": os.path.basename(f).split("_sorted_")[1].replace(".sort.mdup.bam", ""),
        }
        for f in os.listdir(ground_truth_path.format(index=hgsvc_batch))
        if f.endswith(".bam")
    ]
)
merge_df = pd.merge(merge_df, david_df, on="QUERY_SAMPLE_CELL", how="left")

########

# Output file

merge_df.to_excel(snakemake.output.stats, index=False)
print(merge_df)
