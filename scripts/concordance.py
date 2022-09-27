import pandas as pd
import sys, os
import gzip
import collections
import parmap
import multiprocessing as mp

m = mp.Manager()
m_l = m.list()

# file = "/scratch/tweber/DATA/NYGC_1KGP/GENOTYPING/HG01493/chr10.vcf.gz"
input_folder_gt = "/scratch/tweber/DATA/NYGC_1KGP/GENOTYPING/HG01493"
filelist_gt = [f for f in os.listdir(input_folder_gt) if f.endswith(".vcf.gz")]

input_folder_aggr = "/scratch/tweber/DATA/NYGC_1KGP/AGGR"
filelist_aggr = [f for f in os.listdir(input_folder_aggr) if f.endswith(".tsv.gz")]


skip = len([e.decode("ISO-8859–1").split('\n') for e in gzip.open("{folder}/{file}".format(folder=input_folder_gt, file=filelist_gt[0]), "rb") if e.decode("ISO-8859–1").split('\n')[0].startswith("##")])

df = pd.concat([pd.read_csv("{folder}/{file}".format(folder=input_folder_gt, file=file), skiprows=skip, compression='gzip',sep="\t") for file in filelist_gt])

chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX"]
df["#CHROM"] = pd.Categorical(df["#CHROM"], categories=chroms, ordered=True)
df["POS"] = df["POS"].astype(int)
df = df.sort_values(by=["#CHROM", "POS"])
df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + "_" + df["REF"] + "_" + df["ALT"]
print(df)

genotyped = set(df["ID"].unique().tolist())


def process(file, m_l):
    df_aggr = pd.read_csv("{folder}/{file}".format(folder=input_folder_aggr, file=file), compression='gzip',sep="\t")
    aggr_dict = df_aggr.groupby("SAMPLE")["ID"].apply(list).to_dict()
    # aggr_dict = aggr_dict["HG00512"]
    aggr_dict_stats = collections.defaultdict(dict)
    m_l.extend([
        {
            "Sample" : k,
            "Overlap_count" : len(genotyped.intersection(set(v))),
            "Overlap_list" : genotyped.intersection(set(v))
        }
        for k,v in aggr_dict.items()
    ])


parmap.starmap(process, list(zip(filelist_aggr)), m_l, pm_pbar=True, pm_processes=23)

stats = pd.DataFrame(list(m_l)).sort_values(by="Overlap_count", ascending=False)
stats["%_overlap"] = stats["Overlap_count"] / len(genotyped)
stats["%_overlap"] = stats["%_overlap"].round(4)
# stats = stats.loc[stats["Overlap_count"] > 0]
# stats = stats.loc[stats["Sample"] == "HG01493"]


print(genotyped)

print(stats)
