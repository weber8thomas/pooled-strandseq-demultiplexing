import multiprocessing as mp
import parmap
import sys, os
import pandas as pd

m = mp.Manager() 
m_l = m.list()

# filelist_isec = ["/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_NA19443/0002.vcf", "/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_HG00109/0002.vcf"]
filelist_isec = list(snakemake.input)

def gather_output(vcf, m_l):
    skip = len([e.split('\n') for e in open("{file}".format(file=vcf), "r") if e.split('\n')[0].startswith("##")])
    chrom = vcf.split("/")[-3]
    sample = vcf.split("/")[-2].split("_")[-1]
    df = pd.read_csv(vcf, skiprows=skip, sep="\t")
    # print(df)
    try:
        df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + "_" + df["REF"] + "_" + df["ALT"]
    except:
        print(vcf)
    m_l.extend(
        [{
            "Chrom" : chrom,
            # "Chrom" : "chr2",
            "Sample" : sample,
            # "Sample" : "TEST",
            "SNP_list" : df["ID"].sort_values().tolist(),
            "Overlap_size" : len(df["ID"].sort_values().tolist()),
        }]
    )

parmap.starmap(gather_output, list(zip(filelist_isec)), m_l, pm_pbar=True, pm_processes=snakemake.threads)
# print(list(m_l))
final_df = pd.DataFrame(list(m_l)).sort_values(by="Overlap_size", ascending=False)
# final_df = pd.DataFrame(list(m_l))
final_df.to_csv(snakemake.output.detailed, index=False, sep="\t")
final_df.groupby("Sample")["Overlap_size"].sum().sort_values(ascending=False).to_csv(snakemake.output.detailed, index=False, sep="\t")

