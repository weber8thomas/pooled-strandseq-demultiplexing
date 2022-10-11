import multiprocessing as mp
import parmap
import sys, os
import pandas as pd

# Multiprocessing manager + shared list
m = mp.Manager() 
m_l = m.list()

## List of inputs
# filelist_isec = ["/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_NA19443/0002.vcf", "/scratch/tweber/DATA/NYGC_1KGP/ISEC/chr2/HG00512_HG00109/0002.vcf"]
filelist_isec = list(snakemake.input)

def gather_output(vcf, m_l):
    """
    MP function to read each VCF one-to-one join comparison (only shared variants are keeped)
    Add a dictionnary to the shared list with overlap size and SNP list
    """

    # Nb of rows to skip at the beginning of the VCF
    skip = len([e.split('\n') for e in open("{file}".format(file=vcf), "r") if e.split('\n')[0].startswith("##")])
    # chrom = vcf.split("/")[-3]
    sample = vcf.split("/")[-2].split("_")[-1]
    df = pd.read_csv(vcf, skiprows=skip, sep="\t")
    try:
        df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + "_" + df["REF"] + "_" + df["ALT"]
    except:
        print(vcf)
    m_l.extend(
        [{
            # "Chrom" : chrom,
            "Sample" : sample,
            "SNP_list" : df["ID"].sort_values().tolist(),
            "Overlap_size" : len(df["ID"].sort_values().tolist()),
        }]
    )

# Run the function on the input list
parmap.starmap(gather_output, list(zip(filelist_isec)), m_l, pm_pbar=True, pm_processes=snakemake.threads)

# Create dataframe from list of dict
final_df = pd.DataFrame(list(m_l)).sort_values(by="Overlap_size", ascending=False)
# Output the complete DF
## COMMENT: useful to have complete info at the chrom level
final_df.to_csv(snakemake.output.detailed, index=False, sep="\t")

# Simplify and sum at the sample level + output 
final_df.groupby("Sample")["Overlap_size"].sum().sort_values(ascending=False).to_csv(snakemake.output.summary, index=True, sep="\t")

