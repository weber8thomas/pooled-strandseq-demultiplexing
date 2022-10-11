import pandas as pd
from tqdm import tqdm
import parmap 
import multiprocessing as mp
import os, subprocess, sys

# m = mp.Manager()
# # l = m.list()
# folder = "/scratch/tweber/DATA/1000G_SNV_with_GT"
# filelist = [e for e in os.listdir(folder) if e.endswith(".vcf.gz")]
# header = subprocess.Popen("zcat {folder}/{file} | head -n 5000 | grep -P  '^#CHROM'".format(folder=folder, file=filelist[0]), shell=True, stdout=subprocess.PIPE)
# samples = header.communicate()[0].decode("utf-8").strip().split("\t")[9:]
# dict_samples = {e:j for j,e in enumerate(samples)}
# hgsvc_samples = ["HG00155", "HG00186", "HG00233", "HG00267", "HG00423", "HG00480", "HG00673", "HG00735", "HG00766", "HG01071", "HG01106", "HG01255", "HG01273", "HG01493", "HG01609", "HG01887", "HG01928", "HG01952", "HG01978", "HG02068", "HG02132", "HG02135", "HG02148", "HG02602", "HG02698", "HG02723", "HG02975", "HG03035", "HG03056", "HG03200", "HG03456", "HG03519", "HG03522", "HG03576", "HG03654", "HG03719", "HG03816", "HG03866", "HG03942", "HG03968", "HG04115", "HG04157", "HG04204", "NA10837", "NA10857", "NA18537", "NA18940", "NA18965", "NA19100", "NA19139", "NA19440", "NA19443", "NA19721", "NA19724", "NA19727", "NA19730", "NA19902", "NA20128", "NA20536", "NA20846",]
# pool_samples = ["HG00268", "HG00512", "HG00514", "HG00731", "HG01352", "HG02059", "HG02818", "NA12878", "NA19239", "NA19240",]


# sample_list = list(sorted(list(set([e for e in samples] + hgsvc_samples + pool_samples))))[:950]

# put = pool_samples
# for j, e in tqdm(enumerate(list(snakemake.input))):
l = [pd.read_csv(e, compression='gzip', sep="\t", names=["ID", "AC", "AF", "SAMPLE"]) for j, e in tqdm(enumerate(list(snakemake.input)))]
pd.concat(l).to_csv(snakemake.output[0], sep="\t", compression="gzip", mode="w", index=False)

# l = list()
# i = 300
# for j, e in tqdm(enumerate(list(snakemake.input))):
#     l.append(pd.read_csv(e, compression='gzip', sep="\t", names=["ID", "AC", "AF", "SAMPLE"]))
#     print(j)
#     if j % i == 0 and j != (len(list(snakemake.input)) -1):
#         if j == i:
#             print("J == {}".format(i))
#             pd.concat(l).to_csv(snakemake.output[0], sep="\t", compression="gzip", mode="w", index=False)
#         elif j > i:
#             print("J > {}".format(i))
#             pd.concat(l).to_csv(snakemake.output[0], sep="\t", compression="gzip", mode="a", header=None, index=False)
#         l = list()
#     elif j == (len(list(snakemake.input)) -1):
#         print("LAST")
#         pd.concat(l).to_csv(snakemake.output[0], sep="\t", compression="gzip", mode="a", header=None, index=False)

