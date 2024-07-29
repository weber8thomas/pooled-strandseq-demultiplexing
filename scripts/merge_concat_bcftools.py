import pandas as pd

# fldir = [
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/NA18989.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/NA19320.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/NA19331.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/NA19836.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/NA20355.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG02282.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG02554.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG02666.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG02769.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG02953.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG03452.txt.gz",
# ]
# fldir = [
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/NA19625.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG02442.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/NA18560.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG00288.txt.gz",
#     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/HG01811.txt.gz",
# ]

# print(list(snakemake.input))

l = [
    pd.read_csv(e, compression="gzip", sep="\t")
    # for j, e in enumerate(list(fldir))
    for j, e in enumerate(list(snakemake.input))
]
# pd.concat(l).to_csv(
#     # "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CONCAT_TAB/Random_for_sanity_check/merge.txt.gz",
#     sep="\t",
#     compression="gzip",
#     mode="w",
#     index=False,
# )
pd.concat(l).to_csv(
    snakemake.output[0], sep="\t", compression="gzip", mode="w", index=False
)
