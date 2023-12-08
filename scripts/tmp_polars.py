import polars as pl
import os
import glob

folder_path = (
    "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/COMPARING_STRANDSEQ_TO_1KG"
)
file_pattern = os.path.join(folder_path, "*.tsv")
column_names = [
    "Concatenated",
    "CHROM_1",
    "POS_1",
    "Dot",
    "REF_1",
    "ALT_1",
    "Quality",
    "Dot2",
    "Info",
    "Format",
    "Genotype",
    "CellID",
    "CHROM_2",
    "POS_2",
    "REF_2",
    "ALT_2",
    "SampleGenotype",
]

# Lazy DataFrame list
lazy_dfs = []

for file_path in sorted(glob.glob(file_pattern)):
    print(file_path)
    sample_name = file_path.split("_")[-1].replace(".tsv", "")
    lazy_df = pl.scan_csv(
        file_path, has_header=False, separator="\t", new_columns=column_names
    ).with_columns([pl.lit(sample_name).alias("SampleName")])

    lazy_dfs.append(lazy_df)

# Concatenate and collect
final_df = pl.concat(lazy_dfs).collect()

final_df_sorted = final_df.sort("CellID")
parquet_path = "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/COMPARING_STRANDSEQ_TO_1KG/LanexHGSVCpool2NEW.parquet"
os.makedirs(parquet_path, exist_ok=True)
final_df_sorted.write_parquet(
    parquet_path, use_pyarrow=True, pyarrow_options={"partition_cols": ["CellID"]}
)
