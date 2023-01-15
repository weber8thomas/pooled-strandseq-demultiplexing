import pandas as pd
import sys, os

df = pd.read_excel(
    "/scratch/tweber/DATA/1000G_SNV_with_GT/OTF/ANALYSIS/HGSVCxPool3.xlsx"
)
df["Cell"] = df["QUERY_SAMPLE_CELL"].str.replace(
    "iPSCsxcontrolx01PE20", "", regex=False
)
print(df.loc[df["Used/Not used in MC"] == 1])
print(df.columns)
print(
    df.loc[df["Used/Not used in MC"] == 1][
        ["QUERY_SAMPLE_CELL", "QUAL >= 10", "('SAMPLE', 1)"]
    ]
)

david_df = pd.DataFrame(
    [
        {
            "QUERY_SAMPLE_CELL": os.path.basename(f),
            "Cell": os.path.basename(f).split("_")[0].replace("HGSVCxpool3x01PE20", ""),
        }
        for f in os.listdir(
            "/g/korbel2/weber/MosaiCatcher_files/POOLING/POOLING_POOL3/HGSVCxpool3/all"
        )
        if f.endswith(".bam")
    ]
)
print(david_df)
