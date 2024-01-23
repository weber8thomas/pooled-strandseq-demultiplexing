# In[1]:


import pandas as pd
import os, sys


# In[40]:


l_samples = [
    "2021-07-29-HWYJ2AFX2/HGSVCxpool1x01",
    "2021-08-03-H22VWAFX3/HGSVCxpool2x02",
    "2021-08-03-H22VWAFX3/HGSVCxpool3x01",
    "2023-11-09-HW5NFAFX5/HGSVCpool1NEW",
    "2023-11-09-HW3YVAFX5/LanexHGSVCpool2NEW",
]

for output_dir in [
    "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/FINAL_RESULTS/without_sanity_check",
    "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/FINAL_RESULTS/with_sanity_check",
]:
    print(output_dir)

    # output_dir = (
    #     "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/FINAL_RESULTS/with_sanity_check"
    # )

    # In[3]:

    # sample = "HGSVCpool1NEW"

    for iteration in l_samples:
        print(iteration)
        year = iteration.split("-")[0]
        sample = iteration.split("/")[1]
        print(sample)
        index = "PE20" if year == "2021" else "iTRU"
        os.makedirs(f"{output_dir}/{sample}", exist_ok=True)

        # In[41]:

        ref_path = f"/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz"
        # ref_path = "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CONCAT_TAB/HGSVCxpool3x01/merge.txt.gz"
        # ref_path = "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CONCAT_TAB/PSEUDOPOOL/merge.txt.gz"
        ref = pd.read_csv(
            ref_path,
            compression="gzip",
            sep="\t",
        )
        ref["GlobalSample"] = ref_path.split("/")[-2]

        if "with_sanity_check" in output_dir:
            sanity_check_path = "/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/BCFTOOLS_CONCAT_TAB/Random_for_sanity_check/merge.txt.gz"
            sanity_check_ref = pd.read_csv(
                sanity_check_path, compression="gzip", sep="\t"
            )
            sanity_check_ref["GlobalSample"] = sanity_check_path.split("/")[-2]

            ref = pd.concat([ref, sanity_check_ref])

        ref["ID"] = "chr" + ref["ID"]
        ref

        # In[42]:

        ref.groupby(["GlobalSample", "SAMPLE"])["ID"].count()

        # In[43]:

        ref_count = ref.groupby(["GlobalSample", "SAMPLE"])["ID"].count().reset_index()
        ref_count.to_csv(
            f"{output_dir}/{sample}/reference_file_SNP_counts.tsv", sep="\t"
        )
        ref_count.head()

        # In[44]:

        # sanity_check_ref.SAMPLE.unique()

        # In[45]:

        import glob

        final_vcf = list()
        for vcf_input in glob.glob(
            f"/scratch/tweber/DATA/MC_DATA/DEMULTIPLEXING_POOLS/SNP_DISCOVERY/{sample}/*.vcf"
        ):
            cell_line = vcf_input.split("/")[-1].split("_")[0].replace(".vcf", "")
            vcf = pd.read_csv(vcf_input, skiprows=255, sep="\t")
            vcf["cell_line"] = cell_line
            vcf["ID"] = (
                vcf["#CHROM"]
                + ":"
                + vcf["POS"].astype(str)
                + ":"
                + vcf["REF"]
                + ":"
                + vcf["ALT"]
            )
            vcf = vcf.drop([cell_line.split(f"{index}")[0]], axis=1)
            final_vcf.append(vcf)
        final_vcf = pd.concat(final_vcf)
        final_vcf.to_csv(
            f"{output_dir}/{sample}/full_set_of_SNPs_called.tsv.gz",
            sep="\t",
            compression="gzip",
        )

        final_vcf.head()

        # In[46]:

        merge_df = pd.merge(final_vcf, ref, on="ID", how="inner")
        merge_df.to_csv(
            f"{output_dir}/{sample}/join_ref_to_SNPs_called.tsv.gz",
            sep="\t",
            compression="gzip",
        )
        merge_df.head()

        # In[47]:

        pd.options.display.max_rows = 200
        pd.options.display.max_columns = 200
        mosaicatcher_stats = pd.concat(
            [
                pd.read_csv(
                    "/scratch/tweber/DATA/MC_DATA/STOCKS/{}/counts/{}.info_raw".format(
                        sample, sample.split("/")[-1]
                    ),
                    sep="\t",
                    skiprows=13,
                )
                for sample in l_samples
            ]
        )
        mosaicatcher_stats["cell"] = mosaicatcher_stats["cell"].str.replace(
            ".sort.mdup.bam", ""
        )
        mosaicatcher_stats.head()

        # In[48]:

        pd.options.display.max_rows = 200
        pd.options.display.max_columns = 200
        ashleys_labels = pd.concat(
            [
                pd.read_csv(
                    f"/scratch/tweber/DATA/MC_DATA/STOCKS/{sample}/cell_selection/labels.tsv",
                    sep="\t",
                )
                for sample in l_samples
            ]
        )
        ashleys_labels["cell"] = ashleys_labels["cell"].str.replace(
            ".sort.mdup.bam", ""
        )
        ashleys_labels.head()

        # In[49]:

        combine_ashleys_mc_stats = pd.merge(
            ashleys_labels, mosaicatcher_stats, on=["sample", "cell"], how="inner"
        )
        combine_ashleys_mc_stats.head()

        # In[50]:

        gb_sample_count = (
            merge_df.rename({"cell_line": "cell"}, axis=1)
            .groupby(["cell", "GlobalSample", "SAMPLE"])["ID"]
            .count()
            .reset_index()
        )
        gb_sample_count_stats = gb_sample_count.rename(columns={"cell_line": "cell"})
        gb_sample_count_stats = pd.merge(
            gb_sample_count_stats, combine_ashleys_mc_stats, on=["cell"], how="inner"
        )
        gb_sample_count_stats = gb_sample_count_stats.loc[
            gb_sample_count_stats["prediction"] == 1
        ]
        gb_sample_count_stats["cell"] = gb_sample_count_stats["cell"].apply(
            lambda r: r.split(index)[1]
        )
        gb_sample_count_stats.to_csv(
            f"{output_dir}/{sample}/groupby_sample_counts_with_stats.tsv.gz",
            sep="\t",
            compression="gzip",
        )
        gb_sample_count_stats.head()

        # In[14]:

        gb_sample_count_stats.cell.nunique()

        # In[54]:

        # pivot_table_stats

        # In[51]:

        # pd.options.display.max_columns = None
        # pd.options.display.max_rows = None
        pivot_table_stats = pd.pivot_table(
            gb_sample_count_stats,
            columns=["GlobalSample", "SAMPLE"],
            index="cell",
            values="ID",
        ).fillna(0)
        import scipy.stats as stats

        pivot_table_stats_zscore = pivot_table_stats.apply(
            lambda x: stats.zscore(x), axis=1
        )
        pivot_table_stats_zscore.to_csv(
            f"{output_dir}/{sample}/pivot_table_zscore_norm.tsv.gz",
            sep="\t",
            compression="gzip",
        )
        # pivot_table_stats_zscore.index = [e.split(index)[1] for e in pivot_table_stats_zscore.index]
        pivot_table_stats_zscore.head()

        # In[69]:

        import seaborn as sns
        import matplotlib.pyplot as plt

        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(40, 20))
        ax = sns.heatmap(pivot_table_stats_zscore, cmap="Reds", vmin=-1, vmax=6)
        # ax.set_title("Matched SNP nb in PseudoPool (z-score adjusted)")
        # ax.set_xlabel("Sample")
        # ax.set_ylabel("Cell Line")
        ax.figure.savefig(f"{output_dir}/{sample}/heatmap_zscore_cell_per_sample.png")
        ax

        # In[17]:

        import seaborn as sns
        import matplotlib.pyplot as plt

        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(40, 20))
        ax = sns.clustermap(pivot_table_stats_zscore, cmap="Reds", vmin=-1, vmax=6)
        # ax.set_title("Matched SNP nb in PseudoPool (z-score adjusted)")
        # ax.set_xlabel("Sample")
        # ax.set_ylabel("Cell Line")
        ax.figure.savefig(
            f"{output_dir}/{sample}/clustermap_zscore_cell_per_sample.png"
        )
        ax

        # In[25]:

        ["-".join(e) for e in pivot_table_stats_zscore.columns]

        # In[52]:

        import seaborn as sns
        import matplotlib.pyplot as plt

        sns.set_theme(style="whitegrid")
        # Assuming pivot_table_stats_zscore is your data

        # Create the clustermap
        ax = sns.clustermap(pivot_table_stats_zscore, cmap="Reds", vmin=-1, vmax=6)

        # Set the size of the figure
        width, height = 40, 30  # You can adjust these values as needed
        ax.fig.set_size_inches(width, height)
        plt.setp(
            ax.ax_heatmap.get_xticklabels(),
            rotation=90,
            ha="right",
            rotation_mode="anchor",
            fontsize=10,
        )  # Adjust fontsize as needed

        # Get the number of labels (assuming they are the same as the number of columns in your data)
        num_labels = pivot_table_stats_zscore.shape[1]

        # Set the x-tick labels manually
        ax.ax_heatmap.set_xticks([x + 0.5 for x in range(num_labels)])
        ax.ax_heatmap.set_xticklabels(
            ["-".join(e) for e in pivot_table_stats_zscore.columns],
            rotation=90,
            ha="right",
            rotation_mode="anchor",
            fontsize=16,
        )
        print(ax.ax_heatmap.get_xticklabels())

        # Save the figure
        output_file = f"{output_dir}/{sample}/clustermap_zscore_cell_per_sample.png"
        ax.savefig(output_file)

        # ax.set_title("Matched SNP nb in PseudoPool (z-score adjusted)")
        # ax.set_xlabel("Sample")
        # ax.set_ylabel("Cell Line")
        # ax.figure.savefig(f"{output_dir}/{sample}/clustermap_zscore_cell_per_sample.png")
        ax

        # In[ ]:

        pivot_table_stats_zscore

        # In[53]:

        import dash_bio

        f = dash_bio.Clustergram(
            data=pivot_table_stats_zscore.values,
            column_labels=pivot_table_stats_zscore.columns.tolist(),
            row_labels=pivot_table_stats_zscore.index.tolist(),
            color_threshold={"row": 150, "col": 700},
            height=1800,
            width=1700,
            color_map=["white", "red"],
            display_range=[-1, 6],
        )

        f.write_html(
            f"{output_dir}/{sample}/clustermap_zscore_cell_per_sample_interactive.html"
        )
        f

        # In[63]:

        pivot_table_stats_zscore_wt_multiindex = pivot_table_stats_zscore.copy()
        pivot_table_stats_zscore_wt_multiindex.columns = (
            pivot_table_stats_zscore.columns.droplevel()
        )
        pivot_table_stats_zscore_wt_multiindex = (
            pivot_table_stats_zscore_wt_multiindex.reset_index().rename(
                {"index": "cell"}, axis=1
            )
        )
        pivot_table_stats_zscore_wt_multiindex.head()

        # In[64]:

        # pivot_table_stats_zscore.melt()
        melt_df = pd.melt(
            pivot_table_stats_zscore_wt_multiindex,
            id_vars=["cell"],
            value_vars=[
                e
                for e in pivot_table_stats_zscore_wt_multiindex.columns
                if e not in ["cell"]
            ],
        )
        melt_df.head()

        # In[21]:

        def top_n_samples(group, n=3):
            return group.sort_values(by="ID", ascending=False).head(n)

        # tmp_debug = gb_sample_count_stats.loc[(gb_sample_count_stats["sample"] == "GM19836x01") & (gb_sample_count_stats["prediction"] == 1)].sort_values(by=["cell", "ID"], ascending=[True, False])
        pd.options.display.max_rows = None
        gb_sample_count_stats
        tmp_debug = (
            gb_sample_count_stats.groupby("cell")
            .apply(top_n_samples, n=1)
            .reset_index(drop=True)
        )
        tmp_debug = pd.merge(tmp_debug, melt_df, on=["cell", "SAMPLE"], how="inner")
        tmp_debug[["cell", "SAMPLE", "ID", "value", "probability", "good"]]
        tmp_debug.head(100)
        # # Generate a color for each unique cell
        # unique_cells = tmp_debug['cell'].unique()
        # colors = [f"background-color: rgb({200 + i*20 % 55}, {220 - (i*20) % 55}, {200 + (i*30) % 55})" for i in range(len(unique_cells))]
        # color_map = dict(zip(unique_cells, colors))

        # # Apply the colors
        # def apply_row_colors(row):
        #     return [color_map[row['cell']]] * len(row)

        # styled_df = tmp_debug.drop(["bam"], axis=1).style.apply(apply_row_colors, axis=1)

        # # Display the styled DataFrame in Jupyter Notebook
        # styled_df

        # In[ ]:

        # In[22]:

        tmp_debug.SAMPLE.nunique()

        # In[3]:

        # Rpy2

        get_ipython().run_line_magic("load_ext", "rpy2.ipython")

        # In[66]:

        metadata = pd.read_csv(
            "../20130606_g1k_3202_samples_ped_population.txt", sep=" "
        )
        metadata.head()

        # In[25]:

        pivot_table_stats_zscore.head()

        # In[67]:

        pivot_table_stats_zscore_wt_multiindex = (
            pivot_table_stats_zscore_wt_multiindex.set_index("cell")
        )
        pivot_table_stats_zscore_wt_multiindex.head()

        # In[ ]:

        get_ipython().run_cell_magic(
            "R",
            "-i pivot_table_stats_zscore_wt_multiindex -i metadata -i ref_count -w 1800 -h 1500",
            '\nlibrary(ComplexHeatmap)\nlibrary(circlize)\n\nset.seed(123) # for reproducibility\nordered_metadata <- metadata[match(colnames(pivot_table_stats_zscore_wt_multiindex), metadata$SampleID), ]\n\n# Map GlobalSample to the SAMPLE in pivot_table_stats_zscore\nglobal_sample_annotation <- ref_count[match(colnames(pivot_table_stats_zscore_wt_multiindex), ref_count$SAMPLE), "GlobalSample"]\n\n\n# Create HeatmapAnnotation objects for Population, Superpopulation, and GlobalSample\ncol_annotation <- HeatmapAnnotation(df = ordered_metadata[c("Population", "Superpopulation")],\n                                    GlobalSample = global_sample_annotation)\n\n\n# Convert the pandas DataFrame to an R matrix\nmat <- as.matrix(pivot_table_stats_zscore_wt_multiindex)\n#print(head(mat))\n\n#png("clustermap_zscore_cell_per_sample_with_annotations.png", width = 1800, height = 1500)\n\n# Creating the heatmap\nHeatmap(mat, \n        name = "z-score", \n        col = colorRamp2(c(-1, 6), c("white", "red")),\n        top_annotation = col_annotation,\n        cluster_rows = TRUE, \n        cluster_columns = TRUE,\n        show_row_names = TRUE,\n        show_column_names = TRUE,\n       row_names_gp = gpar(fontsize = 8)) # Adjust fontsize as needed\n\n#dev.off()\n',
        )

        # In[4]:

        get_ipython().run_cell_magic(
            "R",
            "-i pivot_table_stats_zscore_wt_multiindex -i metadata -i ref_count -w 1800 -h 1500",
            '\n# Ensure that the InteractiveComplexHeatmap package is installed\nif (!requireNamespace("InteractiveComplexHeatmap", quietly = TRUE)) {\n    if (!requireNamespace("BiocManager", quietly = TRUE)) {\n        install.packages("BiocManager")\n    }\n    BiocManager::install("InteractiveComplexHeatmap")\n}\n\nlibrary(ComplexHeatmap)\nlibrary(InteractiveComplexHeatmap)\n\n# Assuming pivot_table_stats_zscore_wt_multiindex is loaded in your R environment\n# as well as the metadata and ref_count\n\n# Your existing code for the heatmap setup\nordered_metadata <- metadata[match(colnames(pivot_table_stats_zscore_wt_multiindex), metadata$SampleID), ]\nglobal_sample_annotation <- ref_count[match(colnames(pivot_table_stats_zscore_wt_multiindex), ref_count$SAMPLE), "GlobalSample"]\ncol_annotation <- HeatmapAnnotation(df = ordered_metadata[c("Population", "Superpopulation")],\n                                    GlobalSample = global_sample_annotation)\nmat <- as.matrix(pivot_table_stats_zscore_wt_multiindex)\n\n# Your ComplexHeatmap creation\nht = Heatmap(mat, \n             name = "z-score", \n             col = colorRamp2(c(-1, 6), c("white", "red")),\n             top_annotation = col_annotation,\n             cluster_rows = TRUE, \n             cluster_columns = TRUE,\n             show_row_names = TRUE,\n             show_column_names = TRUE,\n             row_names_gp = gpar(fontsize = 8)) # Adjust fontsize as needed\n\n# Now to make it interactive\nht = draw(ht, heatmap_legend_side = "bot", annotation_legend_side = "bot")\ninteractive_heatmap(ht)\n',
        )

        # In[28]:

        # ashleys_labels.loc[(ashleys_labels["sample"] == "GM19836x01") & (ashleys_labels["prediction"] == 1)]

        # In[29]:

        # pivot_table_stats_zscore.columns = pivot_table_stats_zscore.columns.droplevel()
        # pivot_table_stats_zscore_melt = pivot_table_stats_zscore.melt(
        #    ignore_index=False,
        #    var_name="SAMPLE",
        #    value_name="SNP nb (z-score adjusted)",
        # ).reset_index()
        # pivot_table_stats_zscore_melt["Sample_to_find"] = pivot_table_stats_zscore_melt["cell"].apply(lambda r: r.split("x")[0].replace("GM", "NA"))
        ## retrieve the Sample with highest SNP nb
        # pivot_table_stats_zscore_melt = pivot_table_stats_zscore_melt.loc[pivot_table_stats_zscore_melt.groupby("cell")["SNP nb (z-score adjusted)"].idxmax()]
        ##pivot_table_stats_zscore_melt.loc[pivot_table_stats_zscore_melt["Sample_to_find"] == pivot_table_stats_zscore_melt["SAMPLE"], "Match"] = True
        ##pivot_table_stats_zscore_melt.loc[pivot_table_stats_zscore_melt["Sample_to_find"] != pivot_table_stats_zscore_melt["SAMPLE"], "Match"] = False
        # pivot_table_stats_zscore_melt

        # In[30]:

        # false_match = pivot_table_stats_zscore_melt.loc[pivot_table_stats_zscore_melt["Match"] == False]
        # false_match = pd.merge(false_match, gb_sample_count_stats, on=["cell", "SAMPLE"], how="inner")
        # false_match

        # In[31]:

        # false_match.cell.values

        # In[32]:

        # pd.options.display.max_rows = 210
        # tmp_false = gb_sample_count_stats.loc[gb_sample_count_stats["cell"].isin(false_match.cell.values.tolist())].sort_values(by=["cell", "ID"], ascending=[True, False])
        #
        # def top_n_samples(group, n=3):
        #    return group.sort_values(by='ID', ascending=False).head(n)
        #
        # tmp_false = tmp_false.groupby('cell').apply(top_n_samples, n=3).reset_index(drop=True)
        #
        #
        # def highlight_cells(x):
        #    colors = ['background-color: yellow', 'background-color: lightgreen', 'background-color: lightblue']
        #    return [colors[i % len(colors)] if x.name in top_samples['cell'].unique() else '' for i in range(len(x))]
        #
        ## Generate a color for each unique cell
        # unique_cells = tmp_false['cell'].unique()
        # colors = [f"background-color: rgb({200 + i*20 % 55}, {220 - (i*20) % 55}, {200 + (i*30) % 55})" for i in range(len(unique_cells))]
        # color_map = dict(zip(unique_cells, colors))
        #
        ## Apply the colors
        # def apply_row_colors(row):
        #    return [color_map[row['cell']]] * len(row)
        #
        # styled_df = tmp_false.drop(["bam"], axis=1).style.apply(apply_row_colors, axis=1)
        #
        ## Display the styled DataFrame in Jupyter Notebook
        # styled_df
