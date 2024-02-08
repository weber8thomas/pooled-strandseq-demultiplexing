import subprocess
import os
import collections


# CONFIG PART
configfile: "config/config.yaml"


results_folder = config["results_folder"]
onekgp_vcf_folder = config["onekgp_vcf_folder"]
bam_folder = config["bam_folder"]

# list samples in bam folder
samples = [
    e
    for e in os.listdir(bam_folder)
    if os.path.isdir(bam_folder + "/" + e) and e not in ["config", "log"] and "pool" in e

]



print(samples)

sample_list = config["samples_pooled"]
# Rename 1KG samples from GM to NA
sample_list = {k: [e.replace("GM", "NA") for e in v] for k, v in sample_list.items()}
print(sample_list)
new_check = True if "NEW" in samples[0] else False
if new_check:
    sample_list = {k:v for k,v in sample_list.items() if "NEW" in k}
else:
    sample_list = {k:v for k,v in sample_list.items() if "NEW" not in k}
sample_list = {s:v for k,v in sample_list.items() for s in samples if k in s}
print(sample_list)

# list 1kgp files
onekgp_filelist = [
    e
    for e in os.listdir(
        onekgp_vcf_folder + "/HGSVC2_1KG_3202_VCF_FILES_RENAMED_SYMLINK"
    )
    if e.endswith(".vcf.gz")
]
print(onekgp_vcf_folder + "/HGSVC2_1KG_3202_VCF_FILES_RENAMED_SYMLINK")
# print(onekgp_filelist)

# Retrieve header
header = subprocess.Popen(
    "zcat {folder}/{file} | head -n 100 | grep -P  '^#CHROM'".format(
        folder=onekgp_vcf_folder + "/HGSVC2_1KG_3202_VCF_FILES_RENAMED_SYMLINK",
        file=onekgp_filelist[0],
    ),
    shell=True,
    stdout=subprocess.PIPE,
)

# List of 1kgp samples & dict with index nb
samples = header.communicate()[0].decode("utf-8").strip().split("\t")[9:]
dict_samples = {e: j for j, e in enumerate(samples)}
all_samples_to_genotype = [sub_e for k, v in sample_list.items() for sub_e in v]
# print(samples)
# print(all_samples_to_genotype)


import random

random.seed(43)


# Retrieve list of cells for each Pool sample to be genotyped
cell_dict = collections.defaultdict(list)
for sample in sorted(os.listdir(bam_folder)):
    if sample not in ["config", "log"] and sample in sample_list.keys():
        l_dir = [
            k
            for k in os.listdir(bam_folder + "/" + sample + "/bam")
            if k.endswith(".sort.mdup.bam")
        ]
        random.shuffle(l_dir)
        for j, e in enumerate(l_dir):
            # if j < 100:
            if e.endswith(".bam"):
                cell_dict[sample].append(e.replace(".sort.mdup.bam", ""))


# random.seed(42)
random_selection = random.sample(
    [x for x in samples if x not in all_samples_to_genotype], 10
)
# random_selection = random.sample([x for x in samples if x not in all_samples_to_genotype], 10)


# List of pool samples to process
samples = list(cell_dict.keys())


# sample_list["GM19317x02"] = ["NA19317"] + random_selection
# sample_list["Random_for_sanity_check"] = random_selection


print(sample_list)
# exit()

def get_mem_mb_heavy(wildcards, attempt):
    """
    To adjust resources in the rules
    """
    mem_avail = [32, 64, 96, 128, 192, 256]
    return mem_avail[attempt - 1] * 1000


def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    """
    mem_avail = [2, 4, 8, 16, 32, 64]
    return mem_avail[attempt - 1] * 1000


rule all:
    input:
        # expand(
        #     "{results_folder}/ANALYSIS/{sample}.xlsx",
        #     results_folder=results_folder,
        #     sample=samples,
        # ),
        [
            expand(
                "{results_folder}/SNP_DISCOVERY/{sample}/{cell_id}.vcf",
                # "{results_folder}/GENOTYPING_OTF/{sample}/{cell_id}.vcf.gz",
                results_folder=results_folder,
                sample=sample,
                cell_id=cell_dict[sample],
            )
            for sample in samples
        ],
        # TODO: CHANGE THIS TO A POOL OUTPUT THAT CAN BE REUSED FOR ALL CORRESPONDING PLATES
        [expand("{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz", results_folder=results_folder, sample=list(sample_list.keys())) for sample in samples],
        [
            expand(
                "{results_folder}/COVERAGE/{sample}/MERGE/merge_coverage.txt",
                results_folder=results_folder,
                sample=sample,
            )
            for sample in samples
        ],
        # [expand(
        #     "{bam_folder}/{sample}/cell_selection/labels.tsv",
        #     bam_folder=bam_folder,
        #     sample=sample,
        # ) for sample in samples],
        # [expand(
        #     "{bam_folder}/{sample}/counts/{sample}.info_raw",
        #     bam_folder=bam_folder,
        #     sample=sample,
        # ) for sample in samples],


rule create_single_1kgp_file_rare_het_intermediate:
    input:
        expand(
            "{onekgp_vcf_folder}/{file}",
            onekgp_vcf_folder=config["onekgp_vcf_folder"],
            file=onekgp_filelist,
        ),
    output:
        "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005.vcf.gz",
    conda:
        "envs/snp_genotyping.yaml"
    threads: 10
    resources:
        mem_mb=128000,
        time="30:00:00",
    params:
        maf=0.05,
    shell:
        """
        bcftools concat {input}  | \
        bcftools view --type snps -i "INFO/AC>0 & INFO/AF<{params.maf}" | \
        bcftools sort -T {config[tmp_data_folder]} | \
        bcftools norm -d none --threads {threads} -Oz -o {output}
        """


rule create_single_1kgp_file_rare_het_decompose:
    input:
        "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005.vcf.gz",
    output:
        "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
    conda:
        "envs/snp_genotyping.yaml"
    threads: 1
    resources:
        mem_mb=128000,
        time="30:00:00",
    params:
        maf=0.05,
    shell:
        """
        vt decompose {input} | \
        bcftools annotate -x "^INFO/AC,INFO/AF" -Oz -o {output}
        """


rule retrieve_rare_variants_bcftools:
    input:
        onekgp=ancient(
            expand(
                "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
                onekgp_vcf_folder=onekgp_vcf_folder,
            )
        ),
        onekgp_index=ancient(
            expand(
                "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
                onekgp_vcf_folder=onekgp_vcf_folder,
            )
        ),
    output:
        "{results_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz",
    conda:
        "envs/snp_genotyping.yaml"
    threads: 1
    resources:
        mem_mb=32000,
        time="10:00:00",
    params:
        index=lambda wc: dict_samples[wc.onekgp_sample],
        col_index=lambda wc: dict_samples[wc.onekgp_sample] + 1 + 9,  # +1 for python index +9 for the nine first cols
        af=0.05,
    shell:
        "bcftools view --threads {threads} -i 'INFO/AC>0 & INFO/AF<{params.af} & GT[{params.index}]=\"het\"' --type snps {input.onekgp} | cut -f 1-9,{params.col_index} | bgzip > {output}"


rule concat_filter_bcftools:
    input:
        "{results_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz",
    output:
        "{results_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt.gz",
    log:
        "{results_folder}/log/BCFTOOLS_CLEAN/{onekgp_sample}.log",
    threads: 1
    resources:
        mem_mb=2000,
    run:
        import pandas as pd
        import gzip
        import os, sys

        skip = len(
            [
                e.decode("ISO-8859–1").split("\n")
                for e in gzip.open("{}".format(input[0], "rb"))
                if e.decode("ISO-8859–1").split("\n")[0].startswith("##")
            ]
        )
        sample = os.path.basename(input[0]).replace(".vcf.gz", "")
        df = pd.read_csv(input[0], compression="gzip", skiprows=skip, sep="\t")
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
        df["AC"] = df["INFO"].apply(lambda r: r.split(";")[0].replace("AC=", ""))
        df["AF"] = df["INFO"].apply(lambda r: r.split(";")[1].replace("AF=", ""))
        df["SAMPLE"] = wildcards.onekgp_sample
        df = df[["ID", "AC", "AF", "SAMPLE"]]
        df.to_csv(output[0], compression="gzip", sep="\t", index=False)
        # 'bcftools view --type snps {input} | cut -f 3,8 | sed "s/$/\\t{wildcards.sample}/g;s/AC=//g;s/AF=//g" |  grep -v "^#" | tail -n+2 | tr ";" "\\t" | bgzip > {output}'


localrules: merge_concat_bcftools_tab, merge_concat_bcftools_vcf


rule merge_concat_bcftools_tab:
    input:
        lambda wc: expand(
            "{results_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt.gz",
            results_folder=results_folder,
            onekgp_sample=sample_list[wc.sample],
        ),
    output:
        "{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
    log:
        "{results_folder}/log/BCFTOOLS_CLEAN/{sample}/merge.log",
    conda:
        # "/g/korbel2/weber/miniconda3"
        "envs/python_env.yaml"
    threads: 1
    resources:
        mem_mb=16000,
    script:
        "scripts/merge_concat_bcftools.py"


rule merge_concat_bcftools_vcf:
    input:
        vcf=lambda wc: expand(
            "{results_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz",
            results_folder=results_folder,
            onekgp_sample=sample_list[wc.sample],
        ),
        vcf_index=lambda wc: expand(
            "{results_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz.tbi",
            results_folder=results_folder,
            onekgp_sample=sample_list[wc.sample],
        ),
    output:
        "{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf",
    log:
        "{results_folder}/log/BCFTOOLS_CONCAT/{sample}/merge.log",
    conda:
        "envs/snp_genotyping.yaml"
    threads: 1
    resources:
        mem_mb=get_mem_mb_heavy,
    shell:
        "bcftools merge {input.vcf} | bcftools sort --temp-dir /tmpdata/tweber/ | bcftools norm -d none | bcftools view -h -H > {output}"


rule samtools_coverage:
    input:
        bam=lambda wc: expand(
            "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam",
            bam_folder=bam_folder,
            sample=wc.sample,
        ),
        bai=lambda wc: expand(
            "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam.bai",
            bam_folder=bam_folder,
            sample=wc.sample,
        ),
    output:
        coverage="{results_folder}/COVERAGE/{sample}/{cell_id}.txt",
    # conda:
    #     "envs/snp_genotyping.yaml"
    envmodules:
        "SAMtools/1.14-GCC-11.2.0",
    resources:
        mem_mb=2000,
        time="01:00:00",
    shell:
        "samtools coverage {input.bam} -o {output.coverage}"


rule merge_coverage:
    input:
        lambda wc: expand(
            "{results_folder}/COVERAGE/{sample}/{cell_id}.txt",
            results_folder=results_folder,
            sample=wc.sample,
            cell_id=cell_dict[wc.sample],
        ),
    output:
        "{results_folder}/COVERAGE/{sample}/MERGE/merge_coverage.txt",
    log:
        "{results_folder}/log/COVERAGE/{sample}/merge_coverage.log",
    run:
        import pandas as pd

        l = list()
        chroms = ["chr" + str(e) for e in list(range(1, 23))] + ["chrX"]
        for file in list(input):
            df = pd.read_csv(file, sep="\t")
            df = df.loc[df["#rname"].isin(chroms)]
            df["QUERY_SAMPLE_CELL"] = file.split("/")[-1].replace(".txt", "")
            df["QUERY_SAMPLE"] = file.split("/")[-2]
            l.append(df)
        coverage_df = pd.concat(l)
        coverage_df_summary = (
            coverage_df.groupby(["QUERY_SAMPLE_CELL"])["meandepth"].mean().reset_index()
        )
        coverage_df_summary.to_csv(output[0], index=False, sep="\t")


# NOTE: use custom fasta that covers contigs present in all libraries from pseudopool
# fasta location: /g/korbel2/weber/MosaiCatcher_files/EXTERNAL_DATA/refgenomes_human_local/hg38.complete.fa
rule regenotype_SNVs:
    input:
        bam=lambda wc: expand(
            "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam",
            bam_folder=bam_folder,
            sample=wc.sample,
        ),
        bai=lambda wc: expand(
            "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam.bai",
            bam_folder=bam_folder,
            sample=wc.sample,
        ),
        # sites="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz",
        # sites_index="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz.tbi",
        fasta=config["fasta_ref"],
        fasta_index="{fasta_ref}.fai".format(fasta_ref=config["fasta_ref"]),
        # merge_table="{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
    output:
        vcf="{results_folder}/SNP_DISCOVERY/{sample}/{cell_id}.vcf",
    log:
        vcf="{results_folder}/log/SNP_DISCOVERY/{sample}/{cell_id}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="20:00:00",
    threads: 1
    conda:
        "envs/snp_genotyping.yaml"
    # envmodules:
    #     "freebayes/1.3.6-foss-2021b-R-4.1.2",
    #     "BCFtools/1.16-GCC-11.3.0",
    shell:
        """
        freebayes -f {input.fasta} --exclude-unobserved-genotype {input.bam} > {output.vcf}
        """
        # """
        # freebayes \
        #     -f {input.fasta} \
        #     -@ {input.sites} \
        #     --only-use-input-alleles \
        #     --genotype-qualities \
        #     --exclude-unobserved-genotypes \
        #     {input.bam} \
        # | bcftools view \
        #     --exclude-uncalled \
        #     --types snps \
        #     --genotype het \
        #     -Oz -o {output.vcf} 2> {log}
        # """


# rule analyse_isec:
#     input:
#         vcf=lambda wc: expand(
#             "{results_folder}/GENOTYPING_OTF/{sample}/{cell_id}.vcf.gz",
#             results_folder=results_folder,
#             sample=wc.sample,
#             cell_id=cell_dict[wc.sample],
#         ),
#         ref="{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
#         coverage="{results_folder}/COVERAGE/{sample}/MERGE/merge_coverage.txt",
#         ashleys_predictions=lambda wc: expand(
#             "{bam_folder}/{sample}/cell_selection/labels.tsv",
#             bam_folder=bam_folder,
#             sample=wc.sample,
#         ),
#         mc_predictions=lambda wc: expand(
#             "{bam_folder}/{sample}/counts/{sample}.info_raw",
#             bam_folder=bam_folder,
#             sample=wc.sample,
#         ),
#     output:
#         stats="{results_folder}/ANALYSIS/{sample}.xlsx",
#     conda:
#         "envs/python_env.yaml"
#     threads: 128
#     resources:
#         mem_mb=256000,
#         time="24:00:00",
#     script:
#         "scripts/gt_analysis.py"


rule index_vcf:
    input:
        "{file}.vcf.gz",
    output:
        "{file}.vcf.gz.tbi",
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        "tabix -p vcf {input}"


rule compress_vcf:
    input:
        "{file}.vcf",
    output:
        "{file}.vcf.gz",
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        "bgzip {input}"
