import subprocess
import os
import collections


configfile: "config/config.yaml"

results_folder = config["results_folder"]
onekgp_vcf_folder = config["onekgp_vcf_folder"]
bam_folder = config["bam_folder"]
sample_list = config["samples_pooled"]

# results_folder = "/scratch/tweber/DATA/1000G_SNV_with_GT/OTF"
# onekgp_vcf_folder = "/scratch/tweber/DATA/1000G_SNV_with_GT/VCF"
# bam_folder = "/g/korbel2/weber/MosaiCatcher_files/POOLING/POOL1_RESEQ/HGSVCxPool1/bam"

# list 1kgp files 
onekgp_filelist = [e for e in os.listdir(onekgp_vcf_folder) if e.endswith(".vcf.gz")]

# Retrieve header
header = subprocess.Popen(
    "zcat {folder}/{file} | head -n 5000 | grep -P  '^#CHROM'".format(
        folder=onekgp_vcf_folder, file=onekgp_filelist[0]
    ),
    shell=True,
    stdout=subprocess.PIPE,
)
# List of 1kgp samples & dict with index nb
samples = header.communicate()[0].decode("utf-8").strip().split("\t")[9:]
dict_samples = {e: j for j, e in enumerate(samples)}

# pool1 = ["NA10864", "NA18530", "NA18537", "NA18941", "NA18965", "NA19145", "NA19202", "NA19436", "NA19474", "NA19705", "NA19721", "NA19730", "NA20128", "NA20510", "NA20862", "NA20869", "HG00097", "HG00179", "HG00438", "HG00735", "HG00766", "HG01106", "HG01273", "HG01611", "HG01881", "HG01916", "HG01978", "HG02523", "HG02615", "HG02662", "HG02812", "HG02945", "HG02972", "HG02978", "HG03098", "HG03480", "HG03771", "HG03804", "HG04115", "HG04157",]
# sample_list = pool_samples
# sample_list = pool1


cell_dict = collections.defaultdict(list)
for sample in os.listdir(bam_folder):
    if sample not in ["config", "log"]:
        for e in os.listdir(bam_folder + "/" + sample + "/bam"): 
            if e.endswith(".bam"):
                cell_dict[sample].append(e.replace(".sort.mdup.bam", ""))   
# sample_name = bam_folder.split("/")[-2]
print(cell_dict)
samples = list(cell_dict.keys())
# print(sample_name)

## DEV - PSEUDOPOOL TO GET SAMPLE & CELL NAMES INTO A DICT

# sample_target_list = [e for e in os.listdir(folder_target) if e.endswith(".bam")]
# sample_target_dict = collections.defaultdict(list)
# [
#     sample_target_dict[e.replace(".sort.mdup.bam", "").split("_")[0]].append(
#         e.replace(".sort.mdup.bam", "").split("_")[1]
#     )
#     for e in sample_target_list
# ]




def get_mem_mb_heavy(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [64, 128, 256]
    return mem_avail[attempt - 1] * 1000


def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [2, 4, 8, 16, 32, 64]
    return mem_avail[attempt - 1] * 1000


rule all:
    input:
        [expand(
            "{results_folder}/GENOTYPING_OTF/{sample}/{cell_id}.vcf.gz",
            results_folder=results_folder,
            sample=sample,
            cell_id=cell_dict[sample],
        ) for sample in samples],
        expand(
            "{results_folder}/COVERAGE/{sample}/MERGE/merge_coverage.txt",
            results_folder=results_folder,
            sample=samples,
        ),
        
        
        # [
        #     expand(
        #         "{folder}/GENOTYPING_OTF/{sample_target}_{cell_id}.vcf.gz",
        #         folder=folder,
        #         sample=sample_list,
        #         sample_target=sample_target,
        #         cell_id=sample_target_dict[sample_target],
        #     )
        #     for sample_target in sample_target_dict
        # ],
        # expand("{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz", results_folder=results_folder),
        # expand("{folder}/COVERAGE_POOL1/{cell_id}.txt", folder=folder, cell_id=cell_dict),
        # [
        #     expand("{folder}/COVERAGE_POOL1/{sample_target}_{cell_id}.txt", folder=folder, sample_target=sample_target, cell_id=sample_target_dict[sample_target])
        #     for sample_target in sample_target_dict
        # ]
        # expand("{folder}/BCFTOOLS_CLEAN/merge.vcf.gz", folder=folder)
        # expand("{folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{sample}.txt.gz", folder=folder, sample=sample_list),
        # expand("{folder}/BCFTOOLS_OTF/{sample}.vcf.gz", folder=folder, sample=sample_list),
        # expand("{folder}/BCFTOOLS_OTF/{sample}.vcf.gz.tbi", folder=folder, sample=sample_list)
        # [
        #     expand(
        #     "{folder}/ISEC/{sample_target}_{cell_id}/{sample}/0001.vcf", folder=folder, sample=sample_list, sample_target=sample_target, cell_id=sample_target_dict[sample_target]
        #     ) for sample_target in sample_target_dict
        # ]
        # [
        #     expand(
        #     "{folder}/ANALYSIS/{sample_target}_{cell_id}_detailed.tsv", folder=folder, sample_target=sample_target, cell_id=sample_target_dict[sample_target]
        #     ) for sample_target in sample_target_dict
        # ]


rule create_single_1kgp_file_rare_het_intermediate:
    input:
        expand("{onekgp_folder}/{file}", onekgp_folder=config["onekgp_vcf_folder"], file=onekgp_filelist)
    output:
        "{onekgp_folder}/CUSTOM_VCF/1000G_rare005.vcf.gz"
    conda:
        "envs/snp_genotyping.yaml"
    threads: 10
    resources:
        mem_mb=128000,
        time="30:00:00",
    params:
        maf = 0.05
    shell:
        """
        bcftools concat {input}  | \
        bcftools view --type snps -i "INFO/AC>0 & INFO/AF<{params.maf}" | \
        bcftools sort -T {config[tmp_data_folder]} | \
        bcftools norm -d none --threads {threads} -Oz -o {output}
        """

rule create_single_1kgp_file_rare_het_decompose:
    input:
        "{onekgp_folder}/CUSTOM_VCF/1000G_rare005.vcf.gz"
    output:
        "{onekgp_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz"
    conda:
        "envs/snp_genotyping.yaml"
    threads: 1
    resources:
        mem_mb=128000,
        time="30:00:00",
    params:
        maf = 0.05
    shell:
        """
        vt decompose {input} | \
        bcftools annotate -x "^INFO/AC,INFO/AF" -Oz -o {output}
        """ 


rule retrieve_rare_variants_bcftools:
    input:
        onekgp="{onekgp_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
        onekgp_index="{onekgp_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
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
        import os,sys
        skip = len([e.decode("ISO-8859–1").split('\n') for e in gzip.open("{}".format(input[0], "rb")) if e.decode("ISO-8859–1").split('\n')[0].startswith("##")])
        sample = os.path.basename(input[0]).replace(".vcf.gz", "")
        df = pd.read_csv(input[0], compression="gzip", skiprows=skip, sep="\t")
        df["#CHROM"] = df["#CHROM"].str.replace("chr", "")
        df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str) + ":" + df["REF"] + ":" + df["ALT"]
        df["AC"] = df["INFO"].apply(lambda r: r.split(";")[0].replace("AC=", ""))
        df["AF"] = df["INFO"].apply(lambda r: r.split(";")[1].replace("AF=", ""))
        df["SAMPLE"] = wildcards.onekgp_sample
        df = df[["ID", "AC", "AF", "SAMPLE"]]
        df.to_csv(output[0], compression='gzip', sep='\t', index=False)
        # 'bcftools view --type snps {input} | cut -f 3,8 | sed "s/$/\\t{wildcards.sample}/g;s/AC=//g;s/AF=//g" |  grep -v "^#" | tail -n+2 | tr ";" "\\t" | bgzip > {output}'


rule merge_concat_bcftools_tab:
    input:
        expand(
            "{results_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt.gz",
            results_folder=results_folder,
            onekgp_sample=sample_list,
        ),
    output:
        "{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
    log:
        "{results_folder}/log/BCFTOOLS_CLEAN/{sample}/merge.log",
    conda:
        "mosaicatcher_env"
    threads: 1
    resources:
        mem_mb=16000,
    script:
        "scripts/merge_concat_bcftools.py"


rule merge_concat_bcftools_vcf:
    input:
        vcf=expand(
            "{results_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz", results_folder=results_folder, onekgp_sample=sample_list
        ),
        vcf_index=expand(
            "{results_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz.tbi", results_folder=results_folder, onekgp_sample=sample_list
        ),
    output:
        "{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz",
    log:
        "{results_folder}/log/BCFTOOLS_CONCAT/{sample}/merge.log",
    conda:
        "envs/snp_genotyping.yaml"
    threads: 1
    resources:
        mem_mb=16000,
    shell:
        "bcftools merge {input.vcf} | bcftools sort | bcftools norm -d none | cut -f -8 | bgzip > {output} 2> {log}"


rule samtools_coverage:
    input:
        bam=lambda wc: expand("{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam", bam_folder=bam_folder, sample=wc.sample),
        bai=lambda wc: expand("{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam.bai", bam_folder=bam_folder, sample=wc.sample),
    output:
        coverage="{results_folder}/COVERAGE/{sample}/{cell_id}.txt",
    conda:
        "envs/snp_genotyping.yaml"
    resources:
        mem_mb=2000,
        time="01:00:00",
    shell:
        "samtools coverage {input.bam} -o {output.coverage}"

rule merge_coverage:
    input:
        lambda wc: expand(
            "{results_folder}/COVERAGE/{sample}/{cell_id}.txt", 
            results_folder=results_folder, sample=wc.sample, cell_id=cell_dict[wc.sample]
            )
    output:
        "{results_folder}/COVERAGE/{sample}/MERGE/merge_coverage.txt"
    log:
        "{results_folder}/log/COVERAGE/{sample}/merge_coverage.log"
    run:    
        import pandas as pd
        l = list()
        chroms = ["chr" + str(e) for e in list(range(1,23))] + ["chrX"]
        # print(list(input))
        for file in list(input):
            # print(file)
            df = pd.read_csv(file, sep="\t")
            df = df.loc[df["#rname"].isin(chroms)]
            df["QUERY_SAMPLE_CELL"] = file.split("/")[-1].replace(".txt", "")
            df["QUERY_SAMPLE"] = file.split("/")[-2]
            # print(df)
            l.append(df)
        coverage_df = pd.concat(l)
        coverage_df_summary = coverage_df.groupby(["QUERY_SAMPLE_CELL"])["meandepth"].mean().reset_index()
        coverage_df_summary.to_csv(output[0], index=False, sep="\t")
            



# NOTE: use custom fasta that covers contigs present in all libraries from pseudopool
# fasta location: /g/korbel2/weber/MosaiCatcher_files/EXTERNAL_DATA/refgenomes_human_local/hg38.complete.fa
rule regenotype_SNVs:
    input:
        bam=lambda wc: expand("{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam", bam_folder=bam_folder, sample=wc.sample),
        bai=lambda wc: expand("{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam.bai", bam_folder=bam_folder, sample=wc.sample),
        sites="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz",
        sites_index="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz.tbi",
        fasta=config["fasta_ref"],
        fasta_index="{fasta_ref}.fai".format(fasta_ref=config["fasta_ref"]),
        merge_table="{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
    output:
        vcf="{results_folder}/GENOTYPING_OTF/{sample}/{cell_id}.vcf.gz",
    log:
        vcf="{results_folder}/log/GENOTYPING/{sample}/{cell_id}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="20:00:00",
    threads: 32
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        """
        /g/korbel2/weber/Gits/freebayes/scripts/freebayes-parallel \
            <(/g/korbel2/weber/Gits/freebayes/scripts/fasta_generate_regions.py {input.fasta} 100000) {threads} \
            -f {input.fasta} \
            -@ {input.sites} \
            --only-use-input-alleles {input.bam} \
            --genotype-qualities \
            --exclude-unobserved-genotypes \
        | bcftools view \
            --threads {threads} \
            --exclude-uncalled \
            --types snps \
            --genotype het \
            -Oz -o {output.vcf} 2> {log}
        """
            # --include "QUAL>=10" \


# rule analyse_isec:
#     input:
#         vcf=[
#             expand(
#                 "{folder}/GENOTYPING_OTF/{sample_target}_{cell_id}.vcf.gz",
#                 folder=folder,
#                 sample_target=sample_target,
#                 cell_id=sample_target_dict[sample_target],
#             )
#             for sample_target in sample_target_dict
#         ],
#         ref="{folder}/BCFTOOLS_CONCAT_TAB/merge.txt.gz",
#     output:
#         detailed="{folder}/ANALYSIS/analysis_detailed.tsv",
#         summary="{folder}/ANALYSIS/analysis_summary.tsv",
#     conda:
#         "mosaicatcher_env"
#     threads: 64
#     resources:
#         mem_mb=get_mem_mb_heavy,
#     script:
#         "scripts/gather_isec.py"



rule index_vcf:
    input:
        "{file}.vcf.gz",
    output:
        "{file}.vcf.gz.tbi",
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    envmodules:
        "tabix/0.2.6-GCCcore-11.2.0",
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
    envmodules:
        "tabix/0.2.6-GCCcore-11.2.0",
    shell:
        "bgzip {input}"
