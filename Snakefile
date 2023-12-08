import subprocess
import os
import collections


# CONFIG PART
configfile: "config/config.yaml"


results_folder = config["results_folder"]
onekgp_vcf_folder = config["onekgp_vcf_folder"]
bam_folder = config["bam_folder"]
sample_list = config["samples_pooled"]
# Rename 1KG samples from GM to NA
# sample_list = {k: [e.replace("GM", "NA") for e in v] for k, v in sample_list.items()}
# print(sample_list)

# list 1kgp files
onekgp_filelist = [e for e in os.listdir(onekgp_vcf_folder) if e.endswith(".vcf.gz")]
print(onekgp_filelist)
chroms_list = [e.replace(".vcf.gz", "") for e in onekgp_filelist]
print(chroms_list)

print("{folder}/{file}".format(folder=onekgp_vcf_folder, file=onekgp_filelist[0]))


# Command to extract the header line containing sample names
cmd = f"zcat {onekgp_vcf_folder}/{onekgp_filelist[0]} | grep -m 1 '^#CHROM'"

# Run the command
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()

# Check for errors
if process.returncode != 0:
    sys.exit("Error:", stderr.decode())
else:
    # Decode and split the header line to get sample names
    header_line = stdout.decode().strip()
    samples = header_line.split("\t")[9:]  # Sample names start from the 10th column

    # Print sample names
    print("Sample names:", samples)

    # Creating a dictionary with index numbers
    sample_index_dict = {sample: idx for idx, sample in enumerate(samples, start=9)}
    print("Sample index dictionary:", sample_index_dict)


# # Retrieve header
# header = subprocess.Popen(
#     "zcat {folder}/{file} | head -n 100 | grep -P  '^#CHROM'".format(
#         folder=onekgp_vcf_folder, file=onekgp_filelist[0]
#     ),
#     shell=True,
#     stdout=subprocess.PIPE,
# )

# print(header)

# # List of 1kgp samples & dict with index nb
# samples = header.communicate()[0].decode("utf-8").strip().split("\\t")[9:]
dict_samples = {e: j for j, e in enumerate(samples)}

# Retrieve list of cells for each Pool sample to be genotyped
cell_dict = collections.defaultdict(list)
for sample in os.listdir(bam_folder):
    if os.path.isdir(bam_folder + "/" + sample) and sample in sample_list.keys():
        for e in os.listdir(bam_folder + "/" + sample + "/selected"):
            if e.endswith(".bam"):
                cell_dict[sample].append(e.replace(".sort.mdup.bam", ""))

# List of pool samples to process
samples = list(cell_dict.keys())
# print(samples)
# print([(dict_samples[e]) for s,e_l in sample_list.items() for e in e_l])
print(sample_list)



# sample_list["random"] = []
# exit()


def get_mem_mb_heavy(wildcards, attempt):
    """
    To adjust resources in the rules
    """
    mem_avail = [64, 128, 256]
    return mem_avail[attempt - 1] * 1000


def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    """
    mem_avail = [2, 4, 8, 16, 32, 64]
    return mem_avail[attempt - 1] * 1000


# print([expand("{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt",
#         onekgp_vcf_folder=onekgp_vcf_folder,
#         onekgp_sample=sample_list[k]
#         ) for k in samples] +
#         [expand("{results_folder}/GENOTYPING_OTF/{sample}/{cell_id}.vcf.gz", results_folder=results_folder, sample=k, cell_id=cell_dict[k]) for k in samples]
# )


print(samples)
print(sample_list)
print([expand("{results_folder}/COMPARING_STRANDSEQ_TO_1KG/{strandseq_sample}_{onekgp_sample}.tsv.gz", results_folder=results_folder, strandseq_sample=k, onekgp_sample=sample_list[k]) for k in samples])
print(dict_samples)


rule all:
    input:
        # expand(
        #     "{results_folder}/ANALYSIS/{sample}.xlsx",
        #     results_folder=results_folder,
        #     sample=samples,
        # ),
        # [expand("{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt",
        # onekgp_vcf_folder=onekgp_vcf_folder,
        # onekgp_sample=sample_list[k]
        # ) for k in samples] +
        # [expand("{bam_folder}/{sample}/haplotag/bam_HP_filtered/{cell_id}.bam.htg.HP.bai", bam_folder=bam_folder, sample=k, cell_id=cell_dict[k])[:10] for k in samples]
        # [expand("{results_folder}/GENOTYPING_OTF/{strandseq_sample}.sorted.vcf.gz", results_folder=results_folder, strandseq_sample=k) for k in samples]
        [expand("{results_folder}/SNPs/{strandseq_sample}/{cell_id}.vcf.gz", results_folder=results_folder, strandseq_sample=k, cell_id=cell_dict[k]) for k in samples],
        [expand("{results_folder}/COMPARING_STRANDSEQ_TO_1KG/{strandseq_sample}_{onekgp_sample}.tsv.gz", results_folder=results_folder, strandseq_sample=k, onekgp_sample=sample_list[k]) for k in samples]


# rule create_single_1kgp_file_rare_het_intermediate:
# localrules:
#     create_single_1kgp_file_intermediate,


# rule create_single_1kgp_file_intermediate:
#     input:
#         expand(
#             "{onekgp_vcf_folder}/{file}",
#             onekgp_vcf_folder=config["onekgp_vcf_folder"],
#             file=onekgp_filelist,
#         ),
#     output:
#         # "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005.vcf.gz",
#         "{onekgp_vcf_folder}/CUSTOM_VCF/1000G.vcf.gz",
#     conda:
#         "envs/snp_genotyping.yaml"
#     threads: 10
#     resources:
#         mem_mb=128000,
#         time="30:00:00",
#     # params:
#     #     maf=0.05,
#     shell:
#         """
#         bcftools concat -a {input} | \
#         bcftools view --type snps -i "INFO/AC>0" | \
#         bcftools norm -d none --threads {threads} -Oz -o {output}
#         """
#         # """
#         # bcftools concat {input}  | \
#         # bcftools view --type snps -i "INFO/AC>0 & INFO/AF<{params.maf}" | \
#         # bcftools sort -T {config[tmp_data_folder]} | \
#         # bcftools norm -d none --threads {threads} -Oz -o {output}
#         # """


# # rule create_single_1kgp_file_rare_het_decompose:
# rule create_single_1kgp_file_decompose:
#     input:
#         # "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005.vcf.gz",
#         "{onekgp_vcf_folder}/CUSTOM_VCF/1000G.vcf.gz",
#     output:
#         # "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
#         "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_decompose.vcf.gz",
#     conda:
#         "envs/snp_genotyping.yaml"
#     threads: 1
#     resources:
#         mem_mb=128000,
#         time="30:00:00",
#     shell:
#         """
#         vt decompose {input} | \
#         bcftools annotate -x "^INFO/AC,INFO/AF" -Oz -o {output}
#         """


rule retrieve_variants_bcftools:
    input:
        fasta=expand(
            "{onekgp_vcf_folder}/{{chrom}}.vcf.gz",
            onekgp_vcf_folder=config["onekgp_vcf_folder"],
            # chrom=chroms_list,
        ),
        fasta_index=expand(
            "{onekgp_vcf_folder}/{{chrom}}.vcf.gz.tbi",
            onekgp_vcf_folder=config["onekgp_vcf_folder"],
            # chrom=chroms_list,
        ),
        # onekgp=ancient(
        #     expand(
        #         # "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
        #         "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_decompose.vcf.gz",
        #         onekgp_vcf_folder=onekgp_vcf_folder,
        #     )
        # ),
        # onekgp_index=ancient(
        #     expand(
        #         # "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_rare005_decompose.vcf.gz",
        #         "{onekgp_vcf_folder}/CUSTOM_VCF/1000G_decompose.vcf.gz",
        #         onekgp_vcf_folder=onekgp_vcf_folder,
        #     )
        # ),
    output:
        "{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}/{chrom}.txt.raw",
        # "{onekgp_vcf_folder}/s/{onekgp_sample}/{chrom}.vcf.gz",
    conda:
        "envs/snp_genotyping.yaml"
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    params:
        index=lambda wc: dict_samples[wc.onekgp_sample],
        col_index=lambda wc: dict_samples[wc.onekgp_sample] + 1 + 9,  # +1 for python index +9 for the nine first cols
        # af=0.05,
    shell:
        "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\n' {input.fasta}  -s {wildcards.onekgp_sample}  > {output}"
        # "bcftools view --threads {threads} -i 'INFO/AC>0' --type snps {input.onekgp} | cut -f 1-9,{params.col_index} | bgzip > {output}"
        # "bcftools view --threads {threads} -i 'INFO/AC>0 & INFO/AF<{params.af}' --type snps {input.onekgp} | cut -f 1-9,{params.col_index} | bgzip > {output}"
        # "bcftools view --threads {threads} -i 'INFO/AC>0 & INFO/AF<{params.af} & GT[{params.index}]=\"het\"' --type snps {input.onekgp} | cut -f 1-9,{params.col_index} | bgzip > {output}"

rule filter_variants_sample:
    input:
        "{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}/{chrom}.txt.raw",   
    output:
        "{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}/{chrom}.txt.filtered",   
    shell:
        "awk 'length($3) == 1 && length($4) == 1 && $5 != \"0|0\"' {input} > {output}"



rule aggregate_filter_variants_sample:
    input: 
        expand("{{onekgp_vcf_folder}}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{{onekgp_sample}}/{chrom}.txt.filtered", chrom=chroms_list)
    output:
        "{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt.gz"
    shell:
        """
        cat {input} |  awk -F'\\t' 'BEGIN {{OFS="\\t"}} {{print $1"_"$2"_"$3"_"$4, $0}}'|  sort -k1,1 | gzip > {output}
        """






# rule concat_filter_bcftools:
#     input:
#         "{onekgp_vcf_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz",
#     output:
#         "{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt.gz",
#     log:
#         "{onekgp_vcf_folder}/log/BCFTOOLS_CLEAN/{onekgp_sample}.log",
#     threads: 1
#     resources:
#         mem_mb=2000,
#     run:
#         import pandas as pd
#         import gzip
#         import os, sys

#         skip = len(
#             [
#                 e.decode("ISO-8859–1").split("\n")
#                 for e in gzip.open("{}".format(input[0], "rb"))
#                 if e.decode("ISO-8859–1").split("\n")[0].startswith("##")
#             ]
#         )
#         sample = os.path.basename(input[0]).replace(".vcf.gz", "")
#         df = pd.read_csv(input[0], compression="gzip", skiprows=skip, sep="\\t")
#         df["#CHROM"] = df["#CHROM"].str.replace("chr", "")
#         df["ID"] = (
#             df["#CHROM"].astype(str)
#             + ":"
#             + df["POS"].astype(str)
#             + ":"
#             + df["REF"]
#             + ":"
#             + df["ALT"]
#         )
#         df["AC"] = df["INFO"].apply(lambda r: r.split(";")[0].replace("AC=", ""))
#         df["AF"] = df["INFO"].apply(lambda r: r.split(";")[1].replace("AF=", ""))
#         df["SAMPLE"] = wildcards.onekgp_sample
#         df = df[["ID", "AC", "AF", "SAMPLE"]]
#         df.to_csv(output[0], compression="gzip", sep="\\t", index=False)
#         # 'bcftools view --type snps {input} | cut -f 3,8 | sed "s/$/\\\t{wildcards.sample}/g;s/AC=//g;s/AF=//g" |  grep -v "^#" | tail -n+2 | tr ";" "\\\t" | bgzip > {output}'



# rule merge_concat_bcftools_tab:
#     input:
#         lambda wc: expand(
#             "{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{onekgp_sample}.txt.gz",
#             onekgp_vcf_folder=onekgp_vcf_folder,
#             onekgp_sample=sample_list[wc.sample],
#         ),
#     output:
#         "{onekgp_vcf_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
#     log:
#         "{onekgp_vcf_folder}/log/BCFTOOLS_CLEAN/{sample}/merge.log",
#     conda:
#         "envs/python_env.yaml"
#     threads: 1
#     resources:
#         mem_mb=16000,
#     script:
#         "scripts/merge_concat_bcftools.py"


# rule merge_concat_bcftools_vcf:
#     input:
#         vcf=lambda wc: expand(
#             "{onekgp_vcf_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz",
#             onekgp_vcf_folder=onekgp_vcf_folder,
#             onekgp_sample=sample_list[wc.sample],
#         ),
#         vcf_index=lambda wc: expand(
#             "{onekgp_vcf_folder}/BCFTOOLS_OTF/{onekgp_sample}.vcf.gz.tbi",
#             onekgp_vcf_folder=onekgp_vcf_folder,
#             onekgp_sample=sample_list[wc.sample],
#         ),
#     output:
#         "{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz",
#     log:
#         "{results_folder}/log/BCFTOOLS_CONCAT/{sample}/merge.log",
#     conda:
#         "envs/snp_genotyping.yaml"
#     threads: 1
#     resources:
#         mem_mb=16000,
#     shell:
#         "bcftools merge {input.vcf} | bcftools sort | bcftools norm -d none | cut -f -8 | bgzip > {output} 2> {log}"


# rule samtools_coverage:
#     input:
#         bam=lambda wc: expand(
#             "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam",
#             bam_folder=bam_folder,
#             sample=wc.sample,
#         ),
#         bai=lambda wc: expand(
#             "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam.bai",
#             bam_folder=bam_folder,
#             sample=wc.sample,
#         ),
#     output:
#         coverage="{results_folder}/COVERAGE/{sample}/{cell_id}.txt",
#     conda:
#         "envs/snp_genotyping.yaml"
#     resources:
#         mem_mb=2000,
#         time="01:00:00",
#     shell:
#         "samtools coverage {input.bam} -o {output.coverage}"


# rule merge_coverage:
#     input:
#         lambda wc: expand(
#             "{results_folder}/COVERAGE/{sample}/{cell_id}.txt",
#             results_folder=results_folder,
#             sample=wc.sample,
#             cell_id=cell_dict[wc.sample],
#         ),
#     output:
#         "{results_folder}/COVERAGE/{sample}/MERGE/merge_coverage.txt",
#     log:
#         "{results_folder}/log/COVERAGE/{sample}/merge_coverage.log",
#     run:
#         import pandas as pd

#         l = list()
#         chroms = ["chr" + str(e) for e in list(range(1, 23))] + ["chrX"]
#         for file in list(input):
#             df = pd.read_csv(file, sep="\\t")
#             df = df.loc[df["#rname"].isin(chroms)]
#             df["QUERY_SAMPLE_CELL"] = file.split("/")[-1].replace(".txt", "")
#             df["QUERY_SAMPLE"] = file.split("/")[-2]
#             l.append(df)
#         coverage_df = pd.concat(l)
#         coverage_df_summary = (
#             coverage_df.groupby(["QUERY_SAMPLE_CELL"])["meandepth"].mean().reset_index()
#         )
#         coverage_df_summary.to_csv(output[0], index=False, sep="\\t")


# NOTE: use custom fasta that covers contigs present in all libraries from pseudopool
# fasta location: /g/korbel2/weber/MosaiCatcher_files/EXTERNAL_DATA/refgenomes_human_local/hg38.complete.fa


rule filter_haplotagged_bam:
    input:
        bam="{bam_folder}/{strandseq_sample}/haplotag/bam/{cell_id}.bam.htg",
        bai="{bam_folder}/{strandseq_sample}/haplotag/bam/{cell_id}.bam.htg.bai",
    output:
        "{bam_folder}/{strandseq_sample}/haplotag/bam_HP_filtered/{cell_id}.bam.htg.HP"
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        """
        samtools view -h {input.bam} | awk '$1 ~ /^@/ || $0 ~ /HP:i:[12]/' | samtools view -b -o {output}
        """


rule index_bam_HP_filtered:
    input:
        "{bam_folder}/{strandseq_sample}/haplotag/bam_HP_filtered/{cell_id}.bam.htg.HP"
    output:
        "{bam_folder}/{strandseq_sample}/haplotag/bam_HP_filtered/{cell_id}.bam.htg.HP.bai"
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        "samtools index {input}"

rule call_SNVs:
    input:
        bam=lambda wc: expand(
            "{bam_folder}/{strandseq_sample}/selected/{{cell_id}}.sort.mdup.bam",
            bam_folder=bam_folder,
            strandseq_sample=wc.strandseq_sample,
        ),
        # bam=lambda wc: expand(
        #     "{bam_folder}/{strandseq_sample}/haplotag/bam_HP_filtered/{{cell_id}}.bam.htg.HP",
        #     bam_folder=bam_folder,
        #     sample=wc.strandseq_sample,
        # ),
        # bai=lambda wc: expand(
        #     "{bam_folder}/{strandseq_sample}/haplotag/bam_HP_filtered/{{cell_id}}.bam.htg.HP.bai",
        #     bam_folder=bam_folder,
        #     sample=wc.strandseq_sample,
        # ),
        # sites="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz",
        # sites_index="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz.tbi",
        # fasta=expand(
        #     "{onekgp_vcf_folder}/{{chrom}}.vcf.gz",
        #     onekgp_vcf_folder=config["onekgp_vcf_folder"],
        #     # chrom=chroms_list,
        # ),
        # fasta_index=expand(
        #     "{onekgp_vcf_folder}/{{chrom}}.vcf.gz.tbi",
        #     onekgp_vcf_folder=config["onekgp_vcf_folder"],
        #     # chrom=chroms_list,
        # ),
        fasta=config["fasta_ref"],
        fasta_index="{fasta_ref}.fai".format(fasta_ref=config["fasta_ref"]),
        # merge_table="{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
    output:
        vcf="{results_folder}/SNPs/{strandseq_sample}/{cell_id}.vcf.gz",
    # log:
    #     vcf="{results_folder}/log/SNPs/{strandseq_sample}/{cell_id}.log",
    resources:
        # partition="bigmem",
        mem_mb=get_mem_mb,
        time="20:00:00",
    threads: 1
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        """
        bcftools mpileup -Ou -f {input.fasta} {input.bam} | bcftools call -cv | bcftools view -Oz --types snps > {output.vcf}
        """




rule gather_SNVs_from_strandseq_sample:
    input:
        lambda wc: expand(
            "{{results_folder}}/SNPs/{{strandseq_sample}}/{cell_id}.vcf.gz",
            cell_id=cell_dict[wc.strandseq_sample],
        )
    output:
        "{results_folder}/SNPs/GATHER_{strandseq_sample}.tsv.gz.raw"
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        """
        for file in {input}; do
            cell_name=$(basename "$file" .txt)
            zcat "$file" | awk -v strandseq_sample="$cell_name" 'BEGIN {{OFS="\\t"}} !/^#/ {{print $0, strandseq_sample}}'
        done | gzip > {output}
        """

rule sort_gather_SNVs_from_strandseq_sample:
    input:
        "{results_folder}/SNPs/GATHER_{strandseq_sample}.tsv.gz.raw"
    output:
        "{results_folder}/SNPs/GATHER_{strandseq_sample}.sorted.tsv.gz"
    conda:
        "envs/snp_genotyping.yaml"
    shell:
        """
        zcat {input} | awk -F'\\t' 'BEGIN {{OFS="\\t"}} {{print $1"_"$2"_"$4"_"$5, $0}}' | sort -k1,1 |  gzip > {output}
        """





rule join_strandseq_and_1KG:
    input:
        vcf_strandseq="{results_folder}/SNPs/GATHER_{strandseq_sample}.sorted.tsv.gz",
        tab_onekgp_sample=expand("{onekgp_vcf_folder}/BCFTOOLS_CLEAN_OTF/SAMPLEWISE/{{onekgp_sample}}.txt.gz", onekgp_vcf_folder=onekgp_vcf_folder)

    output:
        "{results_folder}/COMPARING_STRANDSEQ_TO_1KG/{strandseq_sample}_{onekgp_sample}.tsv.gz",
    shell:
        """
        zcat {input.vcf_strandseq} | join -t $'\\t' - <(zcat {input.tab_onekgp_sample}) | gzip > {output}
        """




# rule regenotype_SNVs:
#     input:
#         bam=lambda wc: expand(
#             "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam",
#             bam_folder=bam_folder,
#             sample=wc.sample,
#         ),
#         bai=lambda wc: expand(
#             "{bam_folder}/{sample}/bam/{{cell_id}}.sort.mdup.bam.bai",
#             bam_folder=bam_folder,
#             sample=wc.sample,
#         ),
#         sites="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz",
#         sites_index="{results_folder}/BCFTOOLS_CONCAT_VCF/{sample}/merge.vcf.gz.tbi",
#         fasta=config["fasta_ref"],
#         fasta_index="{fasta_ref}.fai".format(fasta_ref=config["fasta_ref"]),
#         merge_table="{results_folder}/BCFTOOLS_CONCAT_TAB/{sample}/merge.txt.gz",
#     output:
#         vcf="{results_folder}/GENOTYPING_OTF/{sample}/{cell_id}.vcf.gz",
#     log:
#         vcf="{results_folder}/log/GENOTYPING/{sample}/{cell_id}.log",
#     resources:
#         # partition="bigmem",
#         mem_mb=get_mem_mb_heavy,
#         time="20:00:00",
#     threads: 1
#     conda:
#         "envs/snp_genotyping.yaml"
#     shell:
#         """
#         freebayes \
#             -f {input.fasta} \
#             {input.bam} \
#             --genotype-qualities \
#             --exclude-unobserved-genotypes \
#         | bcftools view \
#             --threads {threads} \
#             --exclude-uncalled \
#             --types snps \
#             -Oz -o {output.vcf} 2> {log}
#         """
#         # """
#         # /g/korbel2/weber/Gits/freebayes/scripts/freebayes-parallel \
#         #     <(/g/korbel2/weber/Gits/freebayes/scripts/fasta_generate_regions.py {input.fasta} 100000) {threads} \
#         #     -f {input.fasta} \
#         #     -@ {input.sites} \
#         #     --only-use-input-alleles {input.bam} \
#         #     --genotype-qualities \
#         #     --exclude-unobserved-genotypes \
#         # | bcftools view \
#         #     --threads {threads} \
#         #     --exclude-uncalled \
#         #     --types snps \
#         #     --genotype het \
#         #     -Oz -o {output.vcf} 2> {log}
#         # """


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
#         partition="bigmem",
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
