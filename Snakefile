import subprocess
import os
import collections

configfile: "config.yaml"

# chroms = config["chroms"]
chroms = ["chr" + str(e) for e in list(range(1,23))]
print(chroms)

# folder = "/g/korbel/weber/MosaiCatcher_files/EXTERNAL_DATA/snv_sites_to_genotype/1000G_SNV_with_GT"
# folder = "/scratch/tweber/DATA/NYGC_1KGP"
# folder = "/g/korbel/weber/MosaiCatcher_files/EXTERNAL_DATA/snv_sites_to_genotype/1000G_SNV_with_GT"
folder = "/scratch/tweber/DATA/1000G_SNV_with_GT"



filelist = [e for e in os.listdir(folder) if e.endswith(".vcf.gz")]
header = subprocess.Popen("zcat {folder}/{file} | head -n 5000 | grep -P  '^#CHROM'".format(folder=folder, file=filelist[0]), shell=True, stdout=subprocess.PIPE)
samples = header.communicate()[0].decode("utf-8").strip().split("\t")[9:]
dict_samples = {e:j for j,e in enumerate(samples)}
hgsvc_samples = ["HG00155", "HG00186", "HG00233", "HG00267", "HG00423", "HG00480", "HG00673", "HG00735", "HG00766", "HG01071", "HG01106", "HG01255", "HG01273", "HG01493", "HG01609", "HG01887", "HG01928", "HG01952", "HG01978", "HG02068", "HG02132", "HG02135", "HG02148", "HG02602", "HG02698", "HG02723", "HG02975", "HG03035", "HG03056", "HG03200", "HG03456", "HG03519", "HG03522", "HG03576", "HG03654", "HG03719", "HG03816", "HG03866", "HG03942", "HG03968", "HG04115", "HG04157", "HG04204", "NA10837", "NA10857", "NA18537", "NA18940", "NA18965", "NA19100", "NA19139", "NA19440", "NA19443", "NA19721", "NA19724", "NA19727", "NA19730", "NA19902", "NA20128", "NA20536", "NA20846",]
pool_samples = ["HG00268", "HG00512", "HG00514", "HG00731", "HG01352", "HG02059", "HG02818", "NA12878", "NA19239", "NA19240",]


sample_list = list(set([e for e in samples] + hgsvc_samples + pool_samples))[:10]
folder_target = "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/all"

sample_target_list = [e for e in os.listdir(folder_target) if e.endswith(".bam")] 
print(sample_target_list)
sample_target_dict = collections.defaultdict(list)
[sample_target_dict[e.replace('.sort.mdup.bam', '').split('_')[0]].append(e.replace('.sort.mdup.bam', '').split('_')[1]) for e in sample_target_list]
print(sample_target_dict)

def get_mem_mb_heavy(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [32,64,128,256]
    return mem_avail[attempt - 1] * 1000

def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [2,4,8,16,32,64]
    return mem_avail[attempt - 1] * 1000


rule all:
    input:
            # "{folder}/PANDAS/{sample}.tsv.gz", folder=folder, chrom=chroms, sample=list(set(samples[:100] + ["HG00268","HG00512","HG00514","HG00731","HG01352","HG02059","HG02818","NA12878","NA19239","NA19240"]))
            # "{folder}/PANDAS/{sample}.tsv.gz", folder=folder, chrom=["1"], sample=samples[:2]
            # expand("{folder}/AGGR/1000G_rare_het.vcf.gz", folder=folder)
            # "{folder}/AGGR.tsv.gz", folder=folder, chrom=["chr21"],
            # "{folder}/GENOTYPING/HG00512.vcf.gz", folder=folder, chrom=["chr1"]
            # "{folder}/ISEC/HG00512_{sample}", folder=folder, sample=sample_list, chrom=chroms
        # [
        #     expand(
        #     "{folder}/ANALYSIS/{sample_target}_{cell_id}_detailed.tsv", folder=folder, sample_target=sample_target, cell_id=sample_target_dict[sample_target]
        #     ) for sample_target in sample_target_dict
        # ]
        [
            expand(
            "{folder}/ISEC/{sample_target}_{cell_id}/{sample}/0002.vcf", folder=folder, sample=sample_list, sample_target=sample_target, cell_id=sample_target_dict[sample_target]
            ) for sample_target in sample_target_dict
        ]
        # expand("{folder}/BCFTOOLS/CHRFREE/{sample}.vcf.gz.tbi", folder=folder, sample=sample_list)


# rule retrieve_rare_variants_bcftools:
#     input:
#         # onekgp= "{folder}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
#         # onekgp= "{folder}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.recalibrated_variants.vcf.gz"
#         onekgp= expand("/g/korbel/weber/MosaiCatcher_files/EXTERNAL_DATA/snv_sites_to_genotype/1000G_SNV_with_GT/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz", chrom=chroms)
#     output:
#         "{folder}/BCFTOOLS/{sample}.vcf.gz"
#     conda:
#         "snp.yaml"
#     envmodules:
#         "BCFtools/1.14-GCC-11.2.0"
#     threads: 1
#     resources:
#         mem_mb=get_mem_mb_heavy,
#         time="01:00:00",
#     params:
#         index=lambda wc: dict_samples[wc.sample],
#         col_index=lambda wc: dict_samples[wc.sample] + 1 + 9, # +1 for python index +9 for the nine first cols
#         af=0.01
#     shell:
#         # "bcftools view -i 'INFO/AC>0 & INFO/AF<{params.af} & GT[{params.index}]=\"het\"' {input} | cut -f 1-8 | sed '/^#/! s/$/;SAMPLE={wildcards.sample}/' | bgzip > {output}"
#         "bcftools view --threads {threads} -i 'INFO/AC>0 & INFO/AF<{params.af} & GT[{params.index}]=\"het\"' {input} | cut -f 1-9,{params.col_index} | bgzip > {output}"

# localrules: config_agg_bcftools


# rule config_agg_bcftools:
#     input:
#         vcf = expand("{folder}/BCFTOOLS/{sample}.vcf.gz", folder=folder, sample=sample_list),
#         index = expand("{folder}/BCFTOOLS/{sample}.vcf.gz.tbi", folder=folder, sample=sample_list)
#     output:
#         "{folder}/CONFIG/config.txt"
#     resources:
#         time="01:00:00",
#     script:
#         "config_agg_bcftools.py"
#         # "ls {input.vcf} | grep {wildcards.chrom} > {output}"

# rule merge_bcftools:
#     input:
#         file_list = "{folder}/CONFIG/config.txt"
#     output:
#         "{folder}/AGGR/1000G_rare_het.vcf.gz"
#     log:
#         "{folder}/log/AGGR.log"
#     conda:
#         "snp.yaml"
#     envmodules:
#         "BCFtools/1.14-GCC-11.2.0"
#     threads: 1
#     resources:
#         mem_mb=128000,
#         time="10:00:00",
#     shell:
#         "bcftools merge -l {input.file_list} | bcftools sort | bcftools norm -d none | cut -f -8 | bgzip > {output} 2> {log}"


rule concat_bcftools:
    input:
        file_list = expand("{{folder}}/BCFTOOLS/{chrom}/{{sample}}.vcf.gz", chrom=chroms)
    output:
        "{folder}/BCFTOOLS/CHRFREE/{sample}.vcf.gz"
    log:
        "{folder}/log/BCFTOOLS/CHRFREE/{sample}.vcf.gz"
    conda:
        "snp.yaml"
    envmodules:
        "BCFtools/1.14-GCC-11.2.0"
    threads: 1
    resources:
        mem_mb=2000,
        time="10:00:00",
    shell:
        "bcftools concat --threads {threads}  {input.file_list} | bcftools sort -T /scratch/tweber/TMP -Oz -o {output} 2> {log}"


rule regenotype_SNVs:
    """
    rule fct:
    input:
    output:
    """
    input:
        bam="/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/all/{sample_target}_{cell_id}.sort.mdup.bam",
        bai="/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/all/{sample_target}_{cell_id}.sort.mdup.bam.bai",
        sites="{folder}/AGGR/1000G_rare_het.vcf.gz",
        sites_index="{folder}/AGGR/1000G_rare_het.vcf.gz.tbi",
        fasta="/g/korbel2/weber/MosaiCatcher_files/EXTERNAL_DATA/refgenomes_human/hg38.fa",
        fasta_index="/g/korbel2/weber/MosaiCatcher_files/EXTERNAL_DATA/refgenomes_human/hg38.fa.fai",
    output:
        vcf="{folder}/GENOTYPING/{sample_target}_{cell_id}.vcf",
    log:
        vcf="{folder}/log/GENOTYPING/{sample_target}_{cell_id}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    conda:
        "snp.yaml"
    shell:
        """
        freebayes \
            -f {input.fasta} \
            -@ {input.sites} \
            --only-use-input-alleles {input.bam} \
            --genotype-qualities \
        | bcftools view \
            --types snps \
            --genotype het \
            -Ov -o {output.vcf}
        """
            

rule isec:
    input:
        vcf="{folder}/GENOTYPING/{sample_target}_{cell_id}.vcf.gz",
        vcf_index="{folder}/GENOTYPING/{sample_target}_{cell_id}.vcf.gz.tbi",
        bcftools="{folder}/BCFTOOLS/CHRFREE/{sample}.vcf.gz",
        bcftools_index="{folder}/BCFTOOLS/CHRFREE/{sample}.vcf.gz.tbi"
    output:
        isect=directory("{folder}/ISEC/{sample_target}_{cell_id}/{sample}"),
        isect_file="{folder}/ISEC/{sample_target}_{cell_id}/{sample}/0002.vcf",
        # t=touch("{folder}/ISECT/{sample_target}_{sample}.ok")
    conda:
        "snp.yaml"
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    shell:
        "bcftools isec -p {output.isect} {input.vcf} {input.bcftools}"

rule analyse_isec:
    input:  
        [
            expand("{folder}/ISEC/{sample_target}_{cell_id}/{sample}/0002.vcf", folder=folder, chrom=chroms, sample=sample_list, sample_target=sample_target, cell_id=sample_target_dict[sample_target])
            for sample_target in sample_target_dict
        ]
    output:
        detailed="{folder}/ANALYSIS/{sample_target}_{cell_id}_detailed.tsv",
        summary="{folder}/ANALYSIS/{sample_target}_{cell_id}_summary.tsv"
    threads: 64
    resources:
        mem_mb=64000,
        time="01:00:00",
    script:
        "gather_isec.py"


rule index_vcf:
    input:
        "{file}.vcf.gz"
    output:
        "{file}.vcf.gz.tbi"
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    envmodules:
        "tabix/0.2.6-GCCcore-11.2.0"
    shell:
        "tabix -p vcf {input}"

rule compress_vcf:
    input:
        "{file}.vcf"
    output:
        "{file}.vcf.gz"
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    envmodules:
        "tabix/0.2.6-GCCcore-11.2.0"
    shell:
        "bgzip {input}"
