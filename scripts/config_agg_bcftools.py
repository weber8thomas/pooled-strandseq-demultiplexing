with open(snakemake.output[0], 'w') as o:
    for v in list(snakemake.input.vcf):
        v_chrom = v.split("/")[-2]
        if snakemake.wildcards.chrom == v_chrom:
            o.write(v + '\n')