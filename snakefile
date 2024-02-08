rule bwa_map:
    input:
        "/data/genome.fa",
        "/data/father_R1.fq",
        "/data/father_R2.fq"
    output:
        "/data/mapped/father.bam"
    container:
        "my-snakemake-image"
    shell:
        "bwa mem {input} | samtools view -Shb -o {output}"

rule sort:
    input:
        "/data/mapped/father.bam"
    output:
        "/data/sorted_reads/father.sorted.bam"
    container:
        "my-snakemake-image"
    shell:
        "samtools sort -o {output} {input}"

rule samtools_index:
    input:
        "/data/sorted_reads/father.sorted.bam"
    output:
        "/data/sorted_reads/father.bam.bai"
    container:
        "my-snakemake-image"
    shell:
        "samtools index {input} > {output}"

rule genomic_variants:
    input:
        genome="/data/genome.fa",
        bam=expand("/data/sorted_reads/father.sorted.bam"),
        bai=expand("/data/sorted_reads/father.bam.bai")
    output:
        "/data/calls/all.avf"
    container:
        "my-snakemake-image"
    shell:
        "bcftools mpileup -f {input.genome} {input.bam} | "
        "bcftools call -mv - > {output}"

# Rule to gather all final outputs
rule gather_results:
    input:
        expand("/data/calls/all.avf"),
        expand("/data/plots/quals.svg")
    output:
        "/data/results/all_results.txt"
    shell:
        "echo 'All results are gathered.' > {output}"

# Specify the final target rule
rule all:
    input:
        ["/data/results/all_results.txt"]

