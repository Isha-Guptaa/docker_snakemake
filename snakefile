# List of samples
samples = ["father", "mother", "proband"]

# Rule to generate mapped BAM files
rule all:
    input:
        expand("/data/calls/all.avf", sample=samples)

# Rule to map reads using BWA
rule bwa_map:
    input:
        genome="/data/genome.fa",
        r1="/data/samples/{sample}_R1.fq",
        r2="/data/samples/{sample}_R2.fq"
    output:
        "/data/mapped/{sample}.bam"
    container:
        "my-snakemake-image"
    shell:
        "bwa mem {input.genome} {input.r1} {input.r2} | samtools view -Shb -o {output}"

# Rule to sort BAM files
rule sort:
    input:
        "/data/mapped/{sample}.bam"
    output:
        "/data/sorted_reads/{sample}.sorted.bam"
    container:
        "my-snakemake-image"
    shell:
        "samtools sort -o {output} {input}"

# Rule to index sorted BAM files
rule samtools_index:
    input:
        "/data/sorted_reads/{sample}.sorted.bam"
    output:
        "/data/sorted_reads/{sample}.bam.bai"
    container:
        "my-snakemake-image"
    shell:
        "samtools index {input} > {output}"

# Rule to call genomic variants
rule genomic_variants:
    input:
        genome="/data/genome.fa",
        bam=expand("/data/sorted_reads/{sample}.sorted.bam", sample=samples),
        bai=expand("/data/sorted_reads/{sample}.bam.bai", sample=samples)
    output:
        "/data/calls/all.avf"
    container:
        "my-snakemake-image"
    shell:
        "bcftools mpileup -f {input.genome} {input.bam} | bcftools call -mv - > {output}"
