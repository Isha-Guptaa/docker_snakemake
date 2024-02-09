configfile: "data/conda.yml"

# List of samples
samples = ["father","mother","proband"]

# Rule to generate mapped BAM files
rule all:
    input:
        expand("mapped/{sample}.bam", sample=samples)

# Rule to perform BWA mapping
rule bwa_map:
    input:
        genome="data/genome.fa",
        r1="data/samples/{sample}_R1.fq",
        r2="data/samples/{sample}_R2.fq"
    output:
        "mapped/{sample}.bam"
    conda:
        "data/conda.yml"
    shell:
        "bwa mem {input.genome} {input.r1} {input.r2} | samtools view -Shb -o {output}"

# Rule to sort BAM files
rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        "sorted_reads/{sample}.sorted.bam"
    conda:
        "data/conda.yml"
    shell:
        "samtools sort -o {output} {input}"

# Rule to index sorted BAM files
rule samtools_index:
    input:
        "sorted_reads/{sample}.sorted.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda:
        "data/conda.yml"
    shell:
        "samtools index {input} > {output}"

# Rule to call genomic variants
rule genomic_variants:
    input:
        genome="data/genome.fa",
        bam=expand("sorted_reads/{sample}.sorted.bam", sample=samples)
    output:
        "calls/all.vcf"
    conda:
        "data/conda.yml"
    shell:
        "bcftools mpileup -f {input.genome} {input.bam} | "
        "bcftools call -mv - > {output}"