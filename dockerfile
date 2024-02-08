# Use a more recent Ubuntu base image
FROM ubuntu:focal

# Set noninteractive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    bcftools \
    python3 \
    python3-pip

# Install specific versions of pulp and tabulate
RUN pip3 install pulp==2.4 tabulate==0.8.9

# Install specific version of Snakemake
RUN pip3 install snakemake==6.10.0

# Set the working directory
WORKDIR /data

# Copy the Snakemake workflow files into the container
COPY snakefile /data/snakefile

# Define the entry point for Snakemake
ENTRYPOINT ["snakemake", "--cores", "1"]