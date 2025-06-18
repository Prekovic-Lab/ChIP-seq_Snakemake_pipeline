# File: transform_files.smk

import os
from pathlib import Path

# Load configuration
configfile: "config/transform_files.yaml"

# Define wildcard constraints
wildcard_constraints:
    sample="[^/]+"

# Define the final output files based on the starting file type
if config["starting_file"] == "cram":
    final_output = expand("aligned_reads/{sample}.bam", sample=config["samples"]) + \
                   expand("aligned_reads/{sample}.bam.bai", sample=config["samples"])
elif config["starting_file"] == "bam":
    final_output = expand(os.path.join(config["storing_cram"], "{sample}.cram"), sample=config["samples"]) + \
                   expand(os.path.join(config["storing_cram"], "{sample}.cram.crai"), sample=config["samples"])
else:
    raise ValueError("starting_file in config must be either 'cram' or 'bam'")

# Define the final output files
rule all:
    input:
        final_output

# Convert CRAM to BAM
rule cram_to_bam:
    input:
        cram = os.path.join(config["cram_dir"], "{sample}.cram"),
        ref = config["reference_genome"]
    output:
        bam = temp("aligned_reads/{sample}.unsorted.bam")
    log:
        "logs/cram_to_bam/{sample}.log"
    conda:
        "envs/chip_pipe_env.yml"
    threads: config["threads"]
    shell:
        "samtools view -@ {threads} -b -T {input.ref} -o {output.bam} {input.cram} 2> {log}"

# Sort BAM file
rule sort_bam:
    input:
        "aligned_reads/{sample}.unsorted.bam"
    output:
        "aligned_reads/{sample}.bam"
    log:
        "logs/sort_bam/{sample}.log"
    conda:
        "envs/chip_pipe_env.yml"
    threads: config["threads"]
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"

# Index BAM file
rule index_bam:
    input:
        "aligned_reads/{sample}.bam"
    output:
        "aligned_reads/{sample}.bam.bai"
    log:
        "logs/index_bam/{sample}.log"
    conda:
        "envs/chip_pipe_env.yml"
    threads: config["threads"]
    shell:
        "samtools index -@ {threads} {input} 2> {log}"

# Convert BAM to CRAM for storage
rule bam_to_cram:
    input:
        bam = "aligned_reads/{sample}.bam",
        ref = config["reference_genome"]
    output:
        cram = os.path.join(config["storing_cram"], "{sample}.cram")
    log:
        "logs/bam_to_cram/{sample}.log"
    conda:
        "envs/chip_pipe_env.yml"
    threads: config["threads"]
    shell:
        "samtools view -@ {threads} -C -T {input.ref} -o {output.cram} {input.bam} 2> {log}"

# Index CRAM file
rule index_cram:
    input:
        os.path.join(config["storing_cram"], "{sample}.cram")
    output:
        os.path.join(config["storing_cram"], "{sample}.cram.crai")
    log:
        "logs/index_cram/{sample}.log"
    conda:
        "envs/chip_pipe_env.yml"
    threads: config["threads"]
    shell:
        "samtools index -@ {threads} {input} 2> {log}"