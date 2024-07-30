# Snakefile

import pandas as pd
from pathlib import Path

# Load sample information
samples = pd.read_csv("samples.txt", sep="\t")

# Load configuration
configfile: "config/peaks_n_bigwigs.yaml"

# Get unique conditions
conditions = samples['condition'].unique()

# Create a list of condition pairs for intervene
intervene_pairs = [f"{pair[0]}_vs_{pair[1]}" for pair in config["intervene_compare"]]

# Define the final output files
rule all:
    input:
        macs3_peaks = expand("results/macs3/{sample}_peaks.{peak_type}",
                             sample=samples['filename'],
                             peak_type=config["chosen_peaks"]),
        merged_peaks = expand("results/merged_replicates/{condition}.bed",
                              condition=conditions) if config["merge_replicates"] else [],
        intervene_results = expand("results/intervene/{conditions}/Intervene_venn.pdf",
                                   conditions=intervene_pairs) if config["merge_replicates"] else [],
        merged_bams = expand("results/merged_bams/{condition}.bam", condition=conditions) if config["merge_replicates"] else [],
        merged_bam_indices = expand("results/merged_bams/{condition}.bam.bai", condition=conditions) if config["merge_replicates"] else [],
        bigwigs = expand("results/bigwig/{sample}.bw", 
                         sample=conditions if config["merge_replicates"] else samples['filename'])

rule macs3_callpeak:
    input:
        bam = "aligned_reads/{sample}.bam",
        control = lambda wildcards: f"aligned_reads/{samples.loc[samples['filename'] == wildcards.sample, 'control'].iloc[0]}.bam" if config["control_macs"] else []
    output:
        peaks = "results/macs3/{sample}_peaks.{peak_type}"
    params:
        name = "{sample}",
        genome_size = config["genome_size"],
        control_flag = lambda wildcards, input: f"-c {input.control}" if config["control_macs"] else "",
        format_flag = "-f BAMPE" if config["sequencing_type"] == "paired_end" else "-f BAM",
        peak_type_params = lambda wildcards: "--broad --broad-cutoff 0.1" if wildcards.peak_type == "broadPeak" else "-q 0.01"
    log:
        "logs/macs3/{sample}_{peak_type}.log"
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        macs3 callpeak -t {input.bam} \
            {params.control_flag} \
            -n {params.name} \
            --outdir $(dirname {output.peaks}) \
            -g {params.genome_size} \
            {params.format_flag} \
            {params.peak_type_params} \
            2> {log}
        """

rule mspc_merge_replicates:
    input:
        peaks = lambda wildcards: expand("results/macs3/{sample}_peaks.{peak_type}",
                                         sample=samples.loc[samples['condition'] == wildcards.condition, 'filename'],
                                         peak_type=config["chosen_peaks"])
    output:
        consensus_peaks = "results/merged_replicates/{condition}.bed"
    params:
        config_json = "config/config.json",
        output_dir = lambda wildcards, output: Path(output.consensus_peaks).parent / wildcards.condition
    log:
        "logs/mspc/{condition}.log"
    shell:
        """
        mspc -i {input.peaks} \
            -r bio -w 1e-4 -s 1e-8 -a 0.05 \
            -p {params.config_json} \
            --excludeHeader \
            -o {params.output_dir} 2> {log}
        cp {params.output_dir}/ConsensusPeaks.bed {output.consensus_peaks}
        """

rule intervene_compare:
    input:
        peaks = lambda wildcards: expand("results/merged_replicates/{condition}.bed",
                                         condition=wildcards.conditions.split("_vs_"))
    output:
        venn = "results/intervene/{conditions}/Intervene_venn.pdf"
    params:
        labels = lambda wildcards: ",".join(wildcards.conditions.split("_vs_"))
    log:
        "logs/intervene/{conditions}.log"
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        intervene venn -i {input.peaks} \
            --output $(dirname {output.venn}) \
            --save-overlaps \
            --names {params.labels} \
            2> {log}
        """

rule merge_bams:
    input:
        lambda wildcards: expand("aligned_reads/{sample}.bam",
                                 sample=samples[samples['condition'] == wildcards.condition]['filename'])
    output:
        bam = "results/merged_bams/{condition}.bam"
    log:
        "logs/merge_bams/{condition}.log"
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        samtools merge {output.bam} {input} 2> {log}
        """

rule index_merged_bams:
    input:
        "results/merged_bams/{condition}.bam"
    output:
        "results/merged_bams/{condition}.bam.bai"
    log:
        "logs/index_merged_bams/{condition}.log"
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        samtools index {input} {output} 2> {log}
        """

rule generate_bigwig:
    input:
        bam = lambda wildcards: f"results/merged_bams/{wildcards.sample}.bam" if config["merge_replicates"] else f"aligned_reads/{wildcards.sample}.bam",
        bai = lambda wildcards: f"results/merged_bams/{wildcards.sample}.bam.bai" if config["merge_replicates"] else f"aligned_reads/{wildcards.sample}.bam.bai"
    output:
        bigwig = "results/bigwig/{sample}.bw"
    params:
        bin_size = config.get("bamCoverage_binSize", 1),
        effective_genome_size = config.get("bamCoverage_effectiveGenomeSize", 2913022398),
        normalize_using = config.get("bamCoverage_normalizeUsing", "RPGC")
    threads: config.get("bamCoverage_threads", 5)
    log:
        "logs/generate_bigwig/{sample}.log"
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} \
        --binSize {params.bin_size} \
        --effectiveGenomeSize {params.effective_genome_size} \
        --normalizeUsing {params.normalize_using} \
        -p {threads} \
        2> {log}
        """