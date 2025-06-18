# Snakefile

import pandas as pd
from pathlib import Path

# Load configuration
configfile: "config/peaks_n_bigwigs.yaml"

# Load sample information
samples = pd.read_csv(config["samples_file"], sep="\t")

# Get unique conditions
conditions = samples['condition'].unique()

# Create a list of condition pairs for intervene
intervene_pairs = ["_vs_".join(pair) for pair in config["intervene_compare"]]

# Define the final output files
rule all:
    input:
        macs3_peaks = expand("{output_dir}/macs3/{sample}_peaks.{peak_type}",
                             output_dir=config["output_dir"],
                             sample=samples['filename'],
                             peak_type=config["chosen_peaks"]),
        #merged_peaks = expand("{output_dir}/merged_replicates/{condition}.bed",
        #                      output_dir=config["output_dir"],
        #                      condition=conditions) if config["merge_replicates"] else [],
        intervene_results = expand("{output_dir}/intervene/{conditions}/Intervene_venn.pdf",
                                   output_dir=config["output_dir"],
                                   conditions=intervene_pairs) if config["merge_replicates"] else [],
        #merged_bams = expand("{output_dir}/merged_bams/{condition}.bam",
        #                     output_dir=config["output_dir"],
        #                     condition=conditions) if config["merge_replicates"] else [],
        #merged_bam_indices = expand("{output_dir}/merged_bams/{condition}.bam.bai",
        #                            output_dir=config["output_dir"],
        #                            condition=conditions) if config["merge_replicates"] else [],
        bigwigs = expand("{output_dir}/bigwig/{sample}.bw",
                         output_dir=config["output_dir"],
                         sample=conditions if config["merge_replicates"] else samples['filename']),
        combined_frip = f"{config['output_dir']}/frip/FRIP_scores.txt",
        pca_plot = f"{config['output_dir']}/pca/pca_plot.pdf"

rule macs3_callpeak:
    input:
        bam = config["bam_dir"] + "/{sample}.bam",
        control = lambda wildcards: f"{config['bam_dir']}/{samples.loc[samples['filename'] == wildcards.sample, 'control'].iloc[0]}.bam" if config["control_macs"] else []
    output:
        peaks = "{output_dir}/macs3/{sample}_peaks.{peak_type}"
    params:
        name = "{sample}",
        genome_size = config["genome_size"],
        control_flag = lambda wildcards, input: f"-c {input.control}" if config["control_macs"] else "",
        format_flag = "-f BAMPE" if config["seq_type"] == "paired" else "-f BAM",
        peak_type_params = lambda wildcards: "--broad --broad-cutoff 0.1" if wildcards.peak_type == "broadPeak" else "-q 0.01"
    log:
        "{output_dir}/logs/macs3/{sample}_{peak_type}.log"
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
        peaks = lambda wildcards: expand("{output_dir}/macs3/{sample}_peaks.{peak_type}",
                                         output_dir=config["output_dir"],
                                         sample=samples.loc[samples['condition'] == wildcards.condition, 'filename'],
                                         peak_type=config["chosen_peaks"])
    output:
        consensus_peaks = "{output_dir}/merged_replicates/{condition}.bed"
    params:
        config_json = "config/config.json",
        output_dir = lambda wildcards, output: Path(output.consensus_peaks).parent / wildcards.condition
    log:
        "{output_dir}/logs/mspc/{condition}.log"
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
        peaks = lambda wildcards: expand("{output_dir}/merged_replicates/{condition}.bed",
                                         output_dir=config["output_dir"],
                                         condition=wildcards.conditions.split("_vs_"))
    output:
        venn = "{output_dir}/intervene/{conditions}/Intervene_venn.pdf"
    params:
        labels = lambda wildcards: ",".join(wildcards.conditions.split("_vs_"))
    log:
        "{output_dir}/logs/intervene/{conditions}.log"
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
        lambda wildcards: expand(config["bam_dir"] + "/{sample}.bam",
                                 sample=samples[samples['condition'] == wildcards.condition]['filename'])
    output:
        bam = "{output_dir}/merged_bams/{condition}.bam"
    log:
        "{output_dir}/logs/merge_bams/{condition}.log"
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        samtools merge {output.bam} {input} 2> {log}
        """

rule index_merged_bams:
    input:
        "{output_dir}/merged_bams/{condition}.bam"
    output:
        "{output_dir}/merged_bams/{condition}.bam.bai"
    log:
        "{output_dir}/logs/index_merged_bams/{condition}.log"
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        samtools index {input} {output} 2> {log}
        """

rule generate_bigwig:
    input:
        bam = lambda wildcards: f"{config['output_dir']}/merged_bams/{wildcards.sample}.bam" if config["merge_replicates"] else f"{config['bam_dir']}/{wildcards.sample}.bam",
        bai = lambda wildcards: f"{config['output_dir']}/merged_bams/{wildcards.sample}.bam.bai" if config["merge_replicates"] else f"{config['bam_dir']}/{wildcards.sample}.bam.bai"
    output:
        bigwig = "{output_dir}/bigwig/{sample}.bw"
    params:
        bin_size = config.get("bamCoverage_binSize", 1),
        effective_genome_size = config.get("bamCoverage_effectiveGenomeSize", 2913022398),
        normalize_using = config.get("bamCoverage_normalizeUsing", "RPGC")
    threads: config.get("bamCoverage_threads", 5)
    log:
        "{output_dir}/logs/generate_bigwig/{sample}.log"
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

rule calculate_frip:
    input:
        bam = config["bam_dir"] + "/{sample}.bam",
        peaks = "{output_dir}/macs3/{sample}_peaks.{peak_type}"
    output:
        frip_score = temp("{output_dir}/frip/{sample}_{peak_type}_frip.txt")
    conda:
        "envs/chip_pipe.yml"
    shell:
        """
        total_reads=$(samtools view -c {input.bam})
        reads_in_peaks=$(bedtools intersect -a {input.bam} -b {input.peaks} -c -bed | awk '{{sum+=$NF}} END {{print sum}}')
        frip=$(awk "BEGIN {{print $reads_in_peaks / $total_reads}}")
        echo -e "{wildcards.sample}\t{wildcards.peak_type}\t$frip" > {output.frip_score}
        """

rule combine_frip_scores:
    input:
        frip_scores = expand("{output_dir}/frip/{sample}_{peak_type}_frip.txt",
                             output_dir=config["output_dir"],
                             sample=samples['filename'],
                             peak_type=config["chosen_peaks"])
    output:
        combined_frip = "{output_dir}/frip/FRIP_scores.txt"
    shell:
        """
        echo -e "Sample\tPeak_Type\tFRIP_Score" > {output.combined_frip}
        cat {input.frip_scores} >> {output.combined_frip}
        """

rule create_coverage_matrix:
    input:
        bams = expand(config["bam_dir"] + "/{sample}.bam", sample=samples['filename'])
    output:
        matrix = "{output_dir}/pca/coverage_matrix.gz"
    params:
        labels = " ".join(samples['filename'])
    conda:
        "envs/chip_pipe.yml"  # Make sure deepTools is included in this environment
    threads: config.get("multiBamSummary_threads", 4)
    shell:
        """
        multiBamSummary bins --bamfiles {input.bams} \
            --labels {params.labels} \
            --numberOfProcessors {threads} \
            -o {output.matrix}
        """

rule generate_pca_plot:
    input:
        matrix = "{output_dir}/pca/coverage_matrix.gz"
    output:
        pca_plot = "{output_dir}/pca/pca_plot.pdf"
    conda:
        "envs/chip_pipe.yml"  # Make sure deepTools is included in this environment
    shell:
        """
        plotPCA -in {input.matrix} \
            -o {output.pca_plot} \
            --plotTitle "PCA of ChIP-seq samples" \
            --outFileNameData {output.pca_plot}.txt
        """