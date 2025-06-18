import os

configfile: "config/align.yaml"

# def get_input_fastq(wildcards):
#     if config["seq_type"] == "paired":
#         return [os.path.join(config["fastq_dir"], f"{wildcards.sample}_R1.fastq.gz"),
#                 os.path.join(config["fastq_dir"], f"{wildcards.sample}_R2.fastq.gz")]
#     else:
#         return os.path.join(config["fastq_dir"], f"{wildcards.sample}.fastq.gz")

##################################################################################################
# Helper functions
##################################################################################################
def get_trimmed_fastq(wildcards):
    if config["seq_type"] == "paired":
        return [f"{config['output_dir']}/trimmed/{wildcards.sample}_R1_val_1.fq.gz",
                f"{config['output_dir']}/trimmed/{wildcards.sample}_R2_val_2.fq.gz"]
    else:
        return f"{config['output_dir']}/trimmed/{wildcards.sample}_trimmed.fq.gz"

def get_trimmed_fastqc_output():
    if config["seq_type"] == "paired":
        return expand("{output_dir}/trimmed/{sample}_R{read}_val_{read}_fastqc.zip", 
                      output_dir=config['output_dir'], sample=config["samples"], read=["1", "2"])
    else:
        return expand("{output_dir}/trimmed/{sample}_trimmed_fastqc.zip", 
                      output_dir=config['output_dir'], sample=config["samples"])

def get_fastq_input(wildcards):
    if config["seq_type"] == "paired":
        return [os.path.join(config["fastq_dir"], f"{wildcards.sample}_R1.fastq.gz"),
                os.path.join(config["fastq_dir"], f"{wildcards.sample}_R2.fastq.gz")]
    else:
        return os.path.join(config["fastq_dir"], f"{wildcards.sample}.fastq.gz")

def get_fastqc_output():
    if config["seq_type"] == "paired":
        return expand("{output_dir}/fastqc/{sample}_R{read}_fastqc.{ext}",
                      output_dir=config['output_dir'], sample=config["samples"],
                      read=[1,2], ext=["html", "zip"])
    else:
        return expand("{output_dir}/fastqc/{sample}_fastqc.{ext}",
                      output_dir=config['output_dir'], sample=config["samples"],
                      ext=["html", "zip"])

# def get_aligner_output(wildcards):
#     if config["aligner"] == "bwa":
#         return multiext(config["ref_genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
#     elif config["aligner"] == "bowtie2":
#         return multiext(config["ref_genome"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
#     else:
#         raise ValueError("Invalid aligner specified in config. Choose 'bwa' or 'bowtie2'.")

# All
rule all:
    input:
        #get_fastqc_output(),
        f"{config['output_dir']}/multiqc/original_fastqc_report.html",
        f"{config['output_dir']}/multiqc/trimmed_fastqc_report.html",
        #get_aligner_output(None),
        #expand("{output_dir}/aligned/{sample}.bam", output_dir=config['output_dir'], sample=config["samples"]),
        #expand("{output_dir}/aligned/{sample}.bam.bai", output_dir=config['output_dir'], sample=config["samples"]),
        expand("{output_dir}/deduped/{sample}.dedup.bam", output_dir=config['output_dir'], sample=config["samples"]),
        expand("{output_dir}/deduped/{sample}.dedup.bam.bai", output_dir=config['output_dir'], sample=config["samples"]),
        f"{config['output_dir']}/multiqc/alignment_report.html"

# fastqc
rule fastqc:
    input:
        get_fastq_input
    output:
        html = "{output_dir}/fastqc/{sample}_fastqc.html" if config["seq_type"] == "single" else ["{output_dir}/fastqc/{sample}_R1_fastqc.html", "{output_dir}/fastqc/{sample}_R2_fastqc.html"],
        zip = "{output_dir}/fastqc/{sample}_fastqc.zip" if config["seq_type"] == "single" else ["{output_dir}/fastqc/{sample}_R1_fastqc.zip", "{output_dir}/fastqc/{sample}_R2_fastqc.zip"]
    log:
        "{output_dir}/logs/fastqc/{sample}.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.html[0] if isinstance(output.html, list) else output.html)
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} {input} > {log} 2>&1 || echo "FastQC failed. Check {log} for details."
        """
##################################################################################################
# Trimming using trim galore
##################################################################################################
rule trim_galore_paired:
    input:
        r1 = os.path.join(config["fastq_dir"], "{sample}_R1.fastq.gz"),
        r2 = os.path.join(config["fastq_dir"], "{sample}_R2.fastq.gz")
    output:
        r1 = temp("{output_dir}/trimmed/{sample}_R1_val_1.fq.gz"),
        r2 = temp("{output_dir}/trimmed/{sample}_R2_val_2.fq.gz"),
        r1_fastqc = "{output_dir}/trimmed/{sample}_R1_val_1_fastqc.zip",
        r2_fastqc = "{output_dir}/trimmed/{sample}_R2_val_2_fastqc.zip"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.r1)
    threads: 4
    log:
        "{output_dir}/logs/trim_galore/{sample}.log"
    shell:
        """
        mkdir -p {params.output_dir}
        trim_galore --paired {input.r1} {input.r2} \
        -o {params.output_dir} --cores {threads} \
        --fastqc_args "-o {params.output_dir} -t {threads}" > {log} 2>&1 || echo "Trim Galore failed. Check {log} for details."
        """

rule trim_galore_single:
    input:
        os.path.join(config["fastq_dir"], "{sample}.fastq.gz")
    output:
        trimmed = temp("{output_dir}/trimmed/{sample}_trimmed.fq.gz"),
        fastqc = "{output_dir}/trimmed/{sample}_trimmed_fastqc.zip"
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.trimmed)
    threads: 4
    log:
        "{output_dir}/logs/trim_galore/{sample}.log"
    shell:
        "trim_galore {input} "
        "-o {params.output_dir} --cores {threads} "
        "--fastqc_args '-o {params.output_dir} -t {threads}' > {log} 2>&1"
##################################################################################################
# Multiqc_reports
##################################################################################################
rule multiqc_original:
    input:
        fastqc_original = get_fastqc_output()
    output:
        report="{output_dir}/multiqc/original_fastqc_report.html"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.report)
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {input.fastqc_original} -o {params.outdir} -n original_fastqc_report.html -f
        """

rule multiqc_trimmed:
    input:
        fastqc_trimmed = get_trimmed_fastqc_output()
    output:
        report="{output_dir}/multiqc/trimmed_fastqc_report.html"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.report)
    shell:
        "multiqc {input.fastqc_trimmed} -o {params.outdir} -n trimmed_fastqc_report.html -f"

rule multiqc_alignment:
    input:
        #bams = expand("{output_dir}/aligned/{sample}.bam", output_dir=config['output_dir'], sample=config["samples"]),
        dedup_bams = expand("{output_dir}/deduped/{sample}.dedup.bam", output_dir=config['output_dir'], sample=config["samples"]),
        dedup_metrics = expand("{output_dir}/deduped/{sample}.dedup.metrics.txt", output_dir=config['output_dir'], sample=config["samples"])
    output:
        report="{output_dir}/multiqc/alignment_report.html"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.report)
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {input.dedup_bams} {input.dedup_metrics} -o {params.outdir} -n alignment_report.html -f
        """
##################################################################################################
# Index aligners
##################################################################################################
rule bwa_index:
    input:
        ref = config["ref_genome"]
    output:
        multiext(config["ref_genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    shell:
        "bwa index {input.ref} > {log} 2>&1"

rule bowtie2_index:
    input:
        ref = config["ref_genome"]
    output:
        multiext(config["ref_genome"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    log:
        "logs/bowtie2_index.log"
    shell:
        "bowtie2-build {input.ref} {input.ref} > {log} 2>&1"
##################################################################################################
# Aligners
##################################################################################################
rule bwa_align:
    input:
        fastq = get_trimmed_fastq,
        index = multiext(config["ref_genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam = temp("{output_dir}/aligned/{sample}.bam"),
        bai = temp("{output_dir}/aligned/{sample}.bam.bai")
    threads: 8
    log:
        "{output_dir}/logs/bwa_mem/{sample}.log"
    shell:
        """
        if [ "{config[seq_type]}" = "paired" ]; then
            bwa mem -M -t {threads} {config[ref_genome]} {input.fastq} 2> {log} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        else
            bwa mem -M -t {threads} {config[ref_genome]} {input.fastq} 2> {log} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        fi
        samtools index {output.bam}
        """

rule bowtie2_align:
    input:
        fastq = get_trimmed_fastq,
        index = multiext(config["ref_genome"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    output:
        bam = temp("{output_dir}/aligned/{sample}.bam"),
        bai = temp("{output_dir}/aligned/{sample}.bam.bai")
    threads: 8
    log:
        "{output_dir}/logs/bowtie2/{sample}.log"
    shell:
        """
        if [ "{config[seq_type]}" = "paired" ]; then
            bowtie2 -p {threads} -x {config[ref_genome]} -1 {input.fastq[0]} -2 {input.fastq[1]} 2> {log} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        else
            bowtie2 -p {threads} -x {config[ref_genome]} -U {input.fastq} 2> {log} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        fi
        samtools index {output.bam}
        """
##################################################################################################
# Deduplicate
##################################################################################################
rule picard_dedup:
    input:
        bam = "{output_dir}/aligned/{sample}.bam"
    output:
        dd_bam = "{output_dir}/deduped/{sample}.dedup.bam",
        dd_bai = "{output_dir}/deduped/{sample}.dedup.bam.bai",
        metric = "{output_dir}/deduped/{sample}.dedup.metrics.txt"
    log:
        "{output_dir}/logs/picard_dedup/{sample}.log"
    shell:
        """
        picard MarkDuplicates -I {input.bam} \
        -O {output.dd_bam} \
        -M {output.metric} \
        -REMOVE_DUPLICATES false > {log} 2>&1

        samtools index {output.dd_bam}
        """

##################################################################################################
# Conditional ruleorder
##################################################################################################
if config["aligner"] == "bwa":
    ruleorder: bwa_align > bowtie2_align
elif config["aligner"] == "bowtie2":
    ruleorder: bowtie2_align > bwa_align
else:
    raise ValueError("Invalid aligner specified in config. Choose 'bwa' or 'bowtie2'.")

ruleorder: trim_galore_paired > trim_galore_single