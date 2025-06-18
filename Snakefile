import pandas as pd
import os
from pathlib import Path
import re

# Load configuration
configfile: "config.yaml"

# Read sample information with error handling
try:
    samples = pd.read_csv(config["samples"], sep=r'\s+')
except FileNotFoundError:
    raise FileNotFoundError(f"Sample file not found: {config['samples']}")
except Exception as e:
    raise ValueError(f"Error reading sample file: {e}")

# Validate required columns
required_columns = ['name', 'seq_type', 'file_1', 'condition']
missing_columns = [col for col in required_columns if col not in samples.columns]
if missing_columns:
    raise ValueError(f"Missing required columns in samples file: {missing_columns}")

# Validate sample names are unique
if samples['name'].duplicated().any():
    raise ValueError("Sample names must be unique")

# Add peak_type column if missing
if 'peak_type' not in samples.columns:
    print("Warning: 'peak_type' column not found in samples.txt. Adding column with default value 'narrow'")
    samples['peak_type'] = 'narrow'

# Set peak_type default for missing values
samples['peak_type'] = samples['peak_type'].fillna('narrow')

# Validate seq_type values
valid_seq_types = ['paired', 'single']
invalid_seq_types = samples[~samples['seq_type'].isin(valid_seq_types)]['seq_type'].unique()
if len(invalid_seq_types) > 0:
    raise ValueError(f"Invalid seq_type values: {invalid_seq_types}. Must be one of: {valid_seq_types}")

# Validate peak_type values
valid_peak_types = ['narrow', 'broad']
invalid_peak_types = samples[~samples['peak_type'].str.lower().isin(valid_peak_types)]['peak_type'].unique()
if len(invalid_peak_types) > 0:
    raise ValueError(f"Invalid peak_type values: {invalid_peak_types}. Must be one of: {valid_peak_types}")

# Create necessary directories using pathlib
output_dirs = [
    "trimmed", "logs/trim_galore", "logs/alignment", "logs/dedup", "logs/macs3",
    "QC/after_trimming", "QC/deduped", "aligned", "deduped", "peaks"
]

for dir_name in output_dirs:
    Path(config['output_dir']).joinpath(dir_name).mkdir(parents=True, exist_ok=True)

# Function to handle multiple file extensions
def strip_sample_suffix(sample, suffixes=[".sorted", ".dedup", ".bam"]):
    for suffix in suffixes:
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
    return sample

# Validate input files exist
def validate_input_files():
    """Validate that all input FASTQ files exist"""
    missing_files = []
    for _, row in samples.iterrows():
        file1_path = Path(config['fastq_dir']) / row['file_1']
        if not file1_path.exists():
            missing_files.append(str(file1_path))
        
        if row['seq_type'] == 'paired' and 'file_2' in row and pd.notna(row['file_2']):
            file2_path = Path(config['fastq_dir']) / row['file_2']
            if not file2_path.exists():
                missing_files.append(str(file2_path))
    
    if missing_files:
        raise FileNotFoundError(f"Missing input files: {missing_files}")

# Validate input files at startup
validate_input_files()

# Helper functions
def get_input_files(wildcards):
    """Get input FASTQ files for a sample"""
    row = samples.loc[samples['name'] == strip_sample_suffix(wildcards.sample)].iloc[0]
    files = {'r1': f"{config['fastq_dir']}/{row['file_1']}"}
    
    if row['seq_type'] == 'paired':
        if 'file_2' not in row or pd.isna(row['file_2']):
            raise ValueError(f"Sample {wildcards.sample} is paired but missing file_2")
        files['r2'] = f"{config['fastq_dir']}/{row['file_2']}"
    
    return files

def get_trimmed_files(wildcards):
    """Get trimmed files for a sample"""
    row = samples.loc[samples['name'] == strip_sample_suffix(wildcards.sample)].iloc[0]
    if row['seq_type'] == 'paired':
        return [
            f"{config['output_dir']}/trimmed/{wildcards.sample}_R1_val_1.fq.gz",
            f"{config['output_dir']}/trimmed/{wildcards.sample}_R2_val_2.fq.gz"
        ]
    return [f"{config['output_dir']}/trimmed/{wildcards.sample}_trimmed.fq.gz"]

def get_multiqc_input(wildcards):
    """Get all FastQC files for MultiQC"""
    fastqc_files = []
    for _, row in samples.iterrows():
        if row['seq_type'] == 'paired':
            fastqc_files.extend([
                f"{config['output_dir']}/trimmed/{row['name']}_R1_val_1_fastqc.zip",
                f"{config['output_dir']}/trimmed/{row['name']}_R2_val_2_fastqc.zip"
            ])
        else:
            fastqc_files.append(f"{config['output_dir']}/trimmed/{row['name']}_trimmed_fastqc.zip")
    return fastqc_files

def get_aligner_index(wildcards): 
    """Get aligner index files"""
    aligner = config.get("aligner", "bwa").lower()
    
    if aligner == "bwa":
        return multiext(config["ref_genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    elif aligner == "bowtie2":
        return multiext(config["ref_genome"], ".1.bt2", ".2.bt2", ".3.bt2", 
                       ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    else:
        raise ValueError(f"Unsupported aligner: {aligner}")

def get_macs_input(wildcards):
    """Get input files for MACS3, including control BAM if specified"""
    row = samples.loc[samples['name'] == wildcards.sample].iloc[0]
    
    inputs = {
        'treatment': f"{config['output_dir']}/deduped/{wildcards.sample}.dedup.bam"
    }
    
    # Always add control key, but set to empty string if no control
    if 'control' in row and pd.notna(row['control']) and row['control'] != '-':
        inputs['control'] = f"{config['output_dir']}/deduped/{row['control']}.dedup.bam"
    else:
        inputs['control'] = ""  # Empty string when no control
    
    return inputs

def get_macs_params(wildcards):
    """Get MACS3 parameters based on sample information"""
    row = samples.loc[samples['name'] == wildcards.sample].iloc[0]
    params = []
    
    # Add format parameter
    params.append('-f BAMPE' if row['seq_type'] == 'paired' else '-f BAM')
    
    # Add control parameter if specified
    if 'control' in row and pd.notna(row['control']) and row['control'] != '-':
        params.append(f"-c {config['output_dir']}/deduped/{row['control']}.dedup.bam")
    
    # Add genome size
    params.append(f"-g {config.get('genome_size', 'hs')}")
    
    return " ".join(params)

def get_peak_files():
    """Get the list of expected peak files"""
    peak_files = []
    # Only include samples that are not input controls
    non_input_samples = samples[samples['condition'].str.lower() != 'input']
    
    for _, row in non_input_samples.iterrows():
        if row['peak_type'].lower() == 'narrow':
            peak_files.append(f"{config['output_dir']}/peaks/{row['name']}_peaks.narrowPeak")
        else:
            peak_files.append(f"{config['output_dir']}/peaks/{row['name']}_peaks.broadPeak")
    
    return peak_files

def get_sample_by_seq_type(seq_type):
    """Get samples filtered by sequencing type"""
    return samples[samples['seq_type'] == seq_type]['name'].tolist()

def get_samples_by_peak_type(peak_type):
    """Get samples filtered by peak type and exclude input samples"""
    non_input_samples = samples[samples['condition'].str.lower() != 'input']
    return non_input_samples[non_input_samples['peak_type'].str.lower() == peak_type]['name'].tolist()

# New helper functions
def get_bigwig_input():
    """Get BAM files for bigwig generation, optionally merging replicates"""
    bigwig_inputs = []
    
    if config.get('bigwig', {}).get('merge_replicates', False):
        # Group by condition (excluding Input samples)
        non_input_samples = samples[samples['condition'].str.lower() != 'input']
        conditions = non_input_samples['condition'].unique()
        
        for condition in conditions:
            bigwig_inputs.append(f"{config['output_dir']}/bigwig/{condition}_merged.bw")
    else:
        # Individual bigwig files for non-input samples
        non_input_samples = samples[samples['condition'].str.lower() != 'input']
        for _, row in non_input_samples.iterrows():
            bigwig_inputs.append(f"{config['output_dir']}/bigwig/{row['name']}.bw")
    
    return bigwig_inputs

def get_mspc_input():
    """Get peak files grouped by condition for MSPC"""
    mspc_inputs = []
    non_input_samples = samples[samples['condition'].str.lower() != 'input']
    conditions = non_input_samples['condition'].unique()
    
    for condition in conditions:
        mspc_inputs.append(f"{config['output_dir']}/mspc/{condition}_consensus.bed")
    
    return mspc_inputs

def get_condition_names():
    """Get condition names for venn diagram labels"""
    non_input_samples = samples[samples['condition'].str.lower() != 'input']
    return list(non_input_samples['condition'].unique())

def get_condition_peak_files(condition):
    """Get peak files for a specific condition"""
    condition_samples = samples[
        (samples['condition'] == condition) & 
        (samples['condition'].str.lower() != 'input')
    ]
    
    peak_files = []
    for _, row in condition_samples.iterrows():
        if row['peak_type'].lower() == 'narrow':
            peak_files.append(f"{config['output_dir']}/peaks/{row['name']}_peaks.narrowPeak")
        else:
            peak_files.append(f"{config['output_dir']}/peaks/{row['name']}_peaks.broadPeak")
    
    return peak_files

def get_condition_bam_files(condition):
    """Get BAM files for a specific condition"""
    condition_samples = samples[
        (samples['condition'] == condition) & 
        (samples['condition'].str.lower() != 'input')
    ]
    
    return [f"{config['output_dir']}/deduped/{row['name']}.dedup.bam" 
            for _, row in condition_samples.iterrows()]

def get_venn_comparisons():
    """Get all comparisons to perform"""
    if 'venn_pairwise' in config:
        comparisons = []
        for conditions_list in config['venn_pairwise']:
            if len(conditions_list) > 6:
                raise ValueError(f"InterVene supports maximum 6 conditions, got {len(conditions_list)}")
            if len(conditions_list) < 2:
                raise ValueError(f"Need at least 2 conditions for comparison, got {len(conditions_list)}")
            
            comparison_name = "_vs_".join(conditions_list)
            comparisons.append((comparison_name, conditions_list))
        return comparisons
    else:
        # Fallback if no venn_pairwise defined
        return []

def get_comparison_input(comparison_name, conditions):
    """Get MSPC input files for a specific comparison"""
    return [f"{config['output_dir']}/mspc/{condition}_consensus.bed" 
            for condition in conditions]

def get_comparison_targets():
    """Get all comparison targets safely"""
    comparisons = get_venn_comparisons()
    if not comparisons:
        return []
    
    targets = []
    comparison_names = [comp[0] for comp in comparisons]
    
    # Venn diagram targets
    targets.extend(expand(f"{config['output_dir']}/venn/{{comparison}}/Intervene_venn.pdf",
                         comparison=comparison_names))
    
    # Shared peaks targets
    targets.extend(expand(f"{config['output_dir']}/venn/{{comparison}}/shared_peaks.bed",
                         comparison=comparison_names))
    
    # Heatmap targets  
    targets.extend(expand(f"{config['output_dir']}/heatmaps/{{comparison}}/consensus_peaks_heatmap.pdf",
                         comparison=comparison_names))
    
    return targets

def get_intervene_bed_files(comparison):
    """Get all InterVene BED files from the sets directory"""
    import glob
    
    venn_sets_dir = f"{config['output_dir']}/venn/{comparison}/sets"
    
    # Get all BED files from the sets directory
    bed_files = glob.glob(f"{venn_sets_dir}/*.bed")
    
    # Sort to ensure consistent order
    bed_files.sort()
    
    return bed_files

# Target rule
rule all:
    input:
        # QC reports
        f"{config['output_dir']}/QC/after_trimming/multiqc_report.html",
        f"{config['output_dir']}/QC/deduped/multiqc_report.html",
        # Aligned BAM files
        expand(f"{config['output_dir']}/aligned/{{sample}}.bam", 
               sample=samples['name']),
        # ALL BAM files (including inputs)
        expand(f"{config['output_dir']}/deduped/{{sample}}.dedup.bam", 
               sample=samples['name']),
        # Peak files (exclude input samples)
        get_peak_files(),
        # NEW OUTPUTS:
        # BigWig files
        get_bigwig_input(),
        # MSPC consensus peaks
        get_mspc_input(),
        # Venn diagram
        expand(f"{config['output_dir']}/venn/{{comparison}}/Intervene_venn.pdf",
               comparison=[comp[0] for comp in get_venn_comparisons()]),
        # Heatmaps
        #f"{config['output_dir']}/heatmaps/consensus_peaks_heatmap.pdf"
        #expand(f"{config['output_dir']}/heatmaps/{{comparison}}/consensus_peaks_heatmap.pdf",
        #        comparison=get_condition_names())
        get_comparison_targets()
                
# Trimming rules 
rule trim_galore_paired:
    input:
        unpack(get_input_files)
    output:
        r1 = f"{config['output_dir']}/trimmed/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{config['output_dir']}/trimmed/{{sample}}_R2_val_2.fq.gz",
        qc1 = f"{config['output_dir']}/trimmed/{{sample}}_R1_val_1_fastqc.zip",
        qc2 = f"{config['output_dir']}/trimmed/{{sample}}_R2_val_2_fastqc.zip",
        report1 = f"{config['output_dir']}/trimmed/{{sample}}_R1.fastq_trimming_report.txt",
        report2 = f"{config['output_dir']}/trimmed/{{sample}}_R2.fastq_trimming_report.txt"
    log:
        f"{config['output_dir']}/logs/trim_galore/{{sample}}.log"
    params:
        output_dir = f"{config['output_dir']}/trimmed",
        extra = config.get('trim_galore_extra', '--illumina --stringency 3'),
        quality = config.get('trim_quality', 20),
        min_length = config.get('min_length', 20)
    threads: 
        config.get('trim_threads', 6)
    resources:
        mem_mb = config.get('trim_mem_mb', 8000),
        runtime = config.get('trim_runtime', 120),
        tmpdir = config["tmp_dir"]
    wildcard_constraints:
        sample = "|".join(get_sample_by_seq_type('paired'))
    conda:
        "envs/trim_galore.yaml"
    shell:
        """
        trim_galore --gzip \
            --paired \
            --fastqc \
            --fastqc_args '-t {threads}' \
            --quality {params.quality} \
            --length {params.min_length} \
            {params.extra} \
            --output_dir {params.output_dir} \
            {input.r1} {input.r2} \
            &> {log}
        
        # Rename outputs to match expected standardized names
        r1_base=$(basename {input.r1} .fastq.gz)
        r2_base=$(basename {input.r2} .fastq.gz)
        
        # Function to safely move files (only if different)
        safe_mv() {{
            if [ "$1" != "$2" ] && [ -f "$1" ]; then
                mv "$1" "$2"
            elif [ "$1" = "$2" ]; then
                echo "Source and destination are the same, skipping: $1"
            fi
        }}
        
        safe_mv {params.output_dir}/${{r1_base}}_val_1.fq.gz {output.r1}
        safe_mv {params.output_dir}/${{r2_base}}_val_2.fq.gz {output.r2}
        safe_mv {params.output_dir}/${{r1_base}}_val_1_fastqc.zip {output.qc1}
        safe_mv {params.output_dir}/${{r2_base}}_val_2_fastqc.zip {output.qc2}
        safe_mv {params.output_dir}/${{r1_base}}.fastq.gz_trimming_report.txt {output.report1}
        safe_mv {params.output_dir}/${{r2_base}}.fastq.gz_trimming_report.txt {output.report2}
        """

rule trim_galore_single:
    input:
        unpack(get_input_files)
    output:
        r1 = f"{config['output_dir']}/trimmed/{{sample}}_trimmed.fq.gz",
        qc = f"{config['output_dir']}/trimmed/{{sample}}_trimmed_fastqc.zip",
        report = f"{config['output_dir']}/trimmed/{{sample}}.fastq_trimming_report.txt"
    log:
        f"{config['output_dir']}/logs/trim_galore/{{sample}}.log"
    params:
        output_dir = f"{config['output_dir']}/trimmed",
        extra = config.get('trim_galore_extra', '--illumina --stringency 3'),
        quality = config.get('trim_quality', 20),
        min_length = config.get('min_length', 20)
    threads: 
        config.get('trim_threads', 6)
    resources:
        mem_mb = config.get('trim_mem_mb', 8000),
        runtime = config.get('trim_runtime', 120),
        tmpdir = config["tmp_dir"]
    wildcard_constraints:
        sample = "|".join(get_sample_by_seq_type('single'))
    conda:
        "envs/trim_galore.yaml"
    shell:
        """
        trim_galore --gzip \
            --fastqc \
            --fastqc_args '-t {threads}' \
            --quality {params.quality} \
            --length {params.min_length} \
            {params.extra} \
            --output_dir {params.output_dir} \
            {input.r1} \
            &> {log}
        
        # Rename outputs to match expected standardized names
        input_base=$(basename {input.r1} .fastq.gz)
        
        # Function to safely move files (only if different)
        safe_mv() {{
            if [ "$1" != "$2" ] && [ -f "$1" ]; then
                mv "$1" "$2"
            elif [ "$1" = "$2" ]; then
                echo "Source and destination are the same, skipping: $1"
            fi
        }}
        
        safe_mv {params.output_dir}/${{input_base}}_trimmed.fq.gz {output.r1}
        safe_mv {params.output_dir}/${{input_base}}_trimmed_fastqc.zip {output.qc}
        safe_mv {params.output_dir}/${{input_base}}.fastq.gz_trimming_report.txt {output.report}
        """

# Alignment rules
rule align_bwa:
    input:
        reads = get_trimmed_files,
        index = get_aligner_index
    output:
        bam = f"{config['output_dir']}/aligned/{{sample}}.bam",
        bai = f"{config['output_dir']}/aligned/{{sample}}.bam.bai"
    log:
        f"{config['output_dir']}/logs/alignment/{{sample}}_bwa.log"
    params:
        index = config["ref_genome"],
        extra = config.get('bwa_extra', '-M'),
        rg = lambda w: f"@RG\\tID:{w.sample}\\tSM:{w.sample}\\tPL:illumina"
    threads: 
        config.get('align_threads', 14)
    wildcard_constraints:
        sample = "|".join(samples['name'])
    resources:
        mem_mb = config.get('align_mem_mb', 32000),
        runtime = config.get('align_runtime', 480),
        tmpdir = config["tmp_dir"]
    conda:
        "envs/alignment.yaml"
    shell:
        """
        bwa mem -t {threads} {params.extra} -R '{params.rg}' \
            {params.index} {input.reads} | \
            samtools view -@ {threads} -bS - | \
            samtools sort -@ {threads} -o {output.bam} - &> {log}
        
        samtools index {output.bam} &>> {log}
        """

rule align_bowtie2:
    input:
        reads = get_trimmed_files,
        index = get_aligner_index
    output:
        bam = f"{config['output_dir']}/aligned/{{sample}}.bam",
        bai = f"{config['output_dir']}/aligned/{{sample}}.bam.bai"
    log:
        f"{config['output_dir']}/logs/alignment/{{sample}}_bowtie2.log"
    params:
        index = config["ref_genome"],
        extra = config.get('bowtie2_extra', '--very-sensitive')
    threads: 
        config.get('align_threads', 14)
    wildcard_constraints:
        sample = "|".join(samples['name'])
    resources:
        mem_mb = config.get('align_mem_mb', 32000),
        runtime = config.get('align_runtime', 480),
        tmpdir = config["tmp_dir"]
    conda:
        "envs/alignment.yaml"
    shell:
        """
        if [ $(echo {input.reads} | wc -w) -eq 2 ]; then
            bowtie2 -p {threads} {params.extra} -x {params.index} \
                -1 $(echo {input.reads} | cut -d' ' -f1) \
                -2 $(echo {input.reads} | cut -d' ' -f2) | \
                samtools view -@ {threads} -bS - | \
                samtools sort -@ {threads} -o {output.bam} - 2> {log}
        else
            bowtie2 -p {threads} {params.extra} -x {params.index} \
                -U {input.reads} | \
                samtools view -@ {threads} -bS - | \
                samtools sort -@ {threads} -o {output.bam} - 2> {log}
        fi
        
        samtools index {output.bam} 2>> {log}
        """

rule remove_duplicates:
    input:
        f"{config['output_dir']}/aligned/{{sample}}.bam"
    output:
        bam = f"{config['output_dir']}/deduped/{{sample}}.dedup.bam",
        metrics = f"{config['output_dir']}/deduped/{{sample}}.dedup.metrics.txt"
    log:
        f"{config['output_dir']}/logs/dedup/{{sample}}.log"
    resources:
        mem_mb = config.get('picard_mem_mb', 8000),
        runtime = config.get('picard_runtime', 180),
        tmpdir = config.get("tmp_dir")
    conda:
        "envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --MAX_RECORDS_IN_RAM 2000000 \
            --TMP_DIR "{resources.tmpdir}" \
            --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
            --REMOVE_DUPLICATES false \
            --ASSUME_SORTED true \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT &> {log}
        """

# MACS3 Peak Calling Rules - Separate for narrow and broad peaks
rule macs3_narrow:
    input:
        unpack(get_macs_input)
    output:
        narrowpeak = f"{config['output_dir']}/peaks/{{sample}}_peaks.narrowPeak",
        xls = f"{config['output_dir']}/peaks/{{sample}}_peaks.xls",
        summits = f"{config['output_dir']}/peaks/{{sample}}_summits.bed"
    log:
        f"{config['output_dir']}/logs/macs3/{{sample}}_narrow.log"
    params:
        name = "{sample}",
        outdir = f"{config['output_dir']}/peaks",
        macs_params = get_macs_params,
        extra = config.get('macs3_extra', '--keep-dup all')
    wildcard_constraints:
        sample = "|".join(get_samples_by_peak_type('narrow'))
    resources:
        mem_mb = config.get('macs3_mem_mb', 4000),
        runtime = config.get('macs3_runtime', 60),
        tmpdir = config.get("tmp_dir")
    conda:
        "envs/macs3.yaml"
    shell:
        """
        # Validate treatment BAM file
        samtools quickcheck {input.treatment} || {{ echo "Treatment BAM file is corrupted: {input.treatment}" >> {log}; exit 1; }}
        
        # Check control file if it exists and is not empty
        if [ -n "{input.control}" ] && [ "{input.control}" != "" ]; then
            samtools quickcheck {input.control} || {{ echo "Control BAM file is corrupted: {input.control}" >> {log}; exit 1; }}
        fi
        
        macs3 callpeak \
            -t {input.treatment} \
            {params.macs_params} \
            -q 0.01 \
            --outdir {params.outdir} \
            -n {params.name} \
            {params.extra} \
            2> {log}
        """

rule macs3_broad:
    input:
        unpack(get_macs_input)
    output:
        broadpeak = f"{config['output_dir']}/peaks/{{sample}}_peaks.broadPeak",
        xls = f"{config['output_dir']}/peaks/{{sample}}_peaks.xls",
        gappedpeak = f"{config['output_dir']}/peaks/{{sample}}_peaks.gappedPeak"
    log:
        f"{config['output_dir']}/logs/macs3/{{sample}}_broad.log"
    params:
        name = "{sample}",
        outdir = f"{config['output_dir']}/peaks",
        macs_params = get_macs_params,
        extra = config.get('macs3_extra', '--keep-dup all')
    wildcard_constraints:
        sample = "|".join(get_samples_by_peak_type('broad'))
    resources:
        mem_mb = config.get('macs3_mem_mb', 4000),
        runtime = config.get('macs3_runtime', 60),
        tmpdir = config.get("tmp_dir")
    conda:
        "envs/macs3.yaml"
    shell:
        """
        # Validate treatment BAM file
        samtools quickcheck {input.treatment} || {{ echo "Treatment BAM file is corrupted: {input.treatment}" >> {log}; exit 1; }}
        
        # Check control file if it exists and is not empty
        if [ -n "{input.control}" ] && [ "{input.control}" != "" ]; then
            samtools quickcheck {input.control} || {{ echo "Control BAM file is corrupted: {input.control}" >> {log}; exit 1; }}
        fi
        
        macs3 callpeak \
            -t {input.treatment} \
            {params.macs_params} \
            --broad --broad-cutoff 0.1 \
            --outdir {params.outdir} \
            -n {params.name} \
            {params.extra} \
            2> {log}
        """

rule multiqc_trimming:
    input:
        get_multiqc_input
    output:
        f"{config['output_dir']}/QC/after_trimming/multiqc_report.html"
    log:
        f"{config['output_dir']}/logs/multiqc_trimming.log"
    params:
        indir = f"{config['output_dir']}/trimmed",
        outdir = f"{config['output_dir']}/QC/after_trimming"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc {params.indir} -o {params.outdir} &> {log}
        """

rule multiqc_deduped:
    input:
        expand(f"{config['output_dir']}/deduped/{{sample}}.dedup.metrics.txt", 
               sample=samples['name'])
    output:
        f"{config['output_dir']}/QC/deduped/multiqc_report.html"
    log:
        f"{config['output_dir']}/logs/multiqc_deduped.log"
    params:
        indir = f"{config['output_dir']}/deduped",
        outdir = f"{config['output_dir']}/QC/deduped"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc {params.indir} -o {params.outdir} &> {log}
        """

# Rule to merge BAM files by condition (for replicate merging)
rule merge_bams_by_condition:
    input:
        lambda wildcards: get_condition_bam_files(wildcards.condition)
    output:
        bam = f"{config['output_dir']}/merged_bams/{{condition}}_merged.bam",
        bai = f"{config['output_dir']}/merged_bams/{{condition}}_merged.bam.bai"
    log:
        f"{config['output_dir']}/logs/merge_bams/{{condition}}.log"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "envs/alignment.yaml"
    shell:
        """
        if [ $(echo {input} | wc -w) -eq 1 ]; then
            # Only one BAM file, create symlink
            ln -sf $(realpath {input}) {output.bam}
            ln -sf $(realpath {input}).bai {output.bai}
        else
            # Multiple BAM files, merge them
            samtools merge -@ {threads} {output.bam} {input} 2> {log}
            samtools index {output.bam} 2>> {log}
        fi
        """

# Rule to create individual BigWig files
rule create_bigwig_individual:
    input:
        bam = f"{config['output_dir']}/deduped/{{sample}}.dedup.bam",
        bai = f"{config['output_dir']}/deduped/{{sample}}.dedup.bam.bai"
    output:
        f"{config['output_dir']}/bigwig/{{sample}}.bw"
    log:
        f"{config['output_dir']}/logs/bigwig/{{sample}}.log"
    params:
        normalize = config.get('bigwig', {}).get('normalize', 'RPKM')
    threads: 6
    resources:
        mem_mb = config.get('macs3_mem_mb', 10000),
        runtime = config.get('macs3_runtime', 120),
        tmpdir = config.get("tmp_dir")
    wildcard_constraints:
        sample = "|".join(samples[samples['condition'].str.lower() != 'input']['name'])
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output} \
            --numberOfProcessors {threads} \
            --normalizeUsing {params.normalize} \
            --binSize 10 \
            --smoothLength 30 \
            --extendReads 150 \
            &> {log}
        """

# Rule to create merged BigWig files by condition
rule create_bigwig_merged:
    input:
        bam = f"{config['output_dir']}/merged_bams/{{condition}}_merged.bam",
        bai = f"{config['output_dir']}/merged_bams/{{condition}}_merged.bam.bai"
    output:
        f"{config['output_dir']}/bigwig/{{condition}}_merged.bw"
    log:
        f"{config['output_dir']}/logs/bigwig/{{condition}}_merged.log"
    params:
        normalize = config.get('bigwig', {}).get('normalize', 'RPKM')
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 180
    wildcard_constraints:
        condition = "|".join(samples[samples['condition'].str.lower() != 'input']['condition'].unique())
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output} \
            --numberOfProcessors {threads} \
            --normalizeUsing {params.normalize} \
            --binSize 10 \
            --smoothLength 30 \
            --extendReads 150 \
        &> {log}
        """

# Rule to run MSPC for consensus peaks
rule mspc_consensus:
    input:
        lambda wildcards: get_condition_peak_files(wildcards.condition)
    output:
        consensus = f"{config['output_dir']}/mspc/{{condition}}_consensus.bed"
    log:
        f"{config['output_dir']}/logs/mspc/{{condition}}.log"
    params:
        outdir = f"{config['output_dir']}/mspc/{{condition}}",
        stringency = config.get('mspc', {}).get('stringency', 1e-8),
        weak_threshold = config.get('mspc', {}).get('weak_threshold', 1e-4)
    wildcard_constraints:
        condition = "|".join(samples[samples['condition'].str.lower() != 'input']['condition'].unique())
    resources:
        mem_mb = 4000,
        runtime = 60
    # conda:
    #     "envs/mspc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        ~/mspc/mspc -i {input} \
            -r bio \
            -w {params.weak_threshold} \
            -s {params.stringency} \
            -o {params.outdir} \
            --excludeHeader \
            &> {log}
        
        # Copy consensus peaks
        cp {params.outdir}/ConsensusPeaks.bed {output.consensus}
        """

# Rule to run InterVenn for multiple pairwise comparisons
rule intervene_venn:
    input:
        lambda wildcards: get_comparison_input(wildcards.comparison, 
                                             dict(get_venn_comparisons())[wildcards.comparison])
    output:
        results = f"{config['output_dir']}/venn/{{comparison}}/intervene_results.txt",
        plot = f"{config['output_dir']}/venn/{{comparison}}/Intervene_venn.pdf",
        shared_peaks = f"{config['output_dir']}/venn/{{comparison}}/shared_peaks.bed"
    log:
        f"{config['output_dir']}/logs/intervene_{{comparison}}.log"
    params:
        outdir = f"{config['output_dir']}/venn/{{comparison}}",
        conditions = lambda wildcards: ",".join(dict(get_venn_comparisons())[wildcards.comparison]),
        num_conditions = lambda wildcards: len(dict(get_venn_comparisons())[wildcards.comparison])
    resources:
        mem_mb = 4000,
        runtime = 30
    conda:
        "envs/intervene.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        intervene venn \
            --input {input} \
            --names {params.conditions} \
            --output {params.outdir} \
            --save-overlaps \
            &> {log}
        
        # Create the shared peaks filename pattern based on number of conditions
        # For n conditions, shared peaks file has n ones: e.g., "11" for 2 conditions, "111" for 3 conditions
        shared_pattern=$(python3 -c "print('1' * {params.num_conditions})")
        shared_file="{params.outdir}/sets/${{shared_pattern}}_*.bed"
        
        # Find and copy the shared peaks file
        if ls $shared_file 1> /dev/null 2>&1; then
            cp $shared_file {output.shared_peaks}
            echo "Found shared peaks file: $shared_file" >> {log}
        else
            touch {output.shared_peaks}
            echo "Warning: No shared peaks file found with pattern: $shared_pattern" >> {log}
            echo "Available files in sets directory:" >> {log}
            ls {params.outdir}/sets/ >> {log}
        fi
        
        echo "InterVenn analysis completed successfully" > {output.results}
        echo "Input files: {input}" >> {output.results}
        echo "Conditions compared: {params.conditions}" >> {output.results}
        echo "Shared peaks file: {output.shared_peaks}" >> {output.results}
        """

# Rule to generate all Venn comparisons
rule all_venn_comparisons:
    input:
        expand(f"{config['output_dir']}/venn/{{comparison}}/Intervene_venn.pdf",
               comparison=[comp[0] for comp in get_venn_comparisons()])

# Rule to compute matrix for each comparison
rule compute_matrix_comparison:
    input:
        venn_done = f"{config['output_dir']}/venn/{{comparison}}/intervene_results.txt",
        bigwig_files = lambda wildcards: [f"{config['output_dir']}/bigwig/{condition}_merged.bw" 
                                         for condition in dict(get_venn_comparisons())[wildcards.comparison]]
    output:
        matrix = f"{config['output_dir']}/heatmaps/{{comparison}}/consensus_peaks.matrix.gz"
    log:
        f"{config['output_dir']}/logs/compute_matrix_{{comparison}}.log"
    params:
        upstream = config.get('heatmap', {}).get('upstream', 5000),
        downstream = config.get('heatmap', {}).get('downstream', 5000),
        bin_size = config.get('heatmap', {}).get('bin_size', 50),
        reference_point = config.get('heatmap', {}).get('reference_point', 'center'),
        venn_sets_dir = f"{config['output_dir']}/venn/{{comparison}}/sets",
        num_conditions = lambda wildcards: len(dict(get_venn_comparisons())[wildcards.comparison])
    threads: 8
    conda:
        "envs/deeptools.yaml"
    resources:
        mem_mb = 16000,
        runtime = 240
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        
        # Generate expected InterVene file patterns based on number of conditions
        bed_files=""
        
        # For 2 conditions: 01, 10, 11
        # For 3 conditions: 001, 010, 011, 100, 101, 110, 111
        # etc.
        
        for i in $(seq 1 $((2**{params.num_conditions} - 1))); do
            binary=$(python3 -c "print(format($i, '0{params.num_conditions}b'))")
            pattern="{params.venn_sets_dir}/${{binary}}_*.bed"
            if ls $pattern 1> /dev/null 2>&1; then
                bed_files="$bed_files $pattern"
            fi
        done
        
        echo "Using InterVene BED files: $bed_files" >> {log}
        
        if [ -n "$bed_files" ]; then
            computeMatrix reference-point \
                --regionsFileName $bed_files \
                --scoreFileName {input.bigwig_files} \
                --outFileName {output.matrix} \
                --referencePoint {params.reference_point} \
                --beforeRegionStartLength {params.upstream} \
                --afterRegionStartLength {params.downstream} \
                --binSize {params.bin_size} \
                --numberOfProcessors {threads} \
                --skipZeros \
                &> {log}
        else
            echo "ERROR: No InterVene BED files found!" >> {log}
            exit 1
        fi
        """

# Rule to create heatmap for each comparison
rule plot_heatmap_comparison:
    input:
        matrix = f"{config['output_dir']}/heatmaps/{{comparison}}/consensus_peaks.matrix.gz"
    output:
        heatmap = f"{config['output_dir']}/heatmaps/{{comparison}}/consensus_peaks_heatmap.pdf"
    log:
        f"{config['output_dir']}/logs/plot_heatmap_{{comparison}}.log"
    params:
        colors = config.get('heatmap', {}).get('colors', 'Reds')
    resources:
        mem_mb = 4000,
        runtime = 90
    shell:
        """
        plotHeatmap \
            --matrixFile {input.matrix} \
            -out {output.heatmap} \
            --colorMap {params.colors} \
            --heatmapWidth 10 \
            &> {log}
        """

# Add conditional rule order for BigWig generation
if config.get('bigwig', {}).get('merge_replicates', False):
    ruleorder: create_bigwig_merged > create_bigwig_individual
else:
    ruleorder: create_bigwig_individual > create_bigwig_merged

# Create output directories
bigwig_dirs = ["bigwig", "merged_bams", "mspc", "venn", "heatmaps"]
for dir_name in bigwig_dirs:
    Path(config['output_dir']).joinpath(dir_name).mkdir(parents=True, exist_ok=True)

# Use appropriate alignment rule based on config
if config.get("aligner", "bwa").lower() == "bwa":
    ruleorder: align_bwa > align_bowtie2
else:
    ruleorder: align_bowtie2 > align_bwa