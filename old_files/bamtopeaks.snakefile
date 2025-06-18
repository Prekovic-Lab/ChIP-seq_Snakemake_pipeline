# Load configuration
configfile: 'fromsamples.yaml'

# Import modules
import os
import pandas as pd

# Read samples information from samples.txt
## MAKE THIS BEING ADDED FROM configfile (samples.txt)
samples_df = pd.read_csv("samples.txt", sep="\t")

# Check if the output folder exists, if not, create it
folder_paths = ['peaks/merged_replicates', 'aligned_reads/merged_replicates', 'bigwigs', 'data/tornado']
for folder_path in folder_paths:
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# Get unique conditions for the samples
unique_conditions = samples_df['condition'].unique()

# Rule all to specify the default target(s) to be generated
rule all:
    input:
        # Generate peaks for all samples and chosen peak types
        expand('peaks/{sample}_peaks.{peak_type}', sample=samples_df['filename'], peak_type=config['chosen_peaks']),
        expand('peaks/merged_replicates/{condition}.bed', condition=unique_conditions),
        expand('aligned_reads/merged_replicates/{condition}.bam', condition=unique_conditions),
        expand('bigwigs/{condition}.bw', condition=unique_conditions),
        #expand('peaks/intervene_output/{pair}', pair=[f'{pair[0]}_{pair[1]}' for pair in config['intervene_compare']])
        expand('peaks/intervene_output/{pair}', pair=range(len(config['intervene_compare'])))
        #'peaks/intervene_output'                ## I WILL HAVE TO FIGURE A BETTER WAY TO SOLVE THIS 

# the next two rules do peakcalling for either narrow or broad peaks
## I PROBABLY HAVE TO FIND A BETTER WAY TO INPUT.CONTROL WHEN IT DOESN'T EXIST
## I NEED TO SET THE PEAK TYPE AND REDIRECT IT TO "narrowPeak" OR "broadPeak"
if config['chosen_peaks'] == "narrowPeak":
    rule MACS3_narrow:
        input:
            bam="aligned_reads/{sample}.bam",
            control= lambda wildcards: f"aligned_reads/{samples_df.loc[samples_df['filename'] == wildcards.sample, 'control'].values[0]}.bam"
        output:
            narrowPeak='peaks/{sample}_peaks.narrowPeak',
            model_r='peaks/{sample}_model.r',
            control_lambda='peaks/{sample}_control_lambda.bdg',
            treat_pileup='peaks/{sample}_treat_pileup.bdg',
            xls='peaks/{sample}_peaks.xls',
            summits='peaks/{sample}_summits.bed'
        params:
            output_name=lambda wildcards: f"{wildcards.sample}",
            macs3_dir='peaks'
        run:
            if config['control_macs']:
                shell("""
                macs3 callpeak -t {input.bam} -c {input.control} -n {params.output_name} --outdir {params.macs3_dir} -f BAM -g hs -B -q 0.01
                """)
            else:
                shell("""
                macs3 callpeak -t {input.bam} -n {params.output_name} --outdir {params.macs3_dir} -f BAM -g hs -B -q 0.01
                """)

elif config['chosen_peaks'] == "broadPeak":
    rule MACS3_broad:
        input:
            bam='aligned_reads/{sample}.bam',
            control=lambda wildcards: f"aligned_reads/{samples_df.loc[samples_df['filename'] == wildcards.sample, 'control'].values[0]}.bam"
        output:
            broadPeak='peaks/{sample}_peaks.broadPeak',
            model_r='peaks/{sample}_model.r',
            gappedPeak='peaks/{sample}_peaks.gappedPeak',
            xls='peaks/{sample}_peaks.xls',
        params:
            output_name= lambda wildcards: f"{wildcards.sample}",
            macs3_dir="peaks"
        run:
            if config['control_macs']:
                shell("""
                macs3 callpeak -t {input.bam} -c {input.control} -n {params.output_name} --outdir {params.macs3_dir} -f BAM -g hs --broad --broad-cutoff 0.1
                """)
            else:
                shell("""
                macs3 callpeak -t {input.bam} -n {params.output_name} --outdir {params.macs3_dir} -f BAM -g hs --broad --broad-cutoff 0.1
                """)

#Rule to generate consensus peaks using MSPC. If there is one replicate then the process is skipped. 
rule MSPC:
    input:
        # List of input narrowPeak files for the current condition
        input_files=lambda wildcards: [f"peaks/{filename}_peaks.narrowPeak" for filename in samples_df[samples_df['condition'] == wildcards.condition]['filename']]
    output:
        # Output file for consensus peaks
        consensus_peaks="peaks/merged_replicates/{condition}.bed"
    params:
        config_file="config.json",
        output_dir="peaks/{condition}/"
    run:
        # Check if the length of input_files is greater than 1
        if len(input.input_files) > 1:
            print("Yes:", input.input_files)
            shell(
                """
                mspc -i {input.input_files} -r bio -w 1e-4 -s 1e-8 -a 0.05 -p {params.config_file} --excludeHeader -o {params.output_dir}
                cp peaks/{wildcards.condition}/ConsensusPeaks.bed peaks/merged_replicates/{wildcards.condition}.bed
                """
            )
        else:
            print("No replicates found for:", input.input_files)
            shell("""
                cp {input.input_files} peaks/merged_replicates/{wildcards.condition}.bed
                """)

# intervene venn compares different bed files and output three bed files (inside sets) and a venn figure. It's most useful for 2-4 way comparisons.
rule intervene:
    input:
        #input_files=lambda wildcards: expand('peaks/merged_replicates/{condition}.bed', condition=config['intervene_compare'][wildcards.pair])
        input_files=lambda wildcards: expand('peaks/merged_replicates/{condition}.bed', condition=config['intervene_compare'][int(wildcards.pair)])
        #input_files=lambda wildcards: [f'peaks/merged_replicates/{condition}.bed' for condition in wildcards.pair.split('_')]
    output:
        venn_output=directory('peaks/intervene_output/{pair}')
    shell:
        """
        intervene venn -i {input.input_files} --save-overlaps -o {output.venn_output}
        """

# this rule merges the replicates that have the same condition
## POTENTIAL PROBLEM: SAME CONDITION BUT DIFFERENT REPLICATES/SAMPLES
if config['merge_replicates']:
    rule mergeBams:
        input:
            input_files=lambda wildcards: [f"aligned_reads/{filename}.bam" for filename in samples_df[samples_df['condition'] == wildcards.condition]['filename']]
        output:
            merged_bam="aligned_reads/merged_replicates/{condition}.bam"
        run: 
            if len(input.input_files) > 1:
                print("Yes:", input.input_files)
                shell(
                """
                samtools merge aligned_reads/merged_replicates/{wildcards.condition}.bam {input.input_files}
                sleep 15
                samtools index aligned_reads/merged_replicates/{wildcards.condition}.bam -o aligned_reads/merged_replicates/{wildcards.condition}.bam.bai
                """)
            else:
                print("No replicates found for:", input.input_files)
                shell(
                    """
                    cp {input.input_files} aligned_reads/merged_replicates/{wildcards.condition}.bam
                    sleep 15
                    samtools index aligned_reads/merged_replicates/{wildcards.condition}.bam -o aligned_reads/merged_replicates/{wildcards.condition}.bam.bai
                    """)

# If we don't want to merge replicates we will have to transfer them to a folder and then do the bigwig from there. Build this later
elif config['merge_replicates'] == False:
    print("add this laterz lol")

# Create the bigwig files from the (merged) bam files
rule BamtoBigwig:
    input:
        bam_file = "aligned_reads/merged_replicates/{condition}.bam"
    output:
        bigwig = "bigwigs/{condition}.bw"
    params:
        gen_size = config['gen_size']
    shell:
        """
        bamCoverage --bam {input.bam_file} -o {output.bigwig} \
        --binSize 10 \
        #--normalizeUsing RPGC \
        #--effectiveGenomeSize {params.gen_size} \
        --extendReads
        """