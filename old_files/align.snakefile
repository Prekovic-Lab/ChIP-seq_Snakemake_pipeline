# Load configuration
configfile: 'config.yaml'

# Define chosen_aligner based on the config
chosen_aligner = config['chosen_aligner']
#chosen_peaks = config['chosen_peaks']

# Import modules
import os

# Check if the output folder exists, if not, create it
for folder in ['aligned_reads', 'fastqc']:
    if not os.path.exists(folder):
        os.makedirs(folder)

rule all:
    input:
        expand('fastqc/{sample}_fastqc.html', sample=config['samples']),
        expand('aligned_reads/{sample}.bam', sample=config['samples'])
#        expand('peaks/{sample}_peaks.{peak_type}', sample=config['samples'], peak_type=config['chosen_peaks'])

rule fastqc:
    input:
        'fastq_files/{sample}.fastq.gz'
    output:
        html='fastqc/{sample}_fastqc.html'
        #seq_length_txt='fastqc/{sample}_seq_length.txt' 
    params:
        fastqc_dir='fastqc'
    shell:
        #shell("mkdir -p {params.fastqc_dir}")
        "fastqc -o {params.fastqc_dir} {input}"

        # # Extract sequence length from HTML file
        # shell("""
        #     html_file=fastqc/{wildcards.sample}_fastqc.html
        #     seq_length=$(grep "Sequence length" "$html_file" | awk '{{print $3}}')
        #     echo "$seq_length" > fastqc/{wildcards.sample}_seq_length.txt
        # """)

# command to calculate the sequencing length
#seq_length=$(zcat A2set3_r1.fastq.gz | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}')
# html_file = 
# seq_length = extract_sequence_length_from_file(html_file)
# if seq_length is not None:
#     print("Sequence length:", seq_length)
# else:
#     print("Sequence length not found in HTML file.")

#rule trim_and_filter:
# I have to discuss what to do here. There is a tool (sickle) that trims fastq files using quality
# https://github.com/najoshi/sickle

# make sure this is correct
if chosen_aligner == 'bowtie2':
    rule align_reads:
        input:
            fastq='fastq_files/{sample}.fastq.gz',
            index='reference_genome/bowtie2_index'
        output:
            bam='aligned_reads/{sample}.bam'
        params:
            aligner='bowtie2'
        threads: 8
        shell:
            """
            bowtie2 -x {input.index} -U {input.fastq} -S /dev/stdout | samtools view -bS - | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}
            """

elif chosen_aligner == 'bwa':
    rule align_reads:
        input:
            fastq='fastq_files/{sample}.fastq.gz',
            index='reference_genome/bwa_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'
        output:
            bam='aligned_reads/{sample}.bam'
        params:
            aligner='bwa'
        threads: 8
        shell:
            """
            bwa mem -M -t {threads} {input.index} {input.fastq} | samtools view -bS - | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}
            """

# Add elif blocks for other aligners similarly
#
#
#
#

# if chosen_peaks == "narrowPeak":
#     rule MACS3:
#         input:
#             bam='aligned_reads/{sample}.bam',
#             control=lambda wildcards: f'aligned_reads/{config["control_files"][wildcards.sample]}.bam'
#         output:
#             narrowPeak='peaks/{sample}_peaks.narrowPeak',
#             model_r='peaks/{sample}_model.r',
#             control_lambda='peaks/{sample}_control_lambda.bdg',
#             treat_pileup='peaks/{sample}_treat_pileup.bdg',
#             xls='peaks/{sample}_peaks.xls',
#             summits='peaks/{sample}_summits.bed'
#         params:
#             output_name=lambda wildcards: f"{wildcards.sample}",
#             macs3_dir='peaks'
#         shell:
#             """
#             macs3 callpeak -t {input.bam} -c {input.control} -n {params.output_name} --outdir {params.macs3_dir} -f BAM -g hs -B -q 0.01
#             """

# # make sure this is correct
# elif chosen_peaks == "broadPeak":
#     rule MACS3:
#         input:
#             bam='aligned_reads/{sample}.bam',
#             control=lambda wildcards: f'aligned_reads/{config["control_files"][wildcards.sample]}.bam'
#         output:
#             broadPeak='peaks/{sample}_peaks.broadPeak',
#             model_r='peaks/{sample}_model.r',
#             gappedPeak='peaks/{sample}_peaks.gappedPeak',
#             xls='peaks/{sample}_peaks.xls',
#         params:
#             output_name= lambda wildcards: f"{wildcards.sample}",
#             macs3_dir="peaks"
#         shell:
#             """
#             macs3 callpeak -t {input_bam} -c {input.control} -n {params.output_name} --outdir {params.macs3_dir} -f BAM -g hs --broad --broad-cutoff 0.1
#             """

