#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cram_files = './cram_files'
params.reference = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'

cram_files = Channel.fromPath(params.cram_files + '/*.cram')
reference = file(params.reference)

process ConvertCRAMtoBAM {
    input:
    path cram_file
    path reference

    output:
    path "${cram_file.baseName}.bam"

    script:
    """
    samtools view -b -T $reference -o ${cram_file.baseName}.bam $cram_file
    """
}

process IndexBAM {
    input:
    path bam_file

    output:
    path "${bam_file}.bai"

    script:
    """
    samtools index $bam_file ${bam_file}.bai
    """
}

workflow {
    bam_files = ConvertCRAMtoBAM(cram_files, reference)
    IndexBAM(bam_files)
}

