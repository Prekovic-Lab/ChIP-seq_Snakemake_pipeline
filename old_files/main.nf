#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cram_files = './cram_files'
params.reference = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'

cram_files = Channel.fromPath(params.cram_files + '/*.cram')
reference = file(params.reference)

process ConvertCRAMtoBAM {
    publishDir 'bam_files', mode: 'copy'
    tag "${cram_file.baseName}"
    
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
    publishDir 'bam_files', mode: 'copy'

    input:
    path bam_file

    output:
    path "${bam_file}.bai"

    script:
    """
    samtools index $bam_file ${bam_file}.bai
    """
}

process BAMtoBigwig {
    publishDir 'bigwig_files', mode: 'copy'
    
    input:
    path bam_file
    path index_file

    output:
    path "${bam_file}.bw"

    script:
    """
    bamCoverage -b $bam_file -o ${bam_file}.bw
    """
}

workflow {
    bam_files = ConvertCRAMtoBAM(cram_files, reference)
    index_files = IndexBAM(bam_files)
    BAMtoBigwig(bam_files, index_files)
}

