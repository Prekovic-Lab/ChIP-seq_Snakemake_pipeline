# ChIP_seq_pipeline
WIP: General idea fastq -> visualization + analysis for ChIP data

```mermaid
graph TD;
    fastq --> bam
    bam --> peak calling (MACS3)
    peak calling (MACS3) --> consensus peaks (mscp)
    consensus peaks (mscp) --> peak sites
```
