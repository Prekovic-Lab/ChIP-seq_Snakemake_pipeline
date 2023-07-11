# ChIP_seq_pipeline
WIP: General idea fastq -> visualization + analysis for ChIP data

```mermaid
flowchart TD
    fastq --> bam
    bam --> peak_call
    peak_call --> con_peak
    con_peak --> peak_site
```
