[![Work in Progress](https://img.shields.io/badge/status-work_in_progress-orange.svg)](https://your-project-url)
[![Snakemake Version](https://img.shields.io/badge/snakemake-9.1.3-brightgreen.svg)](https://snakemake.github.io/)

# ChIP sequencing pipeline
WIP: General idea fastq -> visualization + analysis for ChIP data


# How to run the pipeline:

First, you have to make sure a samples file exists (e.g. samples.tsv). This file contains multiple columns and is structured more or less as shown below:
|  name  | condition | replicate |  control  | seq_type | peak_type |      file_1      |      file_2      |
|  ----  | --------- | --------- | --------- | -------- | --------- | ---------------- | ---------------- |
| name1  |   cond1   |     1     |  control  |  single  |   narrow  |  name1.fastq.gz  |        -         |
| name2  |   cond2   |     2     |     -     |  paired  |   broad   | name2_1.fastq.gz | name2_2.fastq.gz |

In the config.yaml, we can edit the directories of multiple files (samples, genome, fastq, output, etc) and provide some extra options.
After modifying those, we have two options to run the pipeline

## Option 1: Running it locally on a computer:
`snakemake --use-conda -c 20 -j 5`

Adjust arguments (cores, jobs) accordingly.

## Option 2: Running it on HPC:
`snakemake --profile my-slurm --use-conda`

However, before that being able to work, you need to set the my-slurm profile.
First, install [snakemake-executor-plugin-slurm](https://bioconda.github.io/recipes/snakemake-executor-plugin-slurm/README.html) in your environment.
Then set up a config.yaml file and save it in `~/.config/snakemake/`.
It should look like this (more or less):
```
executor: slurm
jobs: 50
default-resources:
  - mem_mb=4000
  - runtime=180
  - disk_mb=10000
  - tmpdir="/path/to/temporary/directory"
#set-resources:
#  - myrule:mem_mb=8000
#  - myrule:runtime=04:00:00
latency-wait: 60
rerun-incomplete: true
keep-going: true
printshellcmds: true
use-conda: true
```

For both options, keep in mind that `--use-conda` will take some time in the beginning, by creating the necessary environments. 


To-Do (Snakemake ChIP-Seq pipeline):
- [ ] UMI-tools options, in case they exist.
- [ ] Integrate the R-scripts ( diffbind / ChipSeeker / pathway analysis).
- [ ] Motif analysis (FIMO vs STREME vs HOMER)

In order to make the pipeline work you have to install [mspc](https://genometric.github.io/MSPC/docs/installation/).
Please visit the address and download the package according to the instructions. Put it in the home directory (~ / should be findable when typing cd in linux environment).
You can check if that worked if you typed `~/mspc/mspc` and get the message that it requires more input arguments. 
This is a workaround, until I make it work somehow else. 

### Early template
<img src="/Figures/chip_seq_workflow.drawio.png" width="400">




