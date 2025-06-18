[![Work in Progress](https://img.shields.io/badge/status-work_in_progress-orange.svg)](https://your-project-url)
[![Python Version](https://img.shields.io/badge/python-3.9.18-blue.svg)](https://www.python.org/downloads/release/python-31013/)
[![Snakemake Version](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg)](https://snakemake.github.io/)

# ChIP sequencing pipeline
WIP: General idea fastq -> visualization + analysis for ChIP data

For now: Uploading individual scripts and working on creating a snakemake pipeline. Individual scripts start from bam files. 
### peakcalling.sh
Inside this script, there are commands to run narrow or broad peaks using macs3. Install macs3 with pip install macs3. Modify the script accordingly and run through the command line.
### run_mspc.sh
The command to get the consensus peaks from .narrowPeak or .broadPeaks replicates. Adjust the command accordingly. Download [mspc](https://genometric.github.io/MSPC/docs/installation/) for your operating system and follow the instructions. Then you have to add it to the path to run it from anywhere. For Linux that means go to .bashrc (type cd, the file should be there) and add somewhere in the file `export PATH=$PATH:/path/to/mspc`. In order to run it requires the config.json file. If it doesn't work, modify (the columns of the .narrow/.broadPeak files).

    How to do on Mac:
      -    nano ~/.bash_profile
      -    export PATH=$PATH:/path/to/mspc (put your path here)
      -    source ~/.bash_profile (activate with this command)

### run_intervene.sh
Then we will use [intervene](https://github.com/asntech/intervene) to get the peaks that are lost/gained. I use as input the consensus peak bed file from `run_mspc.sh` comparing it to a different condition (e.g. 0 days vs 7 days / lost gained peaks). It needs Python 3.7 (at least to get the Venn graphs, to get the different sets it works regardless I think; you don't really need the graphs).
### Proximal genes
Now we have the .bed files that correspond to the peaks that were lost/gained/shared. From this point, we can go to [go.cistrome](http://go.cistrome.org/) and find the genes that are close to the peaks. It greatly enhances the output to have differentially expressed gene sites.

### bam to bigwig
From the bam files we have initially, we can merge the biological replicates (merge_bams.sh) and transform the merged .bam files to .bigwig files (using deeptools; bam_to_bigwig.sh). Bigwig files will be used to create the tornado plots. 
### creating the tornado plot
In order to create the tornado plot we have to create the matrix first (matrix_tornado.sh) and then plot it (plot_matrix.sh). Everything is included in the [deeptools website](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html) with great detail. We need the bigwig file(s) that was created earlier plus the bed file, which will be the regions that will be plotted (y-axis).

### doing the analysis in R
In R I am doing the analysis using [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) and [ChIPSeeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html). It enables you to skip the whole process that was previously described, however, it is less customizable. Check it and if you have any question let me know.

<img src="/Figures/chip_seq_workflow.drawio.png" width="400">



# How to run the pipeline:

First, you have to make sure a samples file exists (e.g. samples.tsv). This file contains multiple columns and is structured more or less as shown below:
|  name  | condition | replicate |  control  | seq_type | peak_type |      file_1      |      file_2      |
|  ----  | --------- | --------- | --------- | -------- | --------- | ---------------- | ---------------- |
| name1  |   cond1   |     1     |  control  |  single  |   narrow  |  name1.fastq.gz  |        -         |
| name2  |   cond2   |     2     |     -     |  paired  |   broad   | name2_1.fastq.gz | name2_2.fastq.gz |

In the config.yaml, we can edit the directories of multiple files (samples, genome, fastq, output, etc) and provide some extra options.
After modifying those, we have two options to run the pipeline

## Option 1: Running it locally on a computer:
`snakemake --use-conda --conda-frontend mamba -c 20 -j 5`

Adjust arguments (cores, jobs) accordingly.

## Option 2: Running it on HPC:
`snakemake --profile my-slurm --use-conda --conda-frontend mamba`

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
