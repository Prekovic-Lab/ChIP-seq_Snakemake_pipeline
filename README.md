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

<img src="/chip_seq_workflow.drawio.png" width="400">

```mermaid
flowchart TD
    bam --> bai
    bam --> bigwig
    bai --> bigwig
```

To Do:
- [ ] Fix mermaid flowchart above - update it when I add more things
- [ ] Figure a way to deal with the replicates/input files (config?)
- [ ] Start the pipeline from fastq; if someone has cram/bam/whatever what will they do?
- [ ] Make it more verbose in the command line (tags)
- [ ] Add comments in the script
- [ ] Add thread option in the script
- [ ] After fixing the input and figuring out a way to deal with the replicates, continue with the second workflow (MACS3, mscp, intervene, gseapy)
- [ ] Find an alternative to cistrome go 
