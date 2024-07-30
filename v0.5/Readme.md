[![Work in Progress](https://img.shields.io/badge/status-work_in_progress-orange.svg)](https://your-project-url)
[![Python Version](https://img.shields.io/badge/python-3.9.18-blue.svg)](https://www.python.org/downloads/release/python-31013/)
[![Snakemake Version](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg)](https://snakemake.github.io/)

## Pipeline is divided into 4 pipelines:
3: `peak_n_bigwigs.smk`: By running this pipeline you get peaks and bigwig files starting from bam files.
- Information:
  - bam files should be indexed and sorted in `/aligned_reads`. This is the output of both pipeline 1 & 2.
  - Modify the `config/peak_n_bigwig.config` accordingly.
  - If `merge_replicates = True`, the pipeline will continue per condition. If it is `False`, you get the peak files and bigwigs only.
- Folders of the pipeline:
  - `/aligned_reads`: where the bams are located,
  - `/config`: It contains the `config/peak_n_bigwig.config`, which is needed to modify the pipeline accordingly and `config/config.json` which is needed for getting the consensus peaks (MSPC) **NOTE: this works in Linux/Windows, for mac we need a different `json` file.**
  - `/envs`: Contains the conda environment for reproducible results.
  - `/logs`: The log outputs of the pipeline. Some will be empty, but kept them, because if there is an error, we can closely inspect it.
  - `/results`: This folder will be created and will contain the multiple results from the pipeline.
  - `samples.txt`: a tab separated text file that contains the *bam file* in `/aligned_reads` (without the .bam extension), the *condition*, *replicates* and *inputs/controls* for peak calling if there are any.
- How to run:
  - Download the folders mentioned above. You get use git clone for that or manually install them.
  - Install the necessary packages, through the `.yml` file in `/envs`. Moreover, you need to install (mspc)[https://genometric.github.io/MSPC/docs/installation] and add it to the path (add to `.bashrc` the following `export PATH=$PATH:/path/to/mspc`, where you modify the `path/to/mspc`.
  - Modify the `samples.txt` and `config/peak_n_bigwig.config` according to your data.
  - Run the pipeline `snakemake -s peak_n_bigwigs.smk -c 10 -j 5`, where -c is the amount of cores and -j are the number of jobs that run in parallel (mostly for servers). If you know your setting, change those however you like.
   
If something does not work or for any questions, contact me at [Theo Chalkiadakis](mailto:t.chalkiadakis@umcutrecht.nl).
