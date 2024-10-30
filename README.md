# EVCpipeline 
## Workflow Introduction
<img src="https://github.com/Tina04021997/EVCpipeline/blob/main/workflow_logo/v0.1.png" width="95%" height="95%">

## How to run EVC pipeline
1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment
2. Make sure you have the following folders and files in your working directory:
  - main.nf
  - nextflow.config
  - conf/base.config
  - sample.csv
  - summary.py

  > You can find these files and folder in `/tscc/projects/ps-lalexandrov/shared/EVC_nextflow` but you need your own **sample.csv** file
3. Due to TSCC memory issue, you will need to modify the temp files in:
  - `params.mkdup_temp_dir` in the `main.nf` file to some folder in restricted
  -  **workDir** in the `nextflow.config` file to some folder in restricted
4. Run the following code in your working directory using an interactive node:
```
# Activate your nextflow conda environment
conda activate env_nf

# Export TSCC temp directory to any folder in restricted
export TMPDIR=/some/folder/in/restricted/

# Run nextflow
nextflow run main.nf 
```
5. Every process result and report will be stored in the **RESULT** folder

## Tool Versions

| Tool | Version |
| --- | --- |
| FastQC | v0.12.1 |
| Picard | v2.18.27 |
| samtool | v1.21 |
| bwa-mem2 | v2.2.1 |
| Conpair | v0.2 |
| Picard MarkDuplicates | v3.2.0-1 |
| gatk4 | v4.6.0.0 |
| mosdepth | v0.3.8 |
| Strelka2 | v2.9.10 |
| Mutect2 | v4.6.0.0 (gatk) |
| SAGE | v3.3 |
| MuSE2 | v2.0 |
