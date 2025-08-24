# mlm2

**mlm2** is a multiomics processing pipeline designed for meta-omic data analysis. Built using Snakemake, it orchestrates various tools to efficiently process and analyze large-scale multiomics datasets, with a focus on metagenomics and metatranscriptomics.

## Features

- Automated and reproducible meta-omic data processing
- Integration with Snakemake for scalable and modular workflows
- Built-in support for MetaPhlAn (taxonomic profiling) and HUMAnN (functional profiling)
- SLURM cluster support via Snakemake executor
- Customizable configuration via sample sheets and config files
- Support for conda/mamba environments for robust dependency management

## Requirements

- `conda` >= 24.7.1
- `snakemake` >= 8.25.3
- `mamba` for environment management
- (Optional) SLURM cluster for distributed computing

## Pipeline Processing Steps Overview

This repository implements bioinformatics workflows for three major data types, each with dedicated processing steps defined in the `rules/` directory:
### 0. Environment Setup
- **Automatic software installation**: Automatically install relevant genome assembly, microbiome profiling, and general bioinformatics tools using mamba/conda infrastructure
- **Automatic database downloads**: Pull relevant databases required for BioBakery, Kraken, and other software suits. 


### 1. AMP (Amplicon) Processing
- **QIIME2**: Implement generic setup and processing for QIIME2
  - **Raw Data QC**: Quality control and filtering of raw amplicon sequencing reads.
  - **Primer/Adapter Removal**: Trimming of primers and adapters.
  - **Denoising**: Sequence denoising (e.g., DADA2) to generate Amplicon Sequence Variants (ASVs).
  - **Chimera Removal**: Detection and removal of chimeric sequences.
  - **Taxonomic Assignment**: Classification of ASVs using reference databases.
  - **Diversity Analyses**: Calculation of alpha and beta diversity metrics.
  - **Visualization**: Generation of summary plots and reports.

### 2. MGX (Metagenomics) Processing
- **Read QC and Filtering**: Quality control and removal of low-quality reads using **KneadData**.
  - **Host/Contaminant Removal**: Filtering out human contaminating sequences using T2T genome assembly.
  - **rRNA Gene Removal**: Filtering of ribosomal RNA gene reads using SILVA.
  - **Overrepresentative Sequence Removal**: Filtering out overrepresentative sequences
- **Abundance Estimation**: Quantification of taxa and functions.
  - **Taxonomic Profiling**: Assignment of taxonomy to reads (MetaPhlAn4, Kraken2, Braken, KrakenUniq, mOTUs).
  - **Functional Profiling**: Functional annotation using HUMAnN
- **Assembly**: Metagenome assembly from filtered reads.
  - **Gene Prediction & Annotation**: Identification and annotation of genes and functional elements.


### 3. MTX (Metatranscriptomics) Processing
- **RNA-seq QC**: Quality control and filtering of raw metatranscriptomic reads.
  - **Host/Contaminant Removal**: Filtering out human contaminating sequences using T2T genome assembly and hg38 transcriptome.
  - **rRNA Gene Removal**: Filtering of ribosomal RNA gene reads using SILVA.
- **Taxonomic and Functional Assignment**: Annotation of transcriptionally active microbial taxa (KrakenUniq, Kraken2, Braken, MetaPhlAn4, mOTUs) and microbially expressed genes/transcripts (HUManN.
  - **Paired and Unpaired MGX Processing**: Use paired metagenome taxonomic abundance profiles for processing metatranscriptomic gene expression data.
- **Assembly/Mapping**: Transcript assembly or direct read mapping to references.



Each step above is automated via modular rules in the `rules/` directory, ensuring reproducibility and scalability for large-scale multi-omic studies.

## Installation

1. **Set up the Snakemake environment:**

   ```bash
   mamba create -c conda-forge -c bioconda -p snakemake_2024 snakemake==8.25.2
   mamba activate snakemake_2024
   # (If using modules) module purge all
   mamba install -c conda-forge -c bioconda conda==24.7.1
   ```

2. **Kneaddata Java Fix:**

   If using Kneaddata, activate the Kneaddata conda environment and reinstall Java:

   ```bash
   mamba activate <kneaddata environment>
   mamba install --force-reinstall java-jdk
   mamba deactivate
   ```

## Usage

1. **Prepare your sample sheet:**

   - Assemble a `sample_sheet.tsv` file describing your samples (see `config/test_sample_sheet.tsv` for an example).
   - Place it in the `config/` directory or specify its path on the command line.

2. **Run the pipeline:**

   - **Basic run:**

     ```bash
     snakemake --use-conda --executor slurm --profile cluster_profile -k
     ```

   - **Generate workflow graphs:**

     - Rule graph:
       ```bash
       module load graphviz
       snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
       ```
     - File graph:
       ```bash
       snakemake --forceall --filegraph | dot -Tpdf > dag.pdf
       ```
     - DAG:
       ```bash
       snakemake --forceall --dag | dot -Tpdf > dag.pdf
       ```

   - **Change output directory or configuration:**
     ```bash
     snakemake --use-conda --executor slurm --profile cluster_profile -k --config RESULTS='results/'
     ```


## Main Software Packages

- **MetaPhlAn** (version 4.1.1): Taxonomic profiling
- **HUMAnN**: Functional profiling
- **Kneaddata**: Quality control
- **KrakenUniq**: Taxonomic profiling (good for MTX) 
- **Kraken2/Braken**: Taxonomic profiling
- **mOTUS**: Taxonomic profiling
- **BAQLaVa**: Viral taxonomic profiling
- **MEGAHIT**: Metagenome assembly
- **Snakemake**: Workflow management

See workflow/envs for installation versions.
## Example Sample Sheet

See [`config/test_sample_sheet.tsv`](https://github.com/jtsumner/mlm2/blob/main/config/test_sample_sheet.tsv):

Replace ... with paths to relevant files. Leave blank if file does not exist. This pipeline will automatically subset to each omic-type for relevant processing steps.  

```tsv
sample	AMP_read1	AMP_read2	MGX_read1	MGX_read2	MTX_read1	MTX_read2
Zymo	...	...	...	...	...	...
Zymo2	...	...	...	...	...	...
```


## License

This project is licensed under the MIT License.

---

For more details, see the full code and workflow at the [mlm2 repository](https://github.com/jtsumner/mlm2).  
To see more files and search results, visit the [GitHub code search interface](https://github.com/jtsumner/mlm2/search).
