# BactiPipe  

**BactiPipe** was developed and is actively used by the **NGS section of the Penn State Animal Diagnostic Laboratory** to bring **consistency** to the analysis of **BACTERIAL** whole-genome sequencing (WGS) data. Multiple tools and parameters are available for genomic data analysis but this flexibility often leads to variability in results from run to run. In diagnostic laboratories, such variability may compromize one of the most critical aspects of quality control in a diagnostic setting: **reproducibility and standardization**. BactiPipe addresses this challenge by wrapping together the core tools used in our laboratory with carefully chosen and fixed parameters, ensuring consistent analysis across all sequencing runs.  

BactiPipe provides **end-to-end quality control and genome assembly workflows** for both **Illumina** and **Oxford Nanopore** sequencing data. It automates read QC, genome assembly, organism identification, and assembly quality checks using tools such as **fastp, filtlong, Spades, Unicycler, Flye, KmerFinder, Kraken2,** and **CheckM**. The pipeline supports data from **local directories, network-mounted storage, and Amazon S3**, making it suitable for both on-premises and cloud-based environments. Outputs include detailed logs, quality metrics plots, assembled genomes, and comprehensive QC summary reports.  

For a detailed description of the pipeline, workflows, and usage examples, please see the [**BactiPipe Wiki Page**](https://github.com/ADL-genomics/BactiPipe/wiki).  

The current CLI exposes four analysis commands:

- `bactipipe qc-illumina` and `bactipipe qc-nanopore` perform read QC, assembly,
  identification, and assembly QC.
- `bactipipe relate` performs organism-appropriate typing plus optional ANI,
  SKA, and cgMLST analysis on a related cohort of assembled genomes.
- `bactipipe detect` reports AMR determinants and virulence genes from assembled
  genomes.

Run `bactipipe <command> --help` for the current input contract and options.

For the `relate`/`detect` helper environment, create the compatible base and
then install the exact legacy CGE/SeqSero wrappers without their obsolete
dependency metadata:

```bash
mamba env create -f envs/genepid.yml
bash envs/install_genepid_tools.sh genepid
mamba env create -f envs/viramr.yml
mamba env create -f envs/kleborate.yml
bash envs/install_kleborate_compat.sh kleborate
```

Install or refresh the databases used by downstream typing and detection with:

```bash
bactipipe update-db mlst serotypefinder virulencefinder cgmlstfinder amrfinder abricate
```

The former CGE `mlst_db` Bitbucket repository has been removed. BactiPipe therefore rebuilds the
same directory format expected by CGE `mlst.py` from the official PubMLST and Institut Pasteur BIGSdb
REST APIs. Source endpoints, dates, loci, and server access notices are saved in
`mlst_db/bactipipe_mlst_sources.json`. If authenticated BIGSdb access has been issued to you, set the
complete HTTP Authorization header in the host-specific `BACTIPIPE_PUBMLST_AUTHORIZATION` or
`BACTIPIPE_PASTEUR_AUTHORIZATION` variable so credentials are sent only to the issuing service.
