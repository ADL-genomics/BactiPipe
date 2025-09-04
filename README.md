**BactiPipe: Quality Assessment and Control for Bacterial Whole Genome Data**

BactiPipe is a pipeline that integrates multiple tools to evaluate and improve the quality of raw sequencing reads from Illumina and Nanopore platforms. It assesses read quality based on user-defined thresholds for average Phred scores and genome coverage depth. If necessary, the pipeline applies corrective measures to enhance data quality.

Reads that pass the quality control (QA/QC) step are assembled into genome sequences, and the resulting assemblies undergo further evaluation. This assessment includes verifying taxonomic consistency and determining assembly completeness and contamination using CheckM.

This page provides a detailed overview of the pipeline's steps and usage instructions.
