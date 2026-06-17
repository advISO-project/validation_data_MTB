# validation_tuberculosis

Validation resources for testing and validating *Mycobacterium tuberculosis* (MTB) bioinformatics workflows as part of the advISO (Advancing Bioinformatics through ISO Accreditation) project.

## Background

ISO standards such as ISO 15189 and ISO 17025 require laboratories to demonstrate that analytical methods are fit for purpose. While these standards are well established in clinical and diagnostic laboratories, practical resources for validating pathogen bioinformatics workflows remain limited.

The aim of this repository is to provide a reproducible validation framework for MTB genomic analysis workflows. This includes curated validation datasets, truth sets, validation procedures, and workflow documentation that can be used to assess workflow performance and support quality assurance activities.

## Objectives

The validation framework is designed to assess:

* Antimicrobial resistance (AMR) prediction
* Lineage assignment
* Variant calling
* Species determination
* Detection of mixed or contaminated samples
* End-to-end workflow reproducibility

## Core Validation Components

### AMR Validation

Ensure genotype-to-phenotype mapping accurately predicts resistance.

Validation datasets include major first-line and second-line drugs:

* Rifampicin
* Isoniazid
* Ethambutol
* Pyrazinamide
* Fluoroquinolones
* Bedaquiline
* Linezolid
* Additional drugs where validated data are available

### Lineage Assignment

Validate correct assignment to MTB lineages and sub-lineages.

The framework aims to:

* Verify lineage classification accuracy
* Evaluate sub-lineage assignment
* Detect mixed-lineage samples
* Assess performance across diverse MTB populations

### Variant Calling

Evaluate detection of genomic variation including:

* Single nucleotide polymorphisms (SNPs)
* Insertions and deletions (indels)
* Resistance-associated mutations
* Genome-wide variant detection

### Mixed / Contaminated Samples

Assess workflow robustness when processing:

* Mixed infections
* Low-frequency variants
* Contaminated datasets
* Ambiguous lineage assignments

### Species Determination

Evaluate the ability to distinguish:

* *Mycobacterium tuberculosis*
* Animal-associated members of the MTB complex
* Non-tuberculous mycobacteria (NTMs)

## Repository Structure

```text
workflow/
    Workflow definitions and workflow documentation.

validation_data/
    Sample manifests, accession lists and metadata for validation datasets.

truth_sets/
    Curated reference datasets and expected results used for validation.

docs/
    Validation procedures, reports and supporting documentation.

scripts/
    Utility scripts for dataset generation and validation.

env/
    Conda and pip environments required for reproducing the validation datasets.
```

## Validation Datasets

### AMR Validation Dataset

The AMR validation dataset is derived from publicly available MTB isolates with high-confidence genotype and phenotype information.

Current resources include:

* AMR truth sets
* Sample manifests
* Galaxy collection manifests
* Validation reports
* Example workflow outputs

Primary source:

* tbtAMR (PRJNA857537)

### Lineage Validation Dataset

The lineage validation framework is currently under active development.

Current resources include:

* TBDB lineage barcode markers
* Coll et al. lineage reference datasets
* Candidate lineage validation samples
* Pilot FASTQ validation panel
* TB-Profiler lineage characterization results

Primary sources:

* Coll et al.
* TBDB barcode database

## Trusted Data Sources

| Validation Target     | Source / Dataset               | Notes                                        |
| --------------------- | ------------------------------ | -------------------------------------------- |
| AMR prediction        | tbtAMR (PRJNA857537)           | High-confidence isolates with phenotypic DST |
| Lineage assignment    | Coll et al.; TBDB              | Lineage-defining SNPs and reference genomes  |
| Variant calling       | Public MTB WGS datasets        | Genome-wide SNP validation                   |
| Species determination | MTB and NTM reference datasets | Evaluation of species classification         |
| Mixed samples         | Under development              | Future validation component                  |

## Dataset Generation Workflow

The repository contains scripts used to construct validation datasets from public resources.

### 1. Extract Truth Set

`extract_samples_from_tbtamr.py`

* Reads source validation datasets
* Selects isolates for validation
* Produces curated truth sets

### 2. Download FASTQs

`extract_fastqs_for_tbtamr.py`

* Retrieves paired-end FASTQ files from ENA
* Downloads validation reads
* Maintains reproducible sample selection

### 3. Generate Validation Manifests

`generate_manifest.py`

Produces:

* `truth_set_manifest.csv`
* `validation_set_manifest.yaml`

These manifests can be used directly within Galaxy or Planemo workflows.

### 4. Optional End-to-End Pipeline

`run_pipeline.py`

Executes the complete validation dataset generation process.

## Validation Workflow

The general validation process is:

1. Obtain validation samples using the supplied manifests.
2. Execute the workflow under evaluation.
3. Compare workflow outputs against the supplied truth sets.
4. Calculate concordance and performance metrics.
5. Document validation results according to the provided procedures.

## Environment Setup

### Using Conda

```bash
conda env create -f env/environment.yml
conda activate validation_env
```

### Using Pip

```bash
pip install -r env/requirements.txt
```

## Current Status

### Available

* AMR validation dataset
* Validation manifests
* Dataset generation scripts
* TBDB lineage marker resources
* Pilot lineage validation dataset
* TB-Profiler lineage characterization

### Under Development

* Expanded lineage validation panel
* Indel validation framework
* Mixed infection validation datasets
* Species determination validation datasets
* Galaxy workflow publication through the Intergalactic Workflow Commission (IWC)

## Future Work

Planned additions include:

* Broader lineage coverage 
* Mixed infection validation datasets
* Contamination assessment datasets
* Additional pathogen-specific validation resources
* Automated validation reporting
* Workflow publication through Galaxy community infrastructure

## advISO Project

The advISO project aims to develop practical resources that support quality-assured pathogen bioinformatics and facilitate adoption of ISO accreditation principles within bioinformatics workflows.

Project website:

https://www.cardiff.ac.uk/adviso-bioinformatics-accreditation

## Citation

If you use these validation resources, please cite the original data sources as described in the accompanying documentation and metadata files.
