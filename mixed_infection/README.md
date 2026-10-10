# Mixed-Infection Validation Dataset

## Overview

This directory contains the mixed-infection validation dataset developed for the ISO-in-a-Box validation framework.

The dataset is intended to support validation of bioinformatics workflows that detect mixed *Mycobacterium tuberculosis* infections and, where applicable, estimate the relative proportions of constituent strains.

The reference dataset consists of 48 previously characterised sequencing samples with known composition:

- 36 mixed samples
- 12 non-mixed (pure) controls

The mixed samples represent three known major/minor strain proportions:

- 70/30
- 90/10
- 95/5

The 12 pure controls have a major-strain proportion of 1.00 and an expected minor-strain proportion of 0.00.

Expected results are assigned from the documented known sample composition rather than from the output of a bioinformatics classification method.

---

## 1. Source Dataset

### 1.1 Primary Truth Source

The sample-level expected composition was obtained from Table 1 of:

> Sobkowiak B, et al.
> *Identifying mixed Mycobacterium tuberculosis infections from whole genome sequence data.*
> BMC Genomics. 2018;19:613.
> DOI: 10.1186/s12864-018-4988-z
> PMCID: PMC6092779
> PMID: 30107785

Table 1 reports the characterised samples together with their known major-strain proportions and constituent strain information.

The published table contains 48 samples distributed as follows:

| Known major proportion | Expected minor proportion | Expected status | Samples |
|---|---:|---|---:|
| 0.70 | 0.30 | Mixed | 12 |
| 0.90 | 0.10 | Mixed | 12 |
| 0.95 | 0.05 | Mixed | 12 |
| 1.00 | 0.00 | Pure | 12 |
| **Total** | | | **48** |

The known sample composition is used as the basis for expected results.

Computational results reported in the 2018 study, including heterozygous-SNP measurements, inferred numbers of strains, and computationally estimated major-strain proportions, are not used to define truth in this validation dataset.

### 1.2 Supporting Validation Study

The same characterised cohort was subsequently used in:

> Sobkowiak B, et al.
> *A new method for detecting mixed Mycobacterium tuberculosis infection and reconstructing constituent strains provides insights into transmission.*
> Genome Medicine. 2025;17:8.
> DOI: 10.1186/s13073-025-01430-y
> PMCID: PMC11771024
> PMID: 39871355

The 2025 study used the 36 mixed samples and 12 non-mixed constituent strains to evaluate MixInfect2.

The MixInfect2 results are treated as supporting evidence for the suitability of this characterised cohort for mixed-infection validation. They are not used to assign the expected mixed/pure status in the ISO-in-a-Box truth dataset.

---

## 2. Truth Assignment

Truth assignment is based on the independently known composition of the characterised samples.

For each sample, the following information is retained from the published source:

- run accession
- known major-strain proportion
- major strain lineage
- minor strain lineage, where applicable
- major strain spoligotype
- minor strain spoligotype, where applicable

The expected minor-strain proportion is derived as:

```text
expected minor proportion = 1 - published major proportion
```

Expected mixed-infection status is assigned according to the documented composition:

```text
major proportion 0.70 -> mixed
major proportion 0.90 -> mixed
major proportion 0.95 -> mixed
major proportion 1.00 -> pure
```

This produces:

```text
36 mixed samples
12 pure controls
48 total samples
```

For pure samples, minor-lineage and minor-spoligotype fields are left blank rather than assigning an artificial second constituent.

The resulting truth table is:

```text
reference_dataset/mixed_infection_truth_set.csv
```

This file represents the expected results against which workflow outputs can be evaluated.

---

## 3. Source Freezing and Provenance

The published Table 1 used for truth construction is frozen locally as:

```text
raw/sobkowiak2018_table1.html
```

Source provenance is recorded in:

```text
metadata/sobkowiak2018_table1_source.txt
```

The frozen source is protected by a SHA256 checksum recorded in:

```text
metadata/sobkowiak2018_table1.sha256
```

The truth dataset is generated deterministically from this frozen source using:

```text
../scripts/mixed_infection/02_build_truth_set.py
```

The generated truth dataset is also checksummed:

```text
metadata/mixed_infection_truth_set.sha256
```

This allows the truth dataset to be reconstructed from the published source while preserving the provenance of the original evidence.

---

## 4. Public Sequence Data

The 48 run accessions represented in the truth dataset were individually verified against the European Nucleotide Archive (ENA).

All 48 accessions resolved successfully.

The verified records have the following characteristics:

```text
Study accession:       PRJEB2794
Run accessions:        48
Library layout:        PAIRED
Instrument platform:   ILLUMINA
Instrument model:      Illumina HiSeq 2000
FASTQ links available: 48/48
```

ENA metadata is stored in:

```text
metadata/mixed_infection_ena_metadata.tsv
```

and checksummed in:

```text
metadata/mixed_infection_ena_metadata.sha256
```

The ENA metadata provides sequence-data provenance and retrieval information. It is not used to assign expected mixed-infection status.

---

## 5. Validation Manifest

The truth dataset and ENA metadata are joined using the run accession as the common identifier.

The resulting validation manifest is:

```text
reference_dataset/mixed_infection_validation_manifest.csv
```

The manifest contains one row per sample and connects the expected result to the corresponding public sequencing data.

### Expected-result information

- expected mixed/pure status
- known major-strain proportion
- expected minor-strain proportion
- major lineage
- minor lineage
- major spoligotype
- minor spoligotype
- truth basis
- primary truth source

### Sequence-data identity

- ENA study accession
- sample accession
- experiment accession
- run accession
- scientific name
- library layout
- sequencing platform
- instrument model

### FASTQ provenance

- read 1 FASTQ location
- read 2 FASTQ location
- ENA MD5 checksum for each FASTQ
- expected file size for each FASTQ

Validation of the joined manifest confirmed:

```text
Truth records:                 48
ENA records:                   48
Joined records:                48
Missing accessions:             0
Unexpected accessions:          0

Mixed samples:                 36
Pure controls:                 12

0.70 major proportion:         12
0.90 major proportion:         12
0.95 major proportion:         12
1.00 major proportion:         12

Unique sample accessions:      48
Unique experiment accessions:  48

Paired FASTQ sets:             48
FASTQ files represented:       96
Expected total FASTQ size:     approximately 23 GiB
```

The manifest is checksummed in:

```text
metadata/mixed_infection_validation_manifest.sha256
```

---

## 6. FASTQ Retrieval

FASTQ files are not required to construct or inspect the truth dataset.

The validation manifest contains the public locations, expected file sizes and ENA MD5 checksums required to retrieve the sequence data when workflow validation is performed.

A retrieval script is provided:

```text
../scripts/mixed_infection/05_download_fastqs.py
```

The downloader is designed to:

- retrieve the paired FASTQs represented in the validation manifest
- verify downloaded files against ENA MD5 checksums
- retain already verified files when restarted
- identify incomplete or invalid files
- produce a per-file retrieval report

The complete dataset represents approximately 23 GiB of compressed FASTQ data.

FASTQ retrieval is therefore performed on demand rather than being required during reference-dataset construction.

Downloaded sequence files are stored under:

```text
fastqs/
```

The `fastqs/` directory may therefore remain empty until the sequence data are required for workflow execution.

---

## 7. Construction Workflow

The mixed-infection reference dataset is constructed using the following workflow:

```text
Sobkowiak et al. 2018 Table 1
              |
              v
01_retrieve_table1.py
              |
              v
Frozen published Table 1
              |
              v
02_build_truth_set.py
              |
              v
mixed_infection_truth_set.csv
              |
              v
03_retrieve_ena_metadata.py
              |
              v
ENA metadata for exact 48 run accessions
              |
              v
04_build_validation_manifest.py
              |
              v
mixed_infection_validation_manifest.csv
              |
              v
05_download_fastqs.py
              |
              v
FASTQs retrieved on demand when required
```

### `01_retrieve_table1.py`

Retrieves and freezes the published Table 1 source and records its provenance and SHA256 checksum.

### `02_build_truth_set.py`

Extracts the known-composition fields from the frozen source, derives expected minor proportions and mixed/pure status, and constructs the canonical truth dataset.

The script requires:

```text
48 unique run accessions
12 samples at major proportion 0.70
12 samples at major proportion 0.90
12 samples at major proportion 0.95
12 samples at major proportion 1.00
36 expected mixed samples
12 expected pure samples
```

Failure of these integrity checks prevents generation of a valid truth dataset.

### `03_retrieve_ena_metadata.py`

Retrieves ENA metadata for the exact 48 run accessions represented in the truth dataset.

The complete ENA project is not used as the starting point for candidate discovery.

### `04_build_validation_manifest.py`

Performs a one-to-one join between the truth dataset and ENA metadata and validates the resulting sequence-data mappings.

### `05_download_fastqs.py`

Retrieves the sequence files represented in the validated manifest and verifies them against the MD5 checksums supplied by ENA.

FASTQ retrieval is optional until the dataset is required for workflow execution.

---

## 8. Directory Structure

```text
mixed_infection/
├── README.md
├── fastqs/
├── metadata/
│   ├── sobkowiak2018_table1_source.txt
│   ├── sobkowiak2018_table1.sha256
│   ├── mixed_infection_truth_set.sha256
│   ├── mixed_infection_ena_metadata.tsv
│   ├── mixed_infection_ena_metadata.sha256
│   └── mixed_infection_validation_manifest.sha256
├── raw/
│   └── sobkowiak2018_table1.html
├── reference_dataset/
│   ├── mixed_infection_truth_set.csv
│   └── mixed_infection_validation_manifest.csv
└── reports/
    ├── mixed_infection_truth_set_report.txt
    ├── mixed_infection_ena_retrieval_report.txt
    └── mixed_infection_validation_manifest_report.txt
```

Associated construction scripts are stored under:

```text
scripts/mixed_infection/
├── 01_retrieve_table1.py
├── 02_build_truth_set.py
├── 03_retrieve_ena_metadata.py
├── 04_build_validation_manifest.py
└── 05_download_fastqs.py
```

---

## 9. Interpretation of the Dataset Components

The mixed-infection validation resource deliberately separates three concepts.

### 9.1 Truth

```text
reference_dataset/mixed_infection_truth_set.csv
```

This defines the expected result from the independently known sample composition.

It answers:

> What result is expected for this sample?

### 9.2 Sequence-Data Provenance

```text
metadata/mixed_infection_ena_metadata.tsv
```

This identifies the public sequencing records corresponding to the samples represented in the truth dataset.

It answers:

> Which public sequencing data correspond to this sample?

### 9.3 Validation Manifest

```text
reference_dataset/mixed_infection_validation_manifest.csv
```

This links the expected result to the exact sequence data required to test a workflow.

It answers:

> Which sequence files should be analysed, and what result should be expected?

This separation prevents computational predictions from being used to establish the expected result against which another computational workflow is evaluated.

---

## 10. Current Status

The mixed-infection reference dataset has been constructed and linked to public sequencing records.

```text
Published source frozen:              COMPLETE
Source checksum recorded:             COMPLETE

Truth dataset constructed:            COMPLETE
Truth-set integrity checks:           PASS
Samples in truth dataset:             48
Mixed samples:                        36
Pure controls:                        12

ENA accessions resolved:              48/48
ENA metadata retrieval:               PASS

Validation manifest constructed:      COMPLETE
Validation manifest integrity checks: PASS
Paired FASTQ sets represented:        48
FASTQ files represented:              96

FASTQ retrieval:                      ON DEMAND / DEFERRED
```

The reference dataset is therefore ready for use in the design and execution of mixed-infection workflow validation.

---

## 11. Key Files

| File | Purpose |
|---|---|
| `raw/sobkowiak2018_table1.html` | Frozen copy of the published source table |
| `metadata/sobkowiak2018_table1_source.txt` | Source and retrieval provenance |
| `metadata/sobkowiak2018_table1.sha256` | Integrity checksum for frozen source |
| `reference_dataset/mixed_infection_truth_set.csv` | Canonical expected-result truth table |
| `metadata/mixed_infection_truth_set.sha256` | Integrity checksum for truth table |
| `metadata/mixed_infection_ena_metadata.tsv` | ENA metadata for the 48 truth-set accessions |
| `metadata/mixed_infection_ena_metadata.sha256` | Integrity checksum for ENA metadata |
| `reference_dataset/mixed_infection_validation_manifest.csv` | Joined truth and sequence-data manifest |
| `metadata/mixed_infection_validation_manifest.sha256` | Integrity checksum for validation manifest |
| `reports/mixed_infection_truth_set_report.txt` | Truth-set construction and QC report |
| `reports/mixed_infection_ena_retrieval_report.txt` | ENA retrieval and QC report |
| `reports/mixed_infection_validation_manifest_report.txt` | Validation-manifest construction and QC report |

---

## 12. Reproducibility

The canonical mixed-infection validation resource can be reconstructed by running the scripts in numerical order:

```bash
python scripts/mixed_infection/01_retrieve_table1.py
python scripts/mixed_infection/02_build_truth_set.py
python scripts/mixed_infection/03_retrieve_ena_metadata.py
python scripts/mixed_infection/04_build_validation_manifest.py
```

FASTQ retrieval is performed separately, only when sequence files are required:

```bash
python scripts/mixed_infection/05_download_fastqs.py
```

The frozen source files, checksums, deterministic truth-assignment rules, ENA accession mappings and validation manifest together provide the provenance required to reconstruct and audit the mixed-infection validation dataset.
