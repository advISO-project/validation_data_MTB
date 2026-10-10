# Contamination Validation Dataset

## Objective

Assess the workflow's ability to identify and appropriately handle samples containing known levels and types of non-target contamination while maintaining accurate *Mycobacterium tuberculosis* genomic results within the evaluated operating range.

This validation objective focuses specifically on contamination by non-target DNA. MTB–MTB mixtures are evaluated separately under the mixed-infection validation objective.

---

## Source Dataset

The contamination validation dataset is derived from:

**Heupink TH, Verboven L, Warren RM, Van Rie A. (2021).**
*Comprehensive and accurate genetic variant identification from contaminated and low-coverage Mycobacterium tuberculosis whole genome sequencing data.*
Microbial Genomics, 7(11):000689.

- BioProject: `PRJNA706121`
- Study design: simulated MTB whole-genome sequencing data
- Sequencing platform: Illumina
- Public sequencing records retrieved from ENA
- Supplementary Table 1 was used to establish the simulation truth

The study contains 1,200 simulated sequencing datasets:

- 600 simulated culture samples
- 600 simulated direct-from-sputum samples

For the initial contamination validation dataset, only the simulated culture samples at 50× MTB coverage are used.

This separates contamination effects from low sequencing coverage, which should be evaluated as a separate validation variable.

---

## Simulation Truth

The simulated culture dataset was generated from MTB strains with known underlying genomic variation.

Each simulated MTB strain was derived from H37Rv and contained:

- 1,000 SNPs
- 50 single-nucleotide insertions
- 50 single-nucleotide deletions

Because the MTB sequence and introduced variants were defined during simulation, the expected MTB genomic content is known independently of the workflow being validated.

Culture simulations were generated at:

- 5×
- 10×
- 20×
- 30×
- 50×
- 100×

The contamination validation dataset is restricted to **50×**.

---

## Contamination Truth

For each simulated sample, the source metadata identifies:

- simulated MTB strain
- MTB coverage
- contaminant identity
- contamination percentage
- simulation seed

Contamination in the simulated culture dataset ranges from **0% to 5%**.

The contaminant codes represented in the source dataset are:

| Code | Contaminant | Class |
|---|---|---|
| `None` | No contamination | None |
| `Hsap` | *Homo sapiens* | Human |
| `Mabs` | *Mycobacterium abscessus* | NTM |
| `Mavi` | *Mycobacterium avium* | NTM |
| `Mint` | *Mycobacterium intracellulare* | NTM |
| `Paer` | *Pseudomonas aeruginosa* | Bacterial |
| `Sepi` | *Staphylococcus epidermidis* | Bacterial |
| `Mix1` | NTM mixture | Mixed NTM |
| `Mix2` | All-contaminant mixture | Mixed |

The contamination percentage and contaminant identity are treated as ground truth because they were defined during simulation rather than inferred from the output of the workflow under validation.

---

## Validation Endpoints

The contamination dataset can be used to evaluate three related aspects of workflow performance.

### 1. Contamination Detection

Where the workflow reports contamination, assess whether contamination is detected and whether the reported contamination level/type is consistent with the known simulation truth.

### 2. QC Handling

Assess whether contaminated samples are assigned the expected workflow QC status according to predefined workflow acceptance criteria.

The expected `PASS`, `WARN`, or `FAIL` behaviour must be defined from the requirements of the workflow being validated and must not be derived from the behaviour observed during validation.

### 3. MTB Result Integrity

Assess whether the workflow continues to recover the expected MTB genomic results in the presence of non-target contamination.

Depending on the workflow being validated, this may include:

- variant calls
- lineage
- antimicrobial-resistance predictions
- other MTB-specific outputs

The known underlying simulated MTB sequence provides the reference truth for these comparisons.

---

## Master Truth Dataset

Supplementary Table 1 from Heupink et al. was parsed programmatically.

A total of **600 simulated culture sample identifiers** were recovered.

The dataset was then restricted to simulations with:

- dataset coverage = 50×
- MTB coverage = 50×

This produced a master contamination truth set containing:

**100 samples**

The master truth set is stored at:

```text
contamination/reference_dataset/contamination_truth_set.csv
```

The truth set contains fields describing:

- source simulation identifier
- simulated MTB strain
- dataset coverage
- MTB coverage
- contaminant code
- contaminant identity
- contaminant class
- contamination percentage
- simulation seed
- source study
- BioProject

---

## Mapping to Public Sequencing Data

Run-level metadata for BioProject `PRJNA706121` was retrieved from the European Nucleotide Archive (ENA).

A total of:

**1,200 sequencing runs**

were retrieved.

The ENA `sample_alias` contains the simulation identifier used in the published dataset. This enables an exact mapping between the independently established simulation truth and the publicly available sequencing run.

The mapping is therefore:

```text
Published simulation truth
        ↓
source_sample_id
        =
ENA sample_alias
        ↓
SRR run accession
        ↓
paired FASTQ files
```

All 100 COV50 truth records mapped uniquely to ENA:

- 100 truth samples
- 100 ENA matches
- 100 unique run accessions
- 100 paired-end sequencing datasets
- no unmatched truth records
- no duplicate run mappings

The complete mapped dataset is stored at:

```text
contamination/reference_dataset/contamination_validation_manifest.csv
```

The manifest also retains the ENA FASTQ locations, file sizes and MD5 checksums.

---

## Pilot Validation Panel

The 100-sample mapped dataset is retained as the **master validation resource**.

A smaller 30-sample pilot validation panel was selected programmatically for the initial validation exercise.

The pilot contains:

| Contamination | Samples |
|---:|---:|
| 0% | 5 |
| 1% | 5 |
| 2% | 5 |
| 3% | 5 |
| 4% | 5 |
| 5% | 5 |
| **Total** | **30** |

The five 0% samples serve as uncontaminated controls.

For contamination levels 1–5%, selection prioritises:

1. contaminant-code diversity within each contamination level;
2. balanced representation of contaminant codes across the complete panel;
3. deterministic selection using `source_sample_id` to resolve ties.

No samples were selected manually by SRR accession.

The selected contaminated panel contains all eight contaminant codes represented in the master dataset.

The contaminant representation in the selected pilot panel is:

| Contaminant code | Selected samples |
|---|---:|
| `Hsap` | 3 |
| `Mabs` | 4 |
| `Mavi` | 3 |
| `Mint` | 4 |
| `Mix1` | 3 |
| `Mix2` | 2 |
| `Paer` | 3 |
| `Sepi` | 3 |

In addition, five `None` samples are included as uncontaminated controls.

Because contaminant types are not represented equally at every contamination percentage in the source dataset, a complete contaminant × percentage factorial design was not imposed.

The selection policy is documented in:

```text
contamination/metadata/contamination_selection_policy.md
```

The resulting pilot manifest is stored at:

```text
contamination/reference_dataset/contamination_pilot_manifest.csv
```

The full 100-sample master manifest remains unchanged and can be used to expand the validation dataset if additional samples are required.

---

## Dataset Construction Workflow

The contamination validation resource is generated using the following scripts.

### 01 — Build Truth Set

```text
scripts/contamination/01_build_truth_set.py
```

Parses the published supplementary workbook, extracts simulation identifiers and constructs the 100-sample COV50 contamination truth set.

### 02 — Retrieve Sequencing Metadata

```text
scripts/contamination/02_retrieve_sra_metadata.py
```

Retrieves run-level metadata for `PRJNA706121` from ENA.

### 03 — Build Validation Manifest

```text
scripts/contamination/03_build_validation_manifest.py
```

Maps the 100 truth-set samples to their corresponding public sequencing runs using:

```text
source_sample_id = ENA sample_alias
```

### 04 — Build Pilot Manifest

```text
scripts/contamination/04_build_pilot_manifest.py
```

Applies the documented selection policy to construct the reproducible 30-sample pilot panel.

### 05 — Download FASTQs

```text
scripts/contamination/05_download_fastqs.py
```

Downloads the paired FASTQ files specified by the pilot manifest and verifies them against the ENA MD5 checksums.

The download script is restart-safe. Existing FASTQ files are checked against their expected MD5 checksum before being skipped, and downloads are initially written as temporary `.part` files before finalisation.

---

## Generated Files

The contamination validation dataset contains both original source material and files generated programmatically during dataset construction.

### Source and manually maintained files

#### `metadata/heupink2021_supplementary_table1_culture_simulations.xlsx`

Original supplementary workbook from Heupink et al. (2021).

This file is retained unchanged and serves as the source of the published simulation truth used to construct the contamination truth set.

#### `metadata/contamination_selection_policy.md`

Documents the predefined rules used to select the 30-sample pilot validation panel from the 100-sample master validation resource.

This file defines the selection strategy and is not generated by the selection script.

### Generated metadata

#### `metadata/contamination_sra_metadata.tsv`

Generated by:

```text
scripts/contamination/02_retrieve_sra_metadata.py
```

Contains run-level ENA metadata retrieved for BioProject `PRJNA706121`, including sample aliases, run accessions, sequencing metadata, FASTQ locations, file sizes and MD5 checksums.

### Generated reference datasets

#### `reference_dataset/contamination_truth_set.csv`

Generated by:

```text
scripts/contamination/01_build_truth_set.py
```

Contains the 100 simulated culture samples at 50× coverage selected from the published supplementary simulation dataset.

This file represents the master contamination truth set before linkage to public sequencing accessions.

#### `reference_dataset/contamination_validation_manifest.csv`

Generated by:

```text
scripts/contamination/03_build_validation_manifest.py
```

Contains the complete 100-sample master validation dataset after mapping each simulation truth record to its corresponding ENA sequencing run.

The manifest retains simulation truth, ENA/SRA accessions, paired FASTQ locations, file sizes and MD5 checksums.

#### `reference_dataset/contamination_pilot_manifest.csv`

Generated by:

```text
scripts/contamination/04_build_pilot_manifest.py
```

Contains the reproducibly selected 30-sample pilot validation panel.

The complete master validation manifest remains unchanged.

### Generated reports

#### `reports/contamination_truth_set_report.txt`

Generated by `01_build_truth_set.py`.

Summarises truth-set construction, including the number of simulations recovered, COV50 samples retained, contamination-level distribution and contaminant-code distribution.

#### `reports/contamination_sra_metadata_report.txt`

Generated by `02_retrieve_sra_metadata.py`.

Summarises ENA metadata retrieval, including the number of runs retrieved, accession uniqueness, sequencing layout and metadata completeness.

#### `reports/contamination_validation_manifest_report.txt`

Generated by `03_build_validation_manifest.py`.

Documents the mapping between the 100 truth-set records and ENA sequencing runs and reports mapping completeness and contaminant distributions.

#### `reports/contamination_pilot_manifest_report.txt`

Generated by `04_build_pilot_manifest.py`.

Documents the resulting 30-sample pilot panel, including contamination-level representation, contaminant-code representation and the selected sample/run accessions.

#### `reports/contamination_fastq_download_report.txt`

Generated by `05_download_fastqs.py` after successful FASTQ retrieval.

Records the download and checksum-verification status of each FASTQ file.

### Downloaded sequencing data

#### `fastqs/*.fastq.gz`

Downloaded by:

```text
scripts/contamination/05_download_fastqs.py
```

Contains the paired-end FASTQ files corresponding to the samples in `contamination_pilot_manifest.csv`.

Each downloaded file is verified against the MD5 checksum supplied by ENA before being accepted as part of the validation dataset.

FASTQ download and checksum verification are currently pending.

---


## Directory Structure

```text
contamination/
├── README.md
├── metadata/
│   ├── heupink2021_supplementary_table1_culture_simulations.xlsx
│   ├── contamination_sra_metadata.tsv
│   └── contamination_selection_policy.md
├── reference_dataset/
│   ├── contamination_truth_set.csv
│   ├── contamination_validation_manifest.csv
│   └── contamination_pilot_manifest.csv
├── fastqs/
└── reports/
    ├── contamination_truth_set_report.txt
    ├── contamination_sra_metadata_report.txt
    ├── contamination_validation_manifest_report.txt
    ├── contamination_pilot_manifest_report.txt
    └── contamination_fastq_download_report.txt
```

Associated scripts are stored under:

```text
scripts/contamination/
```

---

## Reproducibility and Provenance

The original supplementary workbook is retained unchanged under `metadata/`.

Derived truth sets and validation manifests are generated programmatically rather than edited manually.

The validation resource maintains the following provenance chain:

```text
Published study
      ↓
Supplementary simulation metadata
      ↓
Contamination truth set
      ↓
ENA run metadata
      ↓
Master validation manifest
      ↓
Documented pilot-selection policy
      ↓
Pilot validation manifest
      ↓
Checksum-verified FASTQ files
      ↓
Workflow results
      ↓
Observed vs expected comparison
```

This ensures that expected results are established independently of the workflow being validated and that the origin and selection of each validation sample can be traced.

---

## Current Status

Completed:

- [x] Source study identified
- [x] Supplementary simulation metadata obtained
- [x] 600 culture simulations parsed
- [x] 100 COV50 samples selected as master truth dataset
- [x] ENA metadata retrieved
- [x] Truth samples mapped to public sequencing runs
- [x] Master validation manifest generated
- [x] Pilot selection policy defined
- [x] 30-sample pilot manifest generated

Pending:

- [ ] Confirm sufficient storage for FASTQ dataset
- [ ] Download pilot FASTQs
- [ ] Verify all FASTQ MD5 checksums
- [ ] Run workflow
- [ ] Extract contamination, QC and MTB genomic outputs
- [ ] Compare observed results with expected truth
- [ ] Define and calculate validation performance metrics
- [ ] Generate final validation report