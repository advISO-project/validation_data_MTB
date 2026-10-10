# tbtAMR AMR Validation Reference Dataset

## Purpose

This directory contains the antimicrobial resistance (AMR) reference dataset constructed for validation of MTB bioinformatics workflows within the advISO / ISO-IN-A-BOX validation framework.

The dataset is derived from the tbtAMR reference data and currently focuses on the four first-line anti-tuberculosis drugs:

* Rifampicin (RIF)
* Isoniazid (INH)
* Pyrazinamide (PZA)
* Ethambutol (EMB)

Phenotypic drug susceptibility testing (DST) is used as the expected AMR classification. Genotypic resistance mechanisms and the original tbtAMR classification results are retained as supporting provenance and confidence information rather than being used independently to define truth.

---

## Dataset Construction

The current AMR dataset was constructed from:

`validation_results_tbtamr_master.csv`

The dataset-construction workflow is implemented in:

`../scripts/amr/`

The workflow consists of three main steps:

```text
Source tbtAMR dataset
        │
        ▼
01 — Construct high-confidence cohort and reference panel
        │
        ▼
02 — Retrieve and verify paired-end FASTQ files
        │
        ▼
03 — Assemble and verify the final validation dataset
```

---

## Step 01 — Construct the Reference Panel

**Script:**

`../scripts/amr/01_extract_samples_from_tbtamr.py`

The source dataset contains **3,004 isolates**.

The current implementation:

1. Restricts the source data to records from `Current study`.
2. Includes the four first-line drugs RIF, INH, PZA and EMB.
3. Requires a definitive phenotypic DST result (`R` or `S`) for every included drug.
4. Requires the corresponding tbtAMR source classification to be high-confidence (`TP` or `TN`) for every included drug.
5. Retains isolates meeting these requirements as the complete high-confidence eligible cohort.
6. Classifies isolates according to the number of phenotypically resistant first-line drugs.
7. Selects a deterministic reference panel using random seed `42`.

The resulting high-confidence cohort contains **1,805 isolates**.

The current reference panel contains **10 isolates**:

* 3 fully susceptible
* 2 mono-resistant
* 5 poly-resistant

**Outputs:**

`reference_dataset/tbtamr_current_study_high_confidence.csv`

Complete high-confidence population from which the validation panel is selected.

`reference_dataset/tbtamr_current_study_reference_panel.csv`

Final selected 10-isolate reference panel. This table contains the expected phenotypic DST classifications together with resistance-category and supporting provenance information.

`reference_dataset/tbtamr_current_study_reference_panel_report.txt`

Human-readable record of the dataset filtering and selection parameters, eligible population counts and final panel composition.

---

## Step 02 — Retrieve Sequence Data

**Script:**

`../scripts/amr/02_download_fastqs_tbtamr.py`

Paired-end sequencing reads are retrieved for the 10 isolates selected during Step 01.

The script verifies that paired FASTQ files are present and non-empty and records the retrieval status for each selected isolate.

**Output:**

`reference_dataset/tbtamr_current_study_fastq_manifest.csv`

The FASTQ retrieval manifest records the sequencing files associated with each selected isolate together with file information and retrieval/verification status.

The FASTQ files themselves are stored in:

`fastqs/`

The current validation dataset contains:

* 10 isolates
* 20 FASTQ files
* 10 R1/R2 pairs

FASTQ files are excluded from Git because of their size and can be retrieved again using the download script.

---

## Step 03 — Assemble the Validation Dataset

**Script:**

`../scripts/amr/03_generate_manifest.py`

This step combines the selected reference panel from Step 01 with the verified FASTQ manifest from Step 02.

Before generating the final manifests, the script verifies that:

* the reference panel and FASTQ manifest contain exactly the same samples;
* sample identifiers are unique;
* each selected isolate has both an R1 and R2 FASTQ;
* the FASTQ files exist and are non-empty;
* FASTQ retrieval status is valid; and
* expected RIF, INH, PZA and EMB classifications are valid `R` or `S` values.

The script then generates two representations of the final validation dataset.

### Portable Validation Manifest

`reference_dataset/tbtamr_current_study_validation_manifest.csv`

This is the human-readable operational manifest for the final AMR validation dataset.

It connects each selected isolate to:

* resistance class;
* expected RIF phenotype;
* expected INH phenotype;
* expected PZA phenotype;
* expected EMB phenotype; and
* paired R1/R2 sequencing files.

FASTQ paths are stored as repository-relative paths rather than machine-specific absolute paths.

### Galaxy Collection Manifest

`reference_dataset/tbtamr_current_study_galaxy_collection.yaml`

This represents the same 10-isolate dataset as a Galaxy `list:pair` collection for workflow execution.

Each selected isolate is represented as a paired collection containing its forward and reverse FASTQ files.

---

## Directory Structure

```text
tbtamr/
├── README.md
│   └── Documentation for the AMR reference-dataset construction workflow
│
├── validation_results_tbtamr_master.csv
│   └── Source tbtAMR dataset used to construct the AMR reference dataset
│
├── fastqs/
│   └── Paired-end FASTQ files for the 10 selected reference isolates
│
└── reference_dataset/
    ├── tbtamr_current_study_high_confidence.csv
    │   └── Complete high-confidence eligible cohort (1,805 isolates)
    │
    ├── tbtamr_current_study_reference_panel.csv
    │   └── Final selected 10-isolate reference panel with expected phenotypic DST
    │
    ├── tbtamr_current_study_reference_panel_report.txt
    │   └── Filtering, selection parameters, counts and provenance report
    │
    ├── tbtamr_current_study_fastq_manifest.csv
    │   └── FASTQ retrieval and verification record for the 10 selected isolates
    │
    ├── tbtamr_current_study_validation_manifest.csv
    │   └── Portable manifest linking expected AMR truth to R1/R2 FASTQs
    │
    └── tbtamr_current_study_galaxy_collection.yaml
        └── Galaxy list:pair representation of the 10-isolate validation dataset
```

---

## Reference Dataset vs Validation Manifest

The reference panel and validation manifest serve related but different purposes.

### Reference Panel

`reference_dataset/tbtamr_current_study_reference_panel.csv`

This is the scientific reference table for the selected isolates.

It records the expected phenotypic DST classifications and supporting information used to define and document the selected reference panel.

### Validation Manifest

`reference_dataset/tbtamr_current_study_validation_manifest.csv`

This is the operational representation of the validation dataset.

It connects the expected AMR classifications from the reference panel to the actual paired sequencing files that will be supplied to the workflow being validated.

In simplified form:

```text
Reference panel
(expected AMR classifications)
        +
Verified sequencing data
(R1/R2 FASTQs)
        │
        ▼
Final validation manifest
```

---

## Reproducibility

From the repository root, the current AMR dataset-construction workflow is:

```bash
python3 scripts/amr/01_extract_samples_from_tbtamr.py
python3 scripts/amr/02_download_fastqs_tbtamr.py
python3 scripts/amr/03_generate_manifest.py
```

The expected progression is:

```text
validation_results_tbtamr_master.csv
        │
        │ Step 01
        ▼
tbtamr_current_study_high_confidence.csv
        │
        ▼
tbtamr_current_study_reference_panel.csv
        │
        │ Step 02
        ▼
FASTQ files
+
tbtamr_current_study_fastq_manifest.csv
        │
        │ Step 03
        ▼
tbtamr_current_study_validation_manifest.csv
+
tbtamr_current_study_galaxy_collection.yaml
```

Step 03 should only generate the final validation manifests after successful reference-panel construction and FASTQ retrieval.

---

## Historical Validation Experiment

An earlier exploratory comparison between TBProfiler / MTB workflow output and an initial representation of the tbtAMR truth set is preserved separately under:

`../legacy/amr_validation_experiment/`

This experiment contains earlier work involving:

* normalisation of tbtAMR truth data;
* transformation of TBProfiler output into a comparable representation;
* application of comparison policies;
* schema and value comparisons; and
* classification and summarisation of mismatches.

The experiment is retained because these concepts may be useful when the formal workflow-validation component of the advISO framework is developed.

It is **not part of the current AMR reference-dataset construction workflow**.

See:

`../legacy/amr_validation_experiment/README.md`

for further information about the preser
