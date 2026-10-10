# MTB Lineage Validation Dataset

## Purpose

This directory contains the data and intermediate outputs used to construct a reproducible *Mycobacterium tuberculosis* lineage validation dataset for the advISO project.

The lineage dataset is designed for objective evaluation of MTB lineage-calling workflows.

Published or computational lineage assignments are treated as **candidate hypotheses**, rather than being accepted directly as reference truth. Independent support for candidate lineage assignments is obtained using lineage-defining barcode markers evaluated against sequencing data using a frozen H37Rv reference genome.

The current implementation focuses on single-lineage isolates. Mixed-lineage candidates are identified and retained separately for future validation work.

---

## Dataset Construction Overview

The current lineage dataset construction workflow is:

```text
FastLin study supplementary dataset
        ↓
Integrate FastLin and TB-Profiler classifications
        ↓
Classify candidate samples
        ↓
3,039 concordant single-lineage candidates
        ↓
Summarise lineage representation
        ↓
Deterministic balancing
        ↓
1,074-sample balanced candidate pool
        ↓
20-sample pilot selection
        ↓
Retrieve sequencing reads
        ↓
Generate expected terminal lineage markers
        ↓
Prepare alignment evidence
        ↓
Extract observed terminal-marker evidence
        ↓
Apply documented confirmation policy
        ↓
20 pilot truth candidates
```

Candidate selection and independent lineage confirmation are deliberately separated.

Agreement between FastLin and TB-Profiler is used to identify suitable candidate samples, but this agreement is not treated as the independent truth assignment.

---

# Source Resources

## FastLin Candidate Dataset

The starting candidate collection is derived from supplementary material associated with the FastLin lineage-classification study.

Source file:

```text
metadata/fastlin/Supplementary_material_1.xlsx
```

The workbook contains classifications produced using FastLin and TB-Profiler.

These classifications provide the initial candidate lineage hypotheses used during dataset construction.

---

## TBDB Lineage Barcode

Lineage-defining barcode markers are obtained from a frozen TBDB barcode resource maintained separately under:

```text
../references/tbdb/
```

The reference area contains the original barcode resource, provenance information, checksums, and the processed canonical lineage barcode reference used by this workflow.

---

## H37Rv Reference Genome

H37Rv accession:

```text
NC_000962.3
```

is used as the fixed reference coordinate system.

The frozen reference genome and associated provenance are maintained under:

```text
../references/reference/
```

The reference resource is maintained separately from the lineage candidate data so that reference provenance can be independently documented and controlled.

---

# Step 01 — Integrate FastLin Study Classifications

Script:

```text
../scripts/lineage/01_merge_fastlin_tbprofiler.py
```

This step parses the relevant FastLin supplementary workbook sheets and integrates the FastLin and TB-Profiler classifications into a single candidate dataset.

The purpose of this step is **dataset integration only**.

No samples are selected as validation truth at this stage.

## Outputs

### `metadata/fastlin/candidate_dataset_master.csv`

Master integrated candidate dataset.

Contains the available FastLin and TB-Profiler information for each sample together with derived fields describing dataset presence, lineage classifications, mixture indicators, and agreement between classifications.

This file is the principal input to Step 02.

### `metadata/fastlin/candidate_dataset_master_merge_report.txt`

Human-readable report describing the dataset integration step.

Records information about the source datasets, merge results, sample representation, and integration outcome.

This report provides an audit record of how the master candidate dataset was constructed.

---

# Step 02 — Classify Candidate Samples

Script:

```text
../scripts/lineage/02_classify_candidates.py
```

This step evaluates the integrated candidate dataset and separates samples according to their suitability for the initial single-lineage validation dataset.

Candidates are classified as:

- concordant single-lineage candidates;
- mixed-lineage candidates; or
- excluded candidates.

The current dataset contains:

```text
3,039 concordant single-lineage candidates
2 mixed-lineage candidates
200 excluded candidates
```

The 3,039 single-lineage candidates represent 120 exact lineage labels across 12 broad lineage groups.

## Outputs

### `metadata/fastlin/candidate_classification_master.csv`

Complete classification audit table.

Contains all samples from the integrated candidate dataset together with their assigned candidate classification and classification reason.

This file preserves the full classification decision trail rather than containing only retained samples.

### `metadata/fastlin/single_lineage_candidates.csv`

Subset containing samples classified as suitable concordant single-lineage candidates.

These samples form the candidate population used for the initial lineage validation dataset.

This file is the principal input for subsequent representation analysis and deterministic selection.

### `metadata/fastlin/mixed_lineage_candidates.csv`

Subset containing samples identified as mixed-lineage candidates.

These samples are excluded from the current single-lineage pilot but retained as a separate candidate group for potential future mixed-lineage validation work.

### `metadata/fastlin/excluded_candidates.csv`

Subset containing samples that did not satisfy the candidate-selection requirements.

The associated classification fields preserve the reason each sample was excluded.

### `metadata/fastlin/candidate_classification_report.txt`

Human-readable summary of the candidate classification step.

Records overall classification counts and provides an audit summary of how the integrated candidate dataset was divided into single-lineage, mixed-lineage, and excluded groups.

---

# Step 03 — Build the Canonical Lineage Barcode Reference

Script:

```text
../scripts/lineage/03_build_lineage_barcode_reference.py
```

This step prepares the independent lineage-marker resource used later for candidate confirmation.

The frozen TBDB barcode is standardised against the frozen H37Rv reference genome.

The resulting resource is maintained under:

```text
../references/tbdb/
```

rather than under the lineage candidate directory because it is a reusable reference resource.

## Outputs

### `../references/tbdb/processed/lineage_barcode_reference.tsv`

Canonical processed lineage barcode reference.

Contains the standardised lineage-defining marker information used by subsequent lineage confirmation steps.

This is the barcode resource consumed by Step 08.

### `../references/tbdb/processed/lineage_barcode_reference_report.txt`

Human-readable report describing construction of the canonical barcode reference.

Records information about the source barcode, reference genome, processing checks, and resulting marker resource.

### `../references/tbdb/metadata/lineage_barcode_reference.tsv.sha256`

SHA-256 checksum for the canonical processed barcode reference.

Used to verify the identity and integrity of the processed reference.

### `../references/tbdb/metadata/lineage_barcode_reference_report.txt.sha256`

SHA-256 checksum for the barcode-reference construction report.

Used to preserve integrity and provenance of the associated report.

---

# Step 04 — Summarise Single-Lineage Candidate Representation

Script:

```text
../scripts/lineage/04_summarise_single_lineage_candidates.py
```

This step characterises lineage representation within the 3,039 concordant single-lineage candidates.

Representation is summarised at both exact-lineage and broad-lineage levels.

No samples are selected or assigned truth during this step.

## Outputs

### `metadata/fastlin/single_lineage_candidate_summary.csv`

Summary of candidate representation by exact lineage.

Provides the number of available candidate isolates associated with each exact lineage represented in the candidate collection.

This information is used to understand lineage imbalance before deterministic selection.

### `metadata/fastlin/broad_lineage_candidate_summary.csv`

Summary of candidate representation after grouping samples into broad lineage categories.

Provides a higher-level view of representation across the major lineage groups.

### `metadata/fastlin/single_lineage_candidate_summary_report.txt`

Human-readable report describing the representation analysis.

Records the number of samples, exact lineage groups, broad lineage groups, and other summary information used to assess the candidate population.

---

# Step 05 — Deterministic Candidate Selection

Script:

```text
../scripts/lineage/05_select_single_lineage_candidates.py
```

The candidate collection is highly uneven across lineage groups.

This step applies a deterministic balancing strategy to reduce over-representation of common lineages while retaining representation of less frequent lineages.

The current balancing policy produces:

```text
1,074 balanced candidate isolates
```

A reproducible 20-sample pilot is subsequently selected from this balanced candidate pool.

The selection process uses documented deterministic rules and seeds so that the same input dataset and policy can reproduce the same selection.

## Outputs

### `metadata/fastlin/balanced_single_lineage_candidates.csv`

The complete balanced candidate pool.

Contains the 1,074 isolates retained after applying the deterministic lineage-balancing policy.

This is the population from which the pilot dataset is selected.

### `metadata/fastlin/balanced_single_lineage_selection_summary.csv`

Summary of the balancing operation by lineage.

Documents candidate availability and resulting representation after application of the balancing policy.

This file provides a compact audit of how the 3,039-candidate population was reduced to the 1,074-sample balanced pool.

### `metadata/fastlin/pilot_single_lineage_candidates.csv`

The 20 samples selected for the current pilot validation dataset.

Contains the selected samples together with selection-related fields documenting their lineage classification and pilot-selection status.

This file defines the **current lineage pilot sample set**.

### `metadata/fastlin/single_lineage_selection_policy.json`

Machine-readable record of the deterministic selection policy.

Stores the selection rules, seeds, target sizes, and related parameters required to reproduce the balancing and pilot-selection procedure.

### `metadata/fastlin/single_lineage_selection_report.txt`

Human-readable audit report for the selection procedure.

Summarises the balancing process, resulting candidate pool, pilot selection, and selection parameters.

---

# Step 06 — Generate the Pilot Sequencing-Read Manifest

Script:

```text
../scripts/lineage/06_generate_pilot_fastq_manifest.py
```

This step identifies sequencing-read information for the 20 selected pilot samples.

The selected sample accessions are queried against ENA and evaluated for suitable paired sequencing-read records.

This stage establishes the link between the selected candidate samples and the sequencing data required for independent confirmation.

## Outputs

### `metadata/fastlin/pilot_single_lineage_fastq_manifest.csv`

Sequencing-read manifest for the selected pilot samples.

Contains the sequencing-read information required to retrieve the corresponding data for each selected sample.

The manifest also records the resolution/status of the sequencing-read lookup.

This file is the primary input to Step 07.

### `metadata/fastlin/pilot_single_lineage_fastq_manifest_report.txt`

Human-readable report describing construction of the sequencing-read manifest.

Summarises the number of selected samples successfully resolved to suitable sequencing reads and records relevant manifest-generation information.

---

# Step 07 — Retrieve Pilot Sequencing Reads

Script:

```text
../scripts/lineage/07_download_pilot_fastqs.py
```

This step retrieves the sequencing reads identified in the Step 06 manifest and verifies the downloaded files against the available source information.

The current pilot contains:

```text
20 samples
40 paired FASTQ files
```

## Outputs

### `pilot/fastqs/`

Directory containing the sequencing reads for the current 20-sample pilot.

Each sample is represented by paired FASTQ files.

These reads are used for the independent evidence-generation stages of the lineage workflow.

### `pilot/reports/pilot_single_lineage_fastq_download_manifest.csv`

Download audit manifest.

Records the local sequencing files associated with each pilot sample together with download and verification status information.

This manifest provides the connection between the selected samples and the local files used in downstream processing.

### `pilot/reports/pilot_single_lineage_fastq_download_report.txt`

Human-readable sequencing-read download report.

Summarises the retrieval outcome and records whether the required sequencing data were successfully obtained for the pilot samples.

---

# Step 08 — Generate Expected Terminal Lineage Markers

Script:

```text
../scripts/lineage/08_generate_expected_barcode_markers.py
```

This step uses each sample's candidate lineage assignment together with the canonical TBDB barcode reference to identify the lineage-defining marker or markers expected for that exact candidate lineage.

Only the terminal lineage assignment used for the candidate hypothesis is included in the primary expectation set.

These outputs describe **expected evidence only**.

They do not confirm that the expected marker is present in the sequencing data.

## Outputs

### `pilot/barcode_expectations/pilot_expected_terminal_barcode_markers.tsv`

Detailed sample-marker expectation table.

Contains the expected terminal lineage marker information for each pilot sample.

This table links the candidate lineage hypothesis to the independent barcode position or positions that will subsequently be evaluated in the sequencing data.

It is the primary marker input for Step 10.

### `pilot/barcode_expectations/pilot_expected_terminal_barcode_summary.csv`

Sample-level summary of expected terminal barcode markers.

Provides a compact overview of the expected marker information associated with each pilot sample.

### `pilot/barcode_expectations/pilot_expected_terminal_barcode_report.txt`

Human-readable report describing generation of the expected-marker dataset.

Records input resources, associated checksums, sample counts, marker counts, and relevant expectation-generation information.

---

# Step 09 — Prepare Alignment Evidence

Script:

```text
../scripts/lineage/09_align_pilot_reads.py
```

This step prepares sequencing evidence for the pilot samples using the frozen H37Rv reference genome.

The resulting alignment artifacts are retained under:

```text
pilot/alignment/
```

The current pilot has valid alignment outputs for all 20 selected samples.

## Outputs

### `pilot/alignment/bam/`

Contains the sample-level alignment files and associated alignment indexes used by downstream evidence extraction.

These files are derived from the current pilot sequencing reads.

### `pilot/alignment/logs/`

Contains execution logs generated during reference-index preparation and sample processing.

These logs provide technical traceability for the alignment stage.

### `pilot/alignment/qc/`

Contains sample-level alignment quality-control summaries.

These files provide basic evidence that the generated alignment outputs completed successfully and can be assessed before downstream use.

### `pilot/alignment/reference/`

Contains reference-index artifacts generated for the alignment run.

These are derived runtime artifacts associated with the frozen H37Rv reference rather than an independent reference source.

The authoritative frozen reference and its provenance remain under:

```text
../references/reference/
```

### `pilot/alignment/reports/pilot_alignment_manifest.csv`

Sample-level alignment manifest.

Links each pilot sample to the sequencing inputs and resulting alignment outputs and records the processing status required for downstream evidence extraction.

### `pilot/alignment/reports/pilot_alignment_report.txt`

Human-readable alignment-stage report.

Summarises the alignment run, input resources, reference information, software information, sample processing outcomes, and relevant provenance.

---

# Step 10 — Extract Terminal-Marker Observations

Script:

```text
../scripts/lineage/10_extract_terminal_marker_observations.py
```

This step examines sequencing evidence at the expected terminal lineage-marker positions for each pilot sample.

The observations are classified descriptively according to the evidence present at each expected marker.

Possible observation categories include:

- expected allele observed;
- reference allele observed;
- alternative allele observed;
- mixed alleles;
- ambiguous allele observation;
- low depth; and
- no coverage.

This stage records observed evidence only.

The observation classification is kept separate from the final candidate-confirmation decision.

## Outputs

### `pilot/barcode_observations/pilot_terminal_marker_observations.tsv`

Detailed marker-level evidence table.

Contains one record for each evaluated sample-marker combination.

The table records the expected marker information together with the observed evidence and associated quality/depth information.

This is the principal evidence table consumed by Step 11.

### `pilot/barcode_observations/pilot_terminal_marker_observation_summary.csv`

Sample-level summary of the terminal-marker observations.

Summarises the number of markers falling into the different observation categories for each pilot sample.

This provides a compact view of the evidence available for each candidate.

### `pilot/barcode_observations/pilot_terminal_marker_observation_report.txt`

Human-readable observation-stage report.

Records the input resources, checksums, evidence-extraction settings, thresholds, sample counts, marker counts, and overall observation results.

---

# Step 11 — Apply the Terminal-Marker Confirmation Policy

Script:

```text
../scripts/lineage/11_confirm_terminal_lineage_candidates.py
```

This step applies the documented confirmation policy to the independent terminal-marker observations.

The policy distinguishes between supported candidates, candidates supported but affected by non-callable low-depth evidence, candidates requiring review, and candidates with conflicting evidence.

The current pilot outcome is:

```text
18 confirmed
2 confirmed_low_depth
0 review
0 rejected
```

These samples remain **truth candidates**.

They should not be described as the final frozen lineage truth set until the formal review and freezing step has been completed.

## Outputs

### `pilot/truth_candidates/pilot_terminal_truth_candidates.tsv`

Sample-level truth-candidate decision table.

Contains the confirmation outcome for each pilot sample together with the evidence summary and decision information supporting that classification.

This is the principal output of the current pilot confirmation workflow.

### `pilot/truth_candidates/pilot_terminal_truth_summary.csv`

Summary of truth-candidate confirmation outcomes.

Provides a compact overview of the resulting confirmation categories across the pilot dataset.

### `pilot/truth_candidates/pilot_terminal_truth_report.txt`

Human-readable confirmation report.

Documents the confirmation run, input evidence, outcome counts, and relevant provenance.

The report explicitly records that these samples remain truth candidates until formally reviewed and frozen.

### `pilot/truth_candidates/terminal_confirmation_policy.json`

Machine-readable representation of the confirmation policy.

Records the decision rules used to classify candidates together with relevant input identities and policy information.

This file provides a reproducible specification of how marker-level observations are translated into candidate-level confirmation outcomes.

---

# Current Directory Structure

```text
lineage/
├── README.md
│
├── metadata/
│   └── fastlin/
│       ├── Supplementary_material_1.xlsx
│       │
│       ├── candidate_dataset_master.csv
│       ├── candidate_dataset_master_merge_report.txt
│       │
│       ├── candidate_classification_master.csv
│       ├── candidate_classification_report.txt
│       ├── single_lineage_candidates.csv
│       ├── mixed_lineage_candidates.csv
│       ├── excluded_candidates.csv
│       │
│       ├── single_lineage_candidate_summary.csv
│       ├── broad_lineage_candidate_summary.csv
│       ├── single_lineage_candidate_summary_report.txt
│       │
│       ├── balanced_single_lineage_candidates.csv
│       ├── balanced_single_lineage_selection_summary.csv
│       ├── pilot_single_lineage_candidates.csv
│       ├── single_lineage_selection_policy.json
│       ├── single_lineage_selection_report.txt
│       │
│       ├── pilot_single_lineage_fastq_manifest.csv
│       └── pilot_single_lineage_fastq_manifest_report.txt
│
└── pilot/
    ├── fastqs/
    │
    ├── reports/
    │   ├── pilot_single_lineage_fastq_download_manifest.csv
    │   └── pilot_single_lineage_fastq_download_report.txt
    │
    ├── barcode_expectations/
    │   ├── pilot_expected_terminal_barcode_markers.tsv
    │   ├── pilot_expected_terminal_barcode_summary.csv
    │   └── pilot_expected_terminal_barcode_report.txt
    │
    ├── alignment/
    │   ├── bam/
    │   ├── logs/
    │   ├── qc/
    │   ├── reference/
    │   └── reports/
    │       ├── pilot_alignment_manifest.csv
    │       └── pilot_alignment_report.txt
    │
    ├── barcode_observations/
    │   ├── pilot_terminal_marker_observations.tsv
    │   ├── pilot_terminal_marker_observation_summary.csv
    │   └── pilot_terminal_marker_observation_report.txt
    │
    └── truth_candidates/
        ├── pilot_terminal_truth_candidates.tsv
        ├── pilot_terminal_truth_summary.csv
        ├── pilot_terminal_truth_report.txt
        └── terminal_confirmation_policy.json
```

Reference resources are maintained separately under:

```text
references/
├── reference/
└── tbdb/
```

Workflow scripts are maintained under:

```text
scripts/lineage/
```

---

# Reproducibility and Provenance

The lineage dataset construction workflow maintains separation and traceability between:

```text
source candidate classifications
        ↓
integrated candidate dataset
        ↓
candidate classification
        ↓
candidate representation
        ↓
deterministic selection
        ↓
pilot sequencing data
        ↓
expected independent barcode evidence
        ↓
observed sequencing evidence
        ↓
documented confirmation policy
        ↓
truth candidates
```

Where applicable, reports and policy files preserve information such as:

- source and input files;
- input checksums;
- reference-resource identities;
- deterministic selection parameters;
- random seeds;
- analysis thresholds;
- software information;
- sample-level processing status; and
- decision rules.

This separation allows the provenance of a final validation sample to be traced back through the dataset-construction process.

---

# Current Status

The current lineage workflow has completed:

- construction of the FastLin-derived candidate population;
- classification of 3,039 concordant single-lineage candidates;
- deterministic construction of a 1,074-sample balanced candidate pool;
- selection of a 20-sample pilot;
- sequencing-read retrieval for all 20 pilot samples;
- generation of expected terminal lineage-marker evidence;
- preparation of sequencing evidence for all 20 samples;
- terminal-marker observation; and
- application of the documented confirmation policy.

The current pilot contains:

```text
18 confirmed
2 confirmed_low_depth
```

The resulting files under:

```text
pilot/truth_candidates/
```

represent the **current lineage truth-candidate dataset**.

They are not yet the formally frozen final lineage validation truth set.

The next methodological stage is formal review of the truth candidates and documentation of the final freezing decision before the reference set is used for objective lineage-workflow performance evaluation.
