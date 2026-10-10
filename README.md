# *Mycobacterium Tuberculosis* (MTB) Bioinformatics Workflow Validation Datasets


This repository contains curated reference datasets, validation manifests, supporting metadata, reproducible scripts, and documentation for evaluating *Mycobacterium tuberculosis* (MTB) bioinformatics workflows.

The work forms part of the advISO project, which aims to support the development, validation, and implementation of quality-assured bioinformatics workflows in clinical and public health laboratories.


## 1. Background

Clinical bioinformatics workflows require documented evidence that they are suitable for their intended use. International standards, including ISO 15189, provide a framework for evaluating analytical performance, reproducibility, traceability, and the management of changes to validated methods in laboratories.

This repository supports these requirements by developing structured, traceable validation datasets for MTB genomic analysis.

The emphasis is on establishing reliable expected results, documenting dataset selection and processing, and providing reusable resources against which workflow outputs can be evaluated.

The repository provides validation resources; it does not, by itself, establish ISO accreditation or demonstrate that a particular workflow has met all validation acceptance criteria.


## 2. Validation Objectives

The current work is organized around four validation objectives.

1. Antimicrobial resistance (AMR) prediction
2. Lineage assignment
3. Variant calling
4. Species determination
5. Detection of mixed or contaminated samples
6. End-to-end workflow reproducibility

Each objective has its own data sources, selection methodology, reference datasets, scripts, and supporting documentation.


## 3. Repository Structure

```text
validation_data_MTB/
├── contamination/
│   ├── metadata/
│   ├── reference_dataset/
│   ├── reports/
│   └── README.md
│
├── lineage/
│   ├── metadata/
│   │   └── fastlin/
│   ├── pilot/
│   └── README.md
│
├── mixed_infection/
│   ├── metadata/
│   ├── reference_dataset/
│   ├── reports/
│   └── README.md
│
├── tbtamr/
│   └── reference_dataset/
│
├── scripts/
│   ├── amr/
│   ├── contamination/
│   ├── lineage/
│   └── mixed_infection/
│
├── references/
├── docs/
├── env/
├── workflow/
├── .gitignore
└── README.md
```

The objective-specific directories contain validation resources and supporting evidence. The `scripts/` directory contains utilities for dataset preparation and processing, while `workflow/` contains workflow-related resources.


## 4. AMR Genotype-to-Phenotype Validation

**Primary source:** tbtAMR

This objective focuses on evaluating antimicrobial resistance predictions against reference resistance information derived from published MTB datasets.

Phenotypic DST results serve as the ground truth, against which resistance predictions derived from WGS data are evaluated.The repository contains curated validation results and sample manifests supporting the selection and processing of MTB sequencing datasets.

The validation approach involves:

- Identifying samples with paired WGS and phenotypic DST results.
- Curating reference datasets and generating validation manifests.
- Running the bioinformatics workflow to obtain genotypic resistance predictions.
- Comparing predicted resistance profiles against phenotypic DST results.
- Assessing concordance and relevant performance metrics.



## 5. Lineage Assignment Validation

**Primary source:** Fastlin

This objective evaluates the ability of MTB bioinformatics workflows to assign the expected lineage and sub-lineage.

Fastlin and TB-Profiler classifications were used to identify candidate samples. Rather than accepting these classifications as ground truth, selected lineage assignments were further confirmed using sequencing evidence at established lineage-defining SNP positions.

The preparation process includes:

1. Integrating Fastlin and TB-Profiler classification information.
2. Identifying and selecting suitable single-lineage candidates.
3. Preparing lineage-defining SNP barcode references and expected markers.
4. Aligning sequencing reads and examining SNP-level evidence to confirm expected lineage assignments.
5. Refining the reference panel based on SNP confirmation results.
6. Evaluating lineage-supporting observations.

Supporting metadata, candidate selection tables, pilot manifests, reports, and scripts are maintained in the `lineage/` and `scripts/lineage/` directories.

See `lineage/README.md` for the detailed methodology.


## 6. Mixed-Infection Detection Validation

**Primary source:** Sobkowiak et al. (2018)

This objective evaluates whether MTB workflows can identify samples containing more than one MTB strains.

The reference dataset includes known mixed-sample compositions and expected major and minor lineage information derived from the published study.

The dataset preparation process includes:

- Extracting reference sample information.
- Recording expected mixed-infection status.
- Documenting major and minor lineage assignments.
- Capturing expected mixture proportions.
- Generating a validation manifest.
- Recording data provenance and integrity checks.


These resources support subsequent comparisons between workflow predictions and known mixed-sample compositions.

See `mixed_infection/README.md` for further documentation.


## 7. Contamination Detection Validation

**Primary source:** Heupink et al. (2021)

This objective evaluates workflow behaviour when MTB sequencing data contain reads originating from contaminating organisms or simulated contamination scenarios.

The reference dataset is constructed from published contamination simulation data and associated sequencing metadata.

The current preparation process includes:

1. Extracting expected contamination information from the published dataset.
2. Applying documented sample selection criteria.
3. Linking selected samples to sequencing accessions.
4. Generating a validation manifest.
5. Selecting a balanced pilot dataset across contamination levels.
6. Preparing FASTQ download information for subsequent workflow evaluation.

The current reference panel contains **100 samples** at the selected coverage level, with a **30-sample pilot panel** comprising five samples at each contamination level from 0% to 5%.

The repository includes the truth set, validation manifest, pilot manifest, selection policy, metadata, and preparation reports.

See `contamination/README.md` for the detailed selection methodology and supporting evidence.


## 8. General Validation Approach

Although each objective has its own reference dataset and performance characteristics, the overall validation approach follows a common structure.

**Reference data identification**

Identify published datasets with suitable expected results and sufficient supporting metadata.

**Dataset preparation**

Extract, standardize, and curate reference information while retaining links to the original sources.

**Sample selection**

Apply objective-specific inclusion and exclusion criteria, documenting selection decisions and dataset limitations.

**Validation manifest generation**

Create structured sample manifests linking identifiers, expected results, and sequencing resources.

**Workflow evaluation**

Execute the workflow under evaluation using the selected datasets and collect the relevant outputs.

**Performance assessment**

Compare observed outputs with the defined reference results and calculate appropriate objective-specific performance measures.

**Documentation and traceability**

Maintain scripts, reports, source references, dataset versions, and integrity information needed to reproduce and review the validation process.

Not every objective has completed every stage of this process. The objective-specific documentation should be consulted for current progress and outstanding work.


## 9. Reproducibility and Data Management

The repository is designed to support transparent and reproducible dataset preparation.

Where applicable, this includes:

- Public source references and accession identifiers.
- Dataset selection criteria and supporting metadata.
- Version-controlled scripts.
- Structured reference datasets and manifests.
- Data preparation reports.
- Checksums for integrity verification.
- Documentation of assumptions and limitations.

Large sequencing files, such as raw FASTQs and alignment BAMs, are excluded from routine Git version control. They should be obtained from the referenced data sources using the appropriate manifests and preparation procedures.

Local legacy files and temporary analysis outputs are also excluded where specified in `.gitignore`.

## 10. Current Development Status

| Validation objective | Current status |
|---|---|
| AMR genotype-to-phenotype | Reference dataset and supporting validation resources prepared |
| Lineage assignment | Reference dataset and supporting validation resources prepared |
| Mixed-infection detection | Reference dataset and supporting validation resources prepared |
| Contamination detection | Reference dataset and supporting validation resources prepared |

Additional validation objectives, including speciation, SNP and indel detection, will be incorporated as the framework develops. Current objectives are also continously refined.

## 11. Using the Repository

To work with the validation resources:

1. Clone the repository.
2. Review the README and metadata for the relevant validation objective.
3. Identify the required source data and software dependencies.
4. Obtain sequencing data using the supplied accessions and manifests.
5. Run the relevant dataset preparation or workflow evaluation procedures.
6. Compare observed results with the reference data.
7. Document performance results, limitations, and any deviations from the defined methodology.


The scripts and datasets are intended to be used with their corresponding documentation. Users should verify dependencies, file paths, and required inputs before running individual scripts.

## 12. Data Sources and Attribution

The validation resources are derived from publicly available research datasets.

Key sources include:

- **tbtAMR:** MTB antimicrobial resistance validation data.
- **Fastlin:** MTB lineage assignment and classification resources.
- **Sobkowiak et al. (2018):** Mixed MTB infection reference data.
- **Heupink et al. (2021):** MTB contamination simulation datasets.

Full source information, accession references, and methodological details should be obtained from the relevant objective-specific documentation.

Users should cite the original publications and datasets when reusing these resources.

## 13. advISO Project

The advISO project develops practical resources to support the implementation of quality-assured bioinformatics workflows and the adoption of ISO-aligned validation practices.

Project information:

https://www.cardiff.ac.uk/adviso-bioinformatics-accreditation

## 14. Citation

If you use resources from this repository, please acknowledge the advISO project and cite the original publications and datasets relevant to the validation objectives used.

This repository is under active development. Dataset preparation, pilot analyses, workflow evaluation, and supporting documentation may be updated as validation activities progress.