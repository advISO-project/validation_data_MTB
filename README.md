## **Aim of this Resource**
The goal of this repository is to provide a **high-quality, reproducible MTB validation dataset** that can be used to **validate genomics pipelines**. Validating a custom pipeline is a critical step in ensuring reproducibility, reliability, and eventual accreditation.  

Dataset focuses on **real-world, high-confidence samples** covering the key validation objectives:  

- **AMR genotype → phenotype mapping**  
- **Lineage assignment**  
- **Variant calling (SNPs and indels)**  
- **Species determination (MTB vs NTMs)**  
- **Detection of mixed or contaminated samples**  

---

## **Core Validation Objectives**

The datasets are designed to test pipelines across all critical areas:

**AMR validation**  
- Ensure genotype → phenotype mapping accurately predicts resistance.  
- Include major first-line and second-line drugs: rifampicin, isoniazid, ethambutol, pyrazinamide, fluoroquinolones, etc.  

**Lineage assignment**  
- Validate correct assignment to MTB lineages.  
- Detect mixed-lineage or ambiguous samples.  

**Variant calling**  
- **SNPs:** single nucleotide polymorphisms across the genome.  
- **Indels:** insertions/deletions, including frameshifts in AMR genes.  

**Mixed / contaminated samples**  
- Include real or verified mixed samples where possible.  
- Test ability to detect minor allele fractions.  
- Validate pipeline behavior on low-quality or contaminated reads.  

**Species determination**  
- Correctly classify MTB vs non-tubercular mycobacteria (NTMs).  
- Include samples with potential contamination or co-infection.  

---

## **Trusted Real-World Data Sources**

Dataset will be built from **trusted public resources**, with high-confidence metadata and verified sequencing:  

| Validation Target | Source / Dataset | Notes |
|------------------|-----------------|-------|
| **AMR genotype → phenotype** | tbtAMR (PRJNA857537) | High-quality isolates with verified culture-based DST. |
| **Lineage assignment** | Coll et al., 2018; Fastlin et al., 2023; TBDB barcode database | Lineage calls for major MTB lineages; TBDB barcodes define lineage-specific SNPs. |
| **SNP calling** | Sanger AMR loci (PRJEB15555) | Gold-standard SNPs at AMR loci; genome-wide validation with WGS isolates. |
| **Indel calling** | To be added | Will include validated or curated indels in key AMR genes. |
| **Mixed / contaminated samples** | To be added | Will include validated or simulated mixtures for pipeline robustness. |
| **Species determination** | MTB isolates + NTM genomes from ENA/SRA | Test pipeline’s ability to distinguish MTB from NTMs. |

---

## **Approach for Creating High-Quality Datasets**

1. Pull samples with **high-confidence metadata** (culture-based DST, verified lineage, known SNPs).  
2. Standardize **sample IDs** across different sources for consistent comparison.  
3. Include only samples with **raw FASTQs available** for pipeline input.  
4. Document **source, validation method, and confidence level** for each sample.  

---

## **Datasets Provided**

- `truth_set.csv` – Master table of samples with AMR, lineage, and metadata.  
- `truth_set_manifest.csv` – Table linking **local FASTQ paths** and **ENA URLs**.  
- `validation_set_manifest.yaml` – Galaxy/Planemo paired collection format for pipeline input.  

---

## **Validation Workflow**

This repository includes a reproducible mini-workflow for generating validation datasets for MTB genomics pipelines. The workflow automates selection of truth set samples, retrieval of associated FASTQ files from ENA, and creation of manifests for downstream analysis. This serves as the initial foundation towards a comprehensive validation dataset. This mini-pipeline performs the following steps:

### 1. Extract the Truth Set
`extract_samples_from_tbtamr.py`
- Reads the `truth_set.csv` exported from tbtamr.
- Selects the isolates needed for validation.
- Outputs a cleaned CSV: `truth_set.csv`.

### 2. Fetch FASTQs from ENA (HTTP)
`extract_fastqs_for_tbtamr.py`
- Retrieves paired-end FASTQ files directly from ENA using HTTP download links.  
- Automatically downloads all FASTQs listed in your input manifest (e.g., SRA/ENA accessions).  
- Saves all downloaded files into the directory:  
  `tbtamr/fastqs/`
- Ensures correct naming of paired reads (`_1.fastq.gz` and `_2.fastq.gz`).  
- Skips already-downloaded files to avoid re-fetching.


### 3. Generate Manifest Files
`generate_manifest.py`
Retrieves ENA FTP links using SRR accessions.
- Generates a table containing:
  - Local FASTQ paths
  - ENA URLs
  - Some metadata from truth set
- Produces:
  - `truth_set_manifest.csv`
  - `validation_set_manifest.yaml`  
    (Galaxy/Planemo paired collection format)

### 4. Optional End-to-End Pipeline
`run_pipeline.py`
- Executes all steps in sequence.


---

## How to Run the Workflow

### 1. Create the Environment

**Using Conda:**

```bash
conda env create -f env/environment.yml
conda activate validation_env
```

**Using pip:**

```bash
pip install -r env/requirements.txt
```

**Run the Entire Pipeline:**

```bash
python scripts/run_pipeline.py
```
