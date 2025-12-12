cat << 'EOF' > README.md
# Validation Workflow

Automate sample selection, FASTQ extraction, ENA metadata lookup, and manifest generation for for MTB validation data.
The workflow uses simple Python scripts and a clean directory structure, making it easy to run and reproduce.

---

## 🌟 Workflow Overview

This mini-pipeline performs the following steps:

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

## 🚀 How to Run

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
