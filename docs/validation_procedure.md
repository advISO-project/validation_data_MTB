# MTB Workflow Validation Procedure

## Aim

To validate a Galaxy workflow for *Mycobacterium tuberculosis* sequence analysis using curated validation datasets and truth sets.

## Validation process

1. Obtain validation FASTQ files using the sample manifests.
2. Upload FASTQs to your workflow as required.
3. Run the workflow.
4. Collate workflow outputs.
5. Compare observed workflow outputs against expected truth set results using the validation toolkit in the advISO github.
6. Calculate concordance and investigate discordant results.
7. Record workflow versions, database versions and validation outcomes.

## Current focus

The current validation focus is antimicrobial resistance prediction using the tbtAMR validation dataset.

## Under development

- Lineage validation
- Mixed infection validation
- Species determination validation
- IWC workflow publication
