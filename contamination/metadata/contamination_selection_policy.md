# Contamination validation panel selection policy

## Source dataset

The contamination validation panel is derived from the 50× simulated
culture dataset from Heupink et al. (2021), BioProject PRJNA706121.

The master contamination validation manifest contains 100 simulated
culture samples at 50× MTB coverage with known contamination levels
ranging from 0% to 5%.

## Validation purpose

The panel is intended to assess the ability of an MTB bioinformatics
workflow to:

1. identify or report non-target contamination where supported by the
   workflow;
2. apply appropriate QC handling to contaminated samples; and
3. maintain accurate MTB genomic results in the presence of
   contamination within the evaluated range.

Low sequencing coverage is not evaluated by this panel. All samples
are therefore restricted to the 50× culture simulation dataset.

## Panel size

A 30-sample pilot validation panel is selected from the 100-sample
master dataset.

## Contamination-level representation

Five samples are selected at each contamination level:

- 0%: 5 samples
- 1%: 5 samples
- 2%: 5 samples
- 3%: 5 samples
- 4%: 5 samples
- 5%: 5 samples

Total: 30 samples.

## Clean controls

Five uncontaminated samples are selected from the available 0%
contamination simulations.

These samples provide negative controls for contamination-related
workflow behaviour.

## Contaminated samples

For each contamination level from 1% to 5%, five samples are selected.

Selection prioritises representation of different contaminant codes
within each contamination level rather than repeated selection of the
same contaminant type.

The contaminant categories represented in the master dataset include:

- Hsap: Homo sapiens
- Mabs: Mycobacterium abscessus
- Mavi: Mycobacterium avium
- Mint: Mycobacterium intracellulare
- Paer: Pseudomonas aeruginosa
- Sepi: Staphylococcus epidermidis
- Mix1: NTM mixture
- Mix2: mixed contaminants

Because contaminant codes are not represented equally at every
contamination percentage in the source dataset, a complete factorial
design is not imposed.

## Selection reproducibility

Selection is performed programmatically from the master validation
manifest.

For contaminated samples, the selection algorithm should maximise
contaminant-code diversity at each contamination percentage.

Where multiple samples satisfy the same selection condition, a
deterministic ordering based on source_sample_id is used to resolve
ties.

For uncontaminated controls, samples are selected deterministically
from the available 0% contamination samples using source_sample_id
ordering.

No samples are selected manually by run accession.

## Provenance

The selected pilot manifest retains the original truth-set fields,
ENA/SRA accessions, FASTQ locations, and checksums from the master
validation manifest.

The 100-sample master validation manifest is retained unchanged so
that the validation panel can be expanded in future analyses.
