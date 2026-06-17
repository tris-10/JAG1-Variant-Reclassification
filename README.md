# JAG1-Variant-Reclassification
## Likelihood‑based calibration improves the clinical utility of JAG1 functional data for variant classification

*Authors:* Tristan J. Hayeck, Christopher J. Sottolano, Justin J. Blair, Markos N. Xenakis, Nancy B. Spinner, and Melissa A. Gilbert

*Abstract:* Mutliplexed Assays of Variant Effects (MAVEs) represent a powerful approach to provide functional information for disease genes at scale. To harness the full utility of these systems, assay readout must be translated from laboratory results into a language that is accommodating to the clinical genomics community. We previously performed a MAVE to characterize missense variants in *JAG1*, the primary cause of the autosomal dominant, multisystemic disease, Alagille syndrome. Using this existing dataset, we derived log likelihood ratios of pathogenicity (LLRp) for each variant score, allowing for direct translation of variant data into recognized evidence weights utilized by the American College of Medical Genetics and Genomics (ACMG) and the Association for Molecular Pathology (AMP) during clinical variant classification. Calibration resulted in improved separation of known benign and pathogenic variants and increased the classification rate of abnormal missense variants from 486 to 610, providing clear binning for strong (n=1), moderate (n=340), and supporting (n=269) evidence toward pathogenicity. Retrospective application of this evidence to a cohort of 29 individuals with a *JAG1* VUS identified from clinical diagnostic sequencing yielded nine (31%) variants with abnormal data meeting supporting (n=3) and moderate (n=6) weight for pathogenicity, of which six (21%) were upgraded to likely pathogenic or pathogenic. Calibration of MAVE data to comply with the ACMG/AMP variant classification framework improves the diagnostic yield for *JAG1* variants for Alagille syndrome. Moreover, our results support the broader application of these models to additional MAVEs, suggesting that they are likely to strongly impact variant classification across disease genes. 

-----------------------------------------------------

# mave_reclassification.R

An R script for converting raw Multiplexed Assay of Variant Effect (MAVE) scores into log likelihood ratios of pathogenicity (LLRp) and mapping them to ACMG/AMP evidence categories using the [jweile/maveLLR](https://github.com/jweile/maveLLR) framework, based on fixed thresholds found in [van Loggerenberg et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37729906/).

---

## Overview

This script takes a TSV file containing MAVE functional scores for variant reference sets (pathogenic and benign controls) and a library of variants of interest. It:

1. Estimates class-conditional probability densities for pathogenic and benign reference scores via kernel density estimation (KDE)
2. Calculates an LLRp for each library variant
3. Maps each LLRp to an ACMG/AMP evidence tier based on Tavtigian et al. (2018) thresholds
4. Outputs an annotated TSV with `llrp_score` and `llrp_functional_consequence` columns appended

---

## Dependencies

The following R (version 4.5.3) packages are required and will be auto-installed if missing:

- [`devtools`](https://cran.r-project.org/package=devtools)
- [`optparse`](https://cran.r-project.org/package=optparse)
- [`dplyr`](https://cran.r-project.org/package=dplyr)
- [`maveLLR`](https://github.com/jweile/maveLLR) *(installed from GitHub)*

> **Windows users:** Uncomment the `options(download.file.method = "wininet")` line in the `maveLLR` installation block if you encounter download issues.

---

## Usage

```bash
Rscript mave_reclassification.R \
  --input  <path/to/input.tsv> \
  --output <path/to/output.tsv> \
  --scoreLabel <score_column_name> \
  --refLabel <ref_column_name> \
  [--bw <bandwidth>]
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `-i` / `--input` | ✅ | Input path for TSV file containing positive reference, negative reference, and library MAVE scores |
| `-o` / `--output` | ✅ | Output path for TSV file containing converted LLRp MAVE scores. |
| `-s` / `--scoreLabel` | ✅ | Column name containing the MAVE scores to use |
| `-r` / `--refLabel` | ✅ | Column name containing reference set categories to use (benign: -1, library: 0, pathogenic: 1)|
| `-b` / `--bw` | ❌ | KDE bandwidth (default: `0.1`) |

### Example

Using input columns from [MAVE calibration for JAG1 and Alagille syndrome — PubMed 39043182](https://pubmed.ncbi.nlm.nih.gov/39043182/). Data can be found at [MAVEdb](https://www.mavedb.org/score-sets/urn:mavedb:00001198-a-1).

```bash
Rscript mave_reclassification.R \
  --input mave_scores.tsv \
  --output mave_reclassified.tsv \
  --scoreLabel meanAcrossReps \
  --refLabel Category \
  --bw 0.1
```

---

## Input File Requirements

The input must be a **tab-delimited (TSV)** file with a header row. The following columns are required:

| Column | Description |
|--------|-------------|
| `<scoreLabel>` | The MAVE functional score column specified via `--scoreLabel` (e.g., `meanAcrossReps`) |
| `<refLabel>` | The Variant class label specified via `--refLabel` (e.g., `Category`). Used to separate variants into reference sets for KDE fitting and the library set for LLRp transformation. Ensure all three category values are present in the input: `1` = pathogenic reference, `-1` = benign reference, `0` = library variant. |

All other columns are preserved as-is in the output.

> The `Category` column is used to separate variants into reference sets for KDE fitting and the library set for LLRp transformation. Ensure all three category values (`1`, `-1`, `0`) are present in the input.

### Example Input

| variant | meanAcrossReps | Category |
|---|---|---|
| var_id_1 | 0.752667413 | 0 |
| var_id_39 | 0.753941193 | -1 |
| var_id_2 | 0.656628397 | 0 |
| var_id_41 | 0.746960144 | -1 |
| var_id_3 | 0.448961487 | 0 |

*(`0` = library variant to be classified, `-1` = benign reference, `1` = pathogenic reference — not shown above, but required somewhere in the file)*

---

## Output

The output is the input TSV with two additional columns appended:

| Column | Description |
|--------|-------------|
| `llrp_score` | The computed log likelihood ratio of pathogenicity |
| `llrp_functional_consequence` | ACMG/AMP evidence tier mapped from the LLRp score |

### Evidence Tier Thresholds (Tavtigian et al., 2018)

| LLRp Range | Label | Interpretation |
|------------|-------|----------------|
| > 2.54 | `PVSt` | Very strong pathogenic |
| 1.27 – 2.54 | `PSt` | Strong pathogenic |
| 0.63 – 1.27 | `PM` | Moderate pathogenic |
| 0.31 – 0.63 | `PSu` | Supporting pathogenic |
| -0.31 – 0.31 | `I` | Indeterminate |
| -1.27 – -0.31 | `LNorm` | Supporting benign |
| < -1.27 | `Norm` | Strong benign |

A summary density/LLR plot is also generated to the active R graphics device during execution.

### Example Output

| variant | meanAcrossReps | Category | llrp_score | llrp_functional_consequence |
|---|---|---|---|---|
| var_id_1 | 0.752667413 | 0 | -0.879932254 | LNorm |
| var_id_2 | 0.656628397 | 0 | -0.658756275 | LNorm |
| var_id_3 | 0.448961487 | 0 | -0.018621713 | I |
| var_id_4 | 0.871723154 | 0 | -1.276291047 | Norm |
| var_id_5 | 0.525012535 | 0 | -0.293072139 | I |

A summary density/LLR plot is also generated to the active R graphics device during execution.

---

## Citation

If you use this script, please cite the associated publication:

> Link coming soon

And the underlying `maveLLR` method:

> [van Loggerenberg, et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37729906/): Systematically testing human HMBS missense variants to reveal mechanism and pathogenic variation and [jweile/maveLLR](https://github.com/jweile/maveLLR): Calculate pathogenicity log likelihood ratios (LLRs) for MAVE datasets.
