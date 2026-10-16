# Data

This page describes the data used in the project, where to get it, how it is cleaned, and the rules that come with it.

Sections marked **TODO** are filled in after looking at the files directly (roadmap module M1.1).

## Source

**POG570**: tumour and matched normal sequencing from 570 patients with advanced or metastatic cancer, enrolled in the Personalized OncoGenomics (POG) program at BC Cancer, Vancouver. Each tumour had whole-genome and transcriptome sequencing.

- Download page: <https://www.bcgsc.ca/downloads/POG570/>
- Paper: Pleasance E, Titmuss E, Williamson L, et al. Pan-cancer analysis of advanced patient tumors reveals interactions between therapy and genomic landscapes. *Nat Cancer*. 2020;1:452–468. doi:10.1038/s43018-020-0050-6
- Read `POG570_README.txt` on the download page before using anything.

## Terms of use (summary)

From the POG570 README. The README itself is the binding version.

- Research use only. Commercial or other use needs contact with poginfo@bcgsc.ca first.
- Treat the data as controlled: keep it secure and anonymous.
- Never try to identify or contact participants.
- **Do not redistribute** the data without written permission. This is why `Data/` is kept out of Git and why reports only show summary results, not patient-level data.
- Acknowledge the POG program and its funders in any publication or presentation.

## Files used

| File | Size | What it contains | Reference / annotation |
|---|---|---|---|
| `POG570_TPM_expression.txt.gz` | ~32 MB | Gene expression for every tumour, as TPM | Ensembl 85 gene models |
| `POG570_small_mutations.txt.gz` | ~257 MB | Somatic SNVs and small indels | hg19, annotated with SnpEff (Cingolani et al., 2012) on Ensembl 69 |
| `Table_S1_Demographics.xlsx` | small | One row per patient: cancer type and clinical details | Supplementary Table 1 of the paper |
| `Table_S2_Treatment.xlsx` | small | Treatment history and survival | Optional, for later analyses |

Files in the release that are **not** used for now: copy number (`POG570_copy_number.txt`, CNV segments, APOLLOH and Ploidetect calls), raw read counts, and sequencing parameters (Table S7). Copy number could be added later as a third feature set.

### Things to watch out for

- **Two different Ensembl versions.** Expression uses Ensembl 85 and mutation annotation uses Ensembl 69. Gene names or IDs may not match between the two files, so genes should not be joined across them without checking.
- **hg19, not GRCh38.** Mutation positions are on the older genome build. This doesn't matter for gene-level features, but it does if positions are compared with newer databases.
- **TPM is relative.** TPM values are scaled within each sample, so they are good for comparing genes across samples, but they are not raw counts (Wagner et al., 2012). Raw counts are in a separate file if they are ever needed.
- **Metastatic biopsies.** A biopsy from a metastasis also contains normal cells from the organ it was taken from. Expression can therefore reflect the biopsy site as well as the primary tumour.

## File structure

**TODO**: fill in after opening each file (`zcat FILE | head`, `wc -l`, `head -1 | tr '\t' '\n'`).

| File | One row means | Number of rows | Patient ID column | Key columns |
|---|---|---|---|---|
| TPM expression | TODO | TODO | TODO | TODO |
| Small mutations | TODO | TODO | TODO | TODO (includes `qual`, `af_tumour`) |
| Table S1 | one patient | TODO | TODO | TODO |

- How patient IDs are written in each file, and whether they match directly: **TODO**
- Patients present in one file but missing from another: **TODO**

## Where the files go

```
Data/
├── Raw/          # files exactly as downloaded, never edited
└── Processed/    # cleaned outputs written by the pipeline
```

Raw files are never changed by hand. Everything in `Processed/` can be rebuilt from `Raw/` by running the code.

## Cleaning rules

Every rule is listed here with the reason for it. Counts are logged by the code at each step.

### Labels

| Rule | Reason |
|---|---|
| Drop patients with no cancer type | No label to learn from |
| Keep cancer types with at least ~15 patients (exact cut-off: TODO) | Too few examples per class gives unstable cross-validation |

### Expression

| Rule | Reason |
|---|---|
| Remove duplicate patients | One patient must count once |
| Flag any TPM < 0 | Impossible value, points to a file problem |
| Transform with log2(TPM + 1) | Shrinks the huge range of expression values so a few very active genes don't dominate. The +1 avoids log(0). |
| Keep the most variable genes (number: TODO) | Most genes barely change between tumours and only add noise |

### Mutations

| Rule | Reason |
|---|---|
| Filter low-confidence calls using `qual` and `af_tumour` (thresholds: TODO) | Removes likely sequencing errors |
| Keep protein-changing effects (missense, nonsense, frameshift, splice site) | Silent changes rarely affect how the tumour behaves |
| Patient × gene table for known driver genes, 1 = mutated, 0 = not | Simple, readable feature for the model |
| Tumour mutation burden (TMB) = number of somatic mutations per patient (or per megabase of genome) | Varies a lot between cancer types (e.g. high in melanoma and lung) |

## References

- Cingolani P, Platts A, Wang LL, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. *Fly*. 2012;6(2):80–92. doi:10.4161/fly.19695
- Pleasance E, Titmuss E, Williamson L, et al. Pan-cancer analysis of advanced patient tumors reveals interactions between therapy and genomic landscapes. *Nat Cancer*. 2020;1:452–468. doi:10.1038/s43018-020-0050-6
- Wagner GP, Kin K, Lynch VJ. Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. *Theory Biosci*. 2012;131:281–285. doi:10.1007/s12064-012-0162-3
