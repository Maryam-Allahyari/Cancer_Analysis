# Variant annotation of misclassified cases: guide

This add-on asks a lab-style question about the model's mistakes: **for patients the classifier got wrong, do their somatic variants support the true cancer type, the predicted one, or neither?**

> Research exercise only. POG570 data is research-use and may not be redistributed. Patient-level worksheets and case reviews stay in `results/annotation/`, which is git-ignored.

## 1. Generate the worksheets

```bash
python src/annotate_cases.py
```

This writes three files to `results/annotation/`:

| File | Contents |
|---|---|
| `cases.csv` | One row per misclassified test patient: label, prediction, histology, biopsy site, tumour content, variant counts |
| `variants_all.csv` | Every QC-passed protein-changing variant for these patients, driver genes first |
| `annotation_worksheet.csv` | Driver-gene variants with empty evidence columns to fill in |

Variant QC matches `clean.py` (VAF ≥ 0.05, ≥ 3 alt reads, HIGH/MODERATE impact). Driver genes come from `resources/driver_genes.txt`.

### Optional: look up the facts automatically

```bash
python src/gather_evidence.py
```

This queries UniProt, CIViC, cBioPortal and ClinVar (open APIs, about 20 seconds; only gene + protein change are sent) and writes:

| File | Contents |
|---|---|
| `annotation_worksheet_evidence.csv` | The worksheet with the factual columns filled in: `civic_evidence`, `clinvar_somatic`, `cosmic_or_cbioportal` (TCGA frequency in the true and predicted type), `typical_cancer_types`, `sources` (including a UniProt check that the reference amino acid matches) |
| `evidence_frequencies.csv` | Gene and exact-variant mutation frequency in each of the 32 TCGA PanCancer Atlas cancer types |

It deliberately leaves `variant_class`, `amp_tier`, `supports_true_type`, `supports_predicted_type` and `curator_notes` empty. Those need judgment. Things it can't know:

- **Gene-level and guideline evidence.** CIViC is searched for the exact variant only, so rules recorded at gene level are missed, e.g. "KRAS mutation predicts no benefit from anti-EGFR therapy in colorectal cancer" (Tier I).
- **Germline vs somatic context.** ClinVar may list a variant as germline pathogenic even when the tumour data show it is somatic.
- **Tumours without a TCGA cohort** (e.g., mesothelioma under LUNG, anal squamous carcinoma under COLO). The comparison uses the label's cohort; check the detailed histology yourself.
- **COSMIC**, which needs a login.

## 2. Things to know about the data

- **Genome build hg19, transcripts from Ensembl v69.** Look variants up by gene and protein change (e.g., *KRAS p.Gly12Arg*) rather than coordinates, or convert with [VariantValidator](https://variantvalidator.org/) or LiftOver first.
- **Old gene names:** Ensembl 69 uses MLL3 (now **KMT2C**) and MLL2 (now **KMT2D**). Search databases with the current HGNC symbol.
- **Somatic only.** Germline variants, fusions and structural variants are not in the public release. Copy number is published separately (`POG570_copy_number.txt`) but is not downloaded yet.
- **Some patients have no driver-gene variants** in the 63-gene list. For those, review `variants_all.csv` (HIGH-impact variants first) and consider genes typical for that tumour type that the list lacks.

## 3. Evidence sources (in this order)

1. **[CIViC](https://civicdb.org/)**: clinical evidence items and variant-level summaries (open)
2. **[cBioPortal](https://www.cbioportal.org/)**: how often the variant occurs, and in which cancer types (use the TCGA PanCancer Atlas studies)
3. **[COSMIC](https://cancer.sanger.ac.uk/cosmic)**: hotspot frequency and tissue distribution (free academic registration)
4. **[OncoKB](https://www.oncokb.org/)**: oncogenicity and therapeutic levels (public website)
5. **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)**: somatic / oncogenicity classifications where present
6. **PubMed**: tumour-type-specific literature; record PMIDs

## 4. Filling in the worksheet

| Column | How to fill it |
|---|---|
| `variant_class` | hotspot / known oncogenic / likely oncogenic / VUS / likely passenger |
| `civic_evidence` | CIViC variant ID or evidence IDs, or "none" |
| `clinvar_somatic` | ClinVar oncogenicity or somatic classification, or "none" |
| `cosmic_or_cbioportal` | e.g., "seen in 3% of PAAD, <0.5% of LUAD (cBioPortal PanCancer)" |
| `amp_tier` | AMP/ASCO/CAP tier: I (strong clinical significance), II (potential), III (unknown), IV (benign/likely benign) |
| `typical_cancer_types` | Where the variant or gene is characteristically found |
| `supports_true_type` / `supports_predicted_type` | yes / no / neutral, judged from the tumour-type distribution |
| `sources` | PMIDs and database links |
| `curator_notes` | Anything that affects interpretation (low VAF, subclonal, tumour content) |

**Judging support for a cancer type:**
- **yes**: the variant is characteristic of that type (clearly enriched compared with others)
- **neutral**: common across many types (e.g., most TP53 missense)
- **no**: rare in that type and common in another

## 5. Per-patient case review

For each patient, copy `docs/templates/case_review_template.md` to `results/annotation/case_<ID>.md` and fill it in. Then check it against `docs/templates/report_checklist.md`.

## 6. Summary for the README

Once all cases are reviewed, summarize **only aggregate findings** (no patient IDs or patient-level variant lists), for example:

- how many errors traced to label groupings, similar biology, rare tumour types or low tumour content
- how many cases had variants that pointed to the true origin
- what additional data (germline, fusions, copy number) would have helped
