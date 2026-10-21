# Case review: patient `<ID>`

> Research exercise on POG570 data. Not a clinical report. Keep filled-in copies in `results/annotation/` (git-ignored); do not publish patient-level details.

## 1. Case summary

| Field | Value |
|---|---|
| Recorded cancer type (label) | |
| Detailed tumour type / histology | |
| Biopsy site | |
| Primary site | |
| Tumour content (%) | |
| Model prediction | |
| Protein-changing variants (QC-passed) | |
| Driver-gene variants | |

## 2. Is the label itself right?

Does the recorded cancer-type group fit the detailed histology? (e.g., a mesothelioma filed under LUNG, or an anal squamous carcinoma filed under COLO.)

- [ ] Label fits the histology
- [ ] Label is a broad grouping that hides a distinct tumour type → explain:

## 3. Key variants

| Gene | Variant (HGVS p.) | VAF | Variant class | AMP tier | Typical cancer types | Supports true / predicted / neutral |
|---|---|---|---|---|---|---|
| | | | | | | |

Evidence for each variant goes in the annotation worksheet (`sources` column).

## 4. Genotype-phenotype interpretation

- Do the variants support the **true** origin, the **predicted** origin, both, or neither?
- Are there variants characteristic of one type? (e.g., KRAS G12R in pancreatic cancer)
- Copy number or fusions that would be expected but can't be assessed with this data:
- Could a **germline** cause be relevant? (e.g., SDHx in paraganglioma) Note that POG570 public data is somatic only.

## 5. Why the model was likely wrong

- [ ] Label grouping (histology differs from the group)
- [ ] Biologically similar tumour types (e.g., squamous carcinomas at different sites)
- [ ] Rare tumour type with few training examples
- [ ] Biopsy-site tissue contamination (e.g., normal liver)
- [ ] Low tumour content
- [ ] Other:

## 6. Conclusion

One or two sentences: what a molecular report would add here beyond the expression-based prediction.

## 7. Report completeness check

See [report_checklist.md](report_checklist.md).

Curator: ______  Date: ______  Second review: ______
