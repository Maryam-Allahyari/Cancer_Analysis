"""Build variant annotation worksheets for patients the model misclassified.

For each misclassified test patient this pulls the protein-changing mutations from the raw
POG570 file (same QC as clean.py) and prepares a worksheet to fill in by hand with
evidence from CIViC, ClinVar, COSMIC/cBioPortal and the literature.

Outputs (in --out-dir, which is git-ignored because it holds POG patient data):
    cases.csv                one row per misclassified patient (true vs predicted type)
    variants_all.csv         every QC-passed protein-changing mutation for these patients
    annotation_worksheet.csv driver-gene variants with empty evidence columns to fill in
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

from clean import clean_mutations
from load import MUTATION_COLS, load_mutations

log = logging.getLogger("annotate_cases")

# the cleaned mutation table drops these, but annotation needs the HGVS change
EXTRA_COLS = ["transcript_id", "aa_change"]

# columns the curator fills in; kept identical to docs/templates/variant_annotation_template.csv
EVIDENCE_COLS = [
    "variant_class",          # hotspot / known oncogenic / likely oncogenic / VUS / likely passenger
    "civic_evidence",         # CIViC variant or evidence IDs
    "clinvar_somatic",        # ClinVar oncogenicity / somatic classification, if any
    "cosmic_or_cbioportal",   # how often seen, and in which cancer types
    "amp_tier",               # AMP/ASCO/CAP tier I-IV
    "typical_cancer_types",   # where this variant is usually found
    "supports_true_type",     # yes / no / neutral
    "supports_predicted_type",
    "sources",                # PMIDs, database links
    "curator_notes",
]


def read_drivers(path):
    """Read the driver gene list (one gene per line, # starts a comment)."""
    lines = Path(path).read_text().splitlines()
    return {g.strip() for g in lines if g.strip() and not g.startswith("#")}


def misclassified_cases(predictions, clinical):
    """Patients where the predicted cancer type differs from the true one, with clinical details."""
    wrong = predictions[predictions["true"] != predictions["predicted"]].copy()
    wrong["patient_id"] = wrong["patient_id"].astype(str)
    clin_cols = ["patient_id", "tumour_type", "histological_type", "biopsy_site",
                 "primary_site", "tumour_content"]
    clin = clinical[[c for c in clin_cols if c in clinical.columns]].copy()
    clin["patient_id"] = clin["patient_id"].astype(str)
    return wrong.merge(clin, on="patient_id", how="left")


def split_hgvs(aa_change):
    """'p.Gly12Asp/c.35G>A' -> ('p.Gly12Asp', 'c.35G>A'). Splice variants have only a c. part or 'None'."""
    hgvs_p, hgvs_c = "", ""
    for part in str(aa_change).split("/"):
        if part.startswith("p."):
            hgvs_p = part
        elif part.startswith("c."):
            hgvs_c = part
    return hgvs_p, hgvs_c


def tidy_variants(coding, drivers):
    """Add HGVS columns and a driver flag, and sort so drivers and HIGH impact come first."""
    coding = coding.copy()
    coding[["hgvs_p", "hgvs_c"]] = coding["aa_change"].apply(lambda a: pd.Series(split_hgvs(a)))
    coding["driver_gene"] = coding["gene_id"].isin(drivers)
    coding = coding.sort_values(
        ["patient_id", "driver_gene", "impact", "af_tumour"],
        ascending=[True, False, True, False],  # HIGH sorts before MODERATE
    )
    cols = ["patient_id", "gene_id", "hgvs_p", "hgvs_c", "effect", "impact", "driver_gene",
            "af_tumour", "tumour_alt_count", "transcript_id", "chrom", "pos", "ref", "alt"]
    return coding[cols].reset_index(drop=True)


def build_worksheet(variants, cases):
    """Driver-gene variants joined to the case details, with empty evidence columns."""
    case_cols = ["patient_id", "true", "predicted", "biopsy_group"]
    sheet = variants[variants["driver_gene"]].merge(cases[case_cols], on="patient_id", how="left")
    sheet = sheet.rename(columns={"true": "true_type", "predicted": "predicted_type"})
    for col in EVIDENCE_COLS:
        sheet[col] = ""
    first = ["patient_id", "true_type", "predicted_type", "biopsy_group", "gene_id", "hgvs_p",
             "hgvs_c", "effect", "af_tumour"]
    return sheet[first + EVIDENCE_COLS]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mutations", default="Data/Raw/POG570_small_mutations.txt.gz")
    parser.add_argument("--predictions", default="results/models/test_predictions.csv")
    parser.add_argument("--clinical", default="Data/Processed/clinical.csv")
    parser.add_argument("--drivers", default="resources/driver_genes.txt")
    parser.add_argument("--out-dir", default="results/annotation")
    parser.add_argument("--min-af", type=float, default=0.05)
    parser.add_argument("--min-alt-reads", type=int, default=3)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    cases = misclassified_cases(pd.read_csv(args.predictions), pd.read_csv(args.clinical))
    log.info("misclassified patients: %d", len(cases))

    muts = load_mutations(args.mutations, columns=MUTATION_COLS + EXTRA_COLS)
    coding, _ = clean_mutations(muts, cases["patient_id"].tolist(), args.min_af, args.min_alt_reads)
    variants = tidy_variants(coding, read_drivers(args.drivers))

    counts = variants.groupby("patient_id").agg(
        n_protein_changing=("gene_id", "size"), n_driver_variants=("driver_gene", "sum"))
    cases = cases.merge(counts, left_on="patient_id", right_index=True, how="left")
    worksheet = build_worksheet(variants, cases)

    cases.to_csv(out / "cases.csv", index=False)
    variants.to_csv(out / "variants_all.csv", index=False)
    worksheet.to_csv(out / "annotation_worksheet.csv", index=False)
    log.info("cases:\n%s", cases[["patient_id", "true", "predicted", "n_protein_changing",
                                  "n_driver_variants"]].to_string(index=False))
    log.info("worksheet rows (driver-gene variants): %d -> %s", len(worksheet), out)


if __name__ == "__main__":
    main()
