"""Quality checks and cleaning of the raw POG570 data.

Outputs (in --out-dir):
    clinical.csv           one row per kept patient, with the label (cancer_type)
    expression.parquet     log2(TPM + 1), patients x genes
    mutations_coding.csv   protein-changing mutations that passed QC
    tmb.csv                tumour mutation burden per patient
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from load import load_clinical, load_expression, load_mutations

log = logging.getLogger("clean")

# MISC is a mix of unrelated rare tumours, so it is not a real class
EXCLUDE_TYPES = ["MISC"]

# SnpEff impact levels that change the protein (missense, nonsense, frameshift, splice)
CODING_IMPACTS = ["HIGH", "MODERATE"]

# rough size of the human genome in megabases, used to express TMB per Mb
GENOME_MB = 3000


def check_duplicates(df, name):
    n_dup = df.index.duplicated().sum()
    if n_dup > 0:
        log.warning("%s: %d duplicated patient IDs, keeping the first", name, n_dup)
        df = df[~df.index.duplicated()]
    return df


def clean_clinical(clin, min_patients):
    """Keep patients with a usable label and cancer types with enough patients."""
    log.info("clinical: %d patients", len(clin))
    clin = clin.rename(columns={"analysis_cohort": "cancer_type",
                                "biopsy_cohort": "biopsy_group"})
    clin = clin.set_index("patient_id")
    clin = check_duplicates(clin, "clinical")

    clin = clin[clin["cancer_type"].notna()]
    clin = clin[~clin["cancer_type"].isin(EXCLUDE_TYPES)]
    log.info("after removing missing / MISC labels: %d", len(clin))

    counts = clin["cancer_type"].value_counts()
    keep_types = counts[counts >= min_patients].index
    clin = clin[clin["cancer_type"].isin(keep_types)]
    log.info("after keeping types with >= %d patients: %d patients, %d types",
             min_patients, len(clin), clin["cancer_type"].nunique())
    log.info("class sizes:\n%s", clin["cancer_type"].value_counts().to_string())
    return clin


def clean_expression(tpm, patients, min_mean_tpm):
    """Log-transform TPM and drop genes that are barely expressed in anyone."""
    tpm = check_duplicates(tpm, "expression")

    if (tpm < 0).any().any():
        raise ValueError("Negative TPM values found, check the input file")
    if tpm.isna().any().any():
        raise ValueError("Missing TPM values found, check the input file")

    missing = set(patients) - set(tpm.index)
    if missing:
        log.warning("%d labelled patients have no expression data", len(missing))
    tpm = tpm.loc[tpm.index.intersection(patients)]

    # this filter does not look at the labels, so it is safe to do before the train/test split
    n_before = tpm.shape[1]
    tpm = tpm.loc[:, tpm.mean() >= min_mean_tpm]
    log.info("expression: kept %d of %d genes with mean TPM >= %s",
             tpm.shape[1], n_before, min_mean_tpm)

    return np.log2(tpm + 1)


def clean_mutations(muts, patients, min_af, min_alt_reads):
    """Remove low-confidence calls, then split into TMB and coding mutations."""
    muts = muts[muts["patient_id"].isin(patients)]
    log.info("mutations for kept patients: %d", len(muts))

    muts = muts.drop_duplicates(subset=["patient_id", "chrom", "pos", "ref", "alt"])
    log.info("after removing duplicate calls: %d", len(muts))

    passed = (muts["af_tumour"] >= min_af) & (muts["tumour_alt_count"] >= min_alt_reads)
    muts = muts[passed]
    log.info("after QC (AF >= %s, alt reads >= %d): %d", min_af, min_alt_reads, len(muts))

    # TMB uses every QC-passed mutation, not only coding ones
    tmb = muts.groupby("patient_id").size().reindex(patients, fill_value=0)
    tmb = (tmb / GENOME_MB).rename("tmb_per_mb").to_frame()

    coding = muts[muts["impact"].isin(CODING_IMPACTS) & (muts["gene_id"] != "None")]
    log.info("protein-changing mutations: %d", len(coding))
    return coding, tmb


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--out-dir", default="Data/Processed")
    parser.add_argument("--min-patients", type=int, default=15)
    parser.add_argument("--min-mean-tpm", type=float, default=1.0)
    parser.add_argument("--min-af", type=float, default=0.05)
    parser.add_argument("--min-alt-reads", type=int, default=3)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    clin = clean_clinical(load_clinical(args.clinical), args.min_patients)

    expr = clean_expression(load_expression(args.expression), clin.index, args.min_mean_tpm)
    clin = clin.loc[expr.index]  # only keep patients that have expression

    coding, tmb = clean_mutations(load_mutations(args.mutations), clin.index,
                                  args.min_af, args.min_alt_reads)

    clin.to_csv(out / "clinical.csv")
    expr.to_parquet(out / "expression.parquet")
    coding.to_csv(out / "mutations_coding.csv", index=False)
    tmb.to_csv(out / "tmb.csv")
    log.info("final: %d patients, %d genes. Files written to %s", len(clin), expr.shape[1], out)


if __name__ == "__main__":
    main()
