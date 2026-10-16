"""Functions to read the raw POG570 files into pandas DataFrames."""

import pandas as pd

# Columns from the mutation file that are used later on.
# The full file has 18 columns and ~8 million rows, so only reading these saves memory.
MUTATION_COLS = [
    "patient_id", "chrom", "pos", "ref", "alt", "gene_id",
    "effect", "impact", "af_tumour", "tumour_alt_count",
]


def load_expression(path):
    """Read the TPM table and return it as patients (rows) x genes (columns).

    In the raw file, genes are rows (Ensembl IDs) and each column is a patient ID.
    """
    tpm = pd.read_csv(path, sep="\t", index_col=0)
    tpm = tpm.T  # flip so that one row = one patient
    tpm.index = tpm.index.astype(str)
    tpm.index.name = "patient_id"
    return tpm


def load_mutations(path, columns=MUTATION_COLS):
    """Read the somatic mutation table (SNVs and indels, hg19)."""
    muts = pd.read_csv(path, sep="\t", usecols=columns,
                       dtype={"patient_id": str, "chrom": str})
    return muts


def load_clinical(path):
    """Read Table S1 (one row per patient) and tidy up the column names."""
    clin = pd.read_excel(path, sheet_name="Sheet1")
    clin.columns = clin.columns.str.lower()
    clin["patient_id"] = clin["patient_id"].astype(str)

    # age has the text "unknown" for some patients
    clin["age"] = pd.to_numeric(clin["age"], errors="coerce")
    return clin
