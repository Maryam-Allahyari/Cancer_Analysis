# Shared tiny fake datasets (5 patients, 10 genes) used by all tests
import numpy as np
import pandas as pd
import pytest

PATIENTS = ["101", "102", "103", "104", "105"]
GENES = [f"ENSG0000000000{i}" for i in range(10)]


@pytest.fixture
def tpm_raw():
    """TPM table in the same layout as the raw file: genes x patients."""
    rng = np.random.default_rng(1)
    values = rng.uniform(0, 50, size=(len(GENES), len(PATIENTS))).round(2)
    return pd.DataFrame(values, index=pd.Index(GENES, name="genes"), columns=PATIENTS)


@pytest.fixture
def tpm_file(tmp_path, tpm_raw):
    path = tmp_path / "tpm.txt.gz"
    tpm_raw.to_csv(path, sep="\t")
    return path


@pytest.fixture
def clinical():
    return pd.DataFrame({
        "patient_id": PATIENTS,
        "age": [50, "unknown", 61, 45, 70],
        "analysis_cohort": ["BRCA", "BRCA", "COLO", "COLO", "MISC"],
        "biopsy_cohort": ["Liver", "Breast", "Liver", "Abdomen", "Liver"],
    })


@pytest.fixture
def mutations():
    return pd.DataFrame({
        "patient_id": ["101", "101", "101", "102", "103", "103"],
        "chrom": ["17", "17", "12", "3", "5", "5"],
        "pos": [100, 100, 200, 300, 400, 500],
        "ref": ["C"] * 6,
        "alt": ["T"] * 6,
        "gene_id": ["TP53", "TP53", "KRAS", "PIK3CA", "APC", "None"],
        "effect": ["missense_variant", "missense_variant", "missense_variant",
                   "missense_variant", "stop_gained", "intergenic_region"],
        "impact": ["MODERATE", "MODERATE", "MODERATE", "MODERATE", "HIGH", "MODIFIER"],
        "af_tumour": [0.4, 0.4, 0.3, 0.02, 0.5, 0.3],
        "tumour_alt_count": [20, 20, 15, 1, 30, 10],
    })
