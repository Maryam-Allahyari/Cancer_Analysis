import numpy as np
import pandas as pd
import pytest

from clean import check_duplicates, clean_clinical, clean_expression, clean_mutations


#1. To check for duplicate cases.
def test_duplicates_removed(tpm_raw):
    expr = tpm_raw.T
    doubled = pd.concat([expr, expr.iloc[[0]]])
    assert check_duplicates(doubled, "test").index.is_unique


#2. TPM (transcripts per million), which shows how active a gene is, shouldn't be negative!
def test_negative_tpm_raises(tpm_raw):
    expr = tpm_raw.T.copy()
    expr.iloc[0, 0] = -1
    with pytest.raises(ValueError):
        clean_expression(expr, expr.index, min_mean_tpm=0)


#3. To  log2(TPM+1) so that none of the genes dominate the model.
def test_log_transform():
    expr = pd.DataFrame({"g1": [0.0, 1.0, 3.0], "g2": [7.0, 15.0, 31.0]},
                        index=["a", "b", "c"])
    out = clean_expression(expr, expr.index, min_mean_tpm=0)
    assert np.allclose(out["g1"], [0, 1, 2])
    assert np.allclose(out["g2"], [3, 4, 5])


def test_low_genes_dropped():
    expr = pd.DataFrame({"low": [0.1, 0.2], "high": [10.0, 20.0]}, index=["a", "b"])
    out = clean_expression(expr, expr.index, min_mean_tpm=1)
    assert list(out.columns) == ["high"]


def test_small_classes_and_misc_removed(clinical):
    clin = clean_clinical(clinical, min_patients=2)
    assert "MISC" not in clin["cancer_type"].values
    assert set(clin["cancer_type"]) == {"BRCA", "COLO"}


def test_mutation_qc(mutations):
    coding, tmb = clean_mutations(mutations, ["101", "102", "103"], min_af=0.05, min_alt_reads=3)
    # the duplicated TP53 call counts once and patient 102's weak call is removed
    assert tmb.loc["101", "tmb_per_mb"] * 3000 == 2
    assert tmb.loc["102", "tmb_per_mb"] == 0
    # intergenic (MODIFIER) mutations are not "coding"
    assert set(coding["gene_id"]) == {"TP53", "KRAS", "APC"}
