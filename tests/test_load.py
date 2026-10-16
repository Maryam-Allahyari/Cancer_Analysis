#Loaders return the expected columns and types
import pandas as pd

from load import load_clinical, load_expression


def test_expression_is_patients_by_genes(tpm_file, tpm_raw):
    expr = load_expression(tpm_file)
    assert expr.shape == (tpm_raw.shape[1], tpm_raw.shape[0])
    assert list(expr.index) == list(tpm_raw.columns)
    assert expr.index.name == "patient_id"


def test_expression_ids_are_strings(tpm_file):
    expr = load_expression(tpm_file)
    assert all(isinstance(i, str) for i in expr.index)


def test_clinical_columns_and_age(tmp_path, clinical):
    path = tmp_path / "s1.xlsx"
    clinical.rename(columns=str.upper).to_excel(path, sheet_name="Sheet1", index=False)
    clin = load_clinical(path)
    assert "analysis_cohort" in clin.columns  # lower-cased
    assert all(isinstance(i, str) for i in clin["patient_id"])
    assert pd.isna(clin.loc[1, "age"])  # "unknown" becomes missing
