from pathlib import Path

import pandas as pd

from annotate_cases import (EVIDENCE_COLS, build_worksheet, misclassified_cases, read_drivers,
                            split_hgvs, tidy_variants)

TEMPLATE = Path(__file__).parents[1] / "docs" / "templates" / "variant_annotation_template.csv"


#1. HGVS strings from SnpEff come as "p.../c..."; splice variants may have only the c. part.
def test_split_hgvs():
    assert split_hgvs("p.Gly12Asp/c.35G>A") == ("p.Gly12Asp", "c.35G>A")
    assert split_hgvs("c.1234+1G>A") == ("", "c.1234+1G>A")
    assert split_hgvs("None") == ("", "")


#2. Only patients whose prediction is wrong should be kept, with their clinical details.
def test_only_misclassified_kept():
    preds = pd.DataFrame({"patient_id": [1, 2, 3], "true": ["BRCA", "COLO", "LUNG"],
                          "predicted": ["BRCA", "PANC", "LUNG"]})
    clin = pd.DataFrame({"patient_id": ["1", "2", "3"], "tumour_type": ["a", "b", "c"]})
    cases = misclassified_cases(preds, clin)
    assert cases["patient_id"].tolist() == ["2"]
    assert cases["tumour_type"].tolist() == ["b"]


#3. Driver-gene variants should be flagged and listed first for each patient.
def test_drivers_first():
    coding = pd.DataFrame({
        "patient_id": ["2", "2"], "gene_id": ["TTN", "KRAS"],
        "aa_change": ["p.Ala1Val/c.2C>T", "p.Gly12Arg/c.34G>C"],
        "effect": ["missense_variant"] * 2, "impact": ["MODERATE"] * 2,
        "af_tumour": [0.5, 0.3], "tumour_alt_count": [10, 8], "transcript_id": ["t1", "t2"],
        "chrom": ["2", "12"], "pos": [1, 2], "ref": ["C", "G"], "alt": ["T", "C"],
    })
    out = tidy_variants(coding, drivers={"KRAS"})
    assert out["gene_id"].tolist() == ["KRAS", "TTN"]
    assert out["driver_gene"].tolist() == [True, False]
    assert out.loc[0, "hgvs_p"] == "p.Gly12Arg"


#4. The worksheet must have exactly the columns of the committed template, so curators can paste rows.
def test_worksheet_matches_template():
    variants = pd.DataFrame({
        "patient_id": ["2"], "gene_id": ["KRAS"], "hgvs_p": ["p.Gly12Arg"], "hgvs_c": ["c.34G>C"],
        "effect": ["missense_variant"], "af_tumour": [0.3], "driver_gene": [True],
    })
    cases = pd.DataFrame({"patient_id": ["2"], "true": ["PANC"], "predicted": ["LUNG"],
                          "biopsy_group": ["Liver"]})
    sheet = build_worksheet(variants, cases)
    template_cols = TEMPLATE.read_text().strip().split(",")
    assert list(sheet.columns) == template_cols
    assert sheet[EVIDENCE_COLS].eq("").all().all()


#5. Comment lines in the driver list are ignored.
def test_read_drivers(tmp_path):
    path = tmp_path / "drivers.txt"
    path.write_text("# comment\nTP53\n\nKRAS\n")
    assert read_drivers(path) == {"TP53", "KRAS"}
