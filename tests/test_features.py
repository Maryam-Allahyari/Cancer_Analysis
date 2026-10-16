#The mutation matrix
import pandas as pd

from features import mutation_matrix, read_gene_list


def test_matrix_is_binary(mutations):
    mat = mutation_matrix(mutations, ["101", "102", "103"], ["TP53", "KRAS", "APC"])
    # patient 101 has two TP53 rows but should still get a 1, not a 2
    assert set(mat.values.ravel()) <= {0, 1}
    assert mat.loc["101", "mut_TP53"] == 1


def test_matrix_keeps_all_patients_and_genes(mutations):
    patients = ["101", "102", "103", "999"]
    genes = ["TP53", "BRAF"]  # BRAF is never mutated here
    mat = mutation_matrix(mutations, patients, genes)
    assert list(mat.index) == patients
    assert list(mat.columns) == ["mut_TP53", "mut_BRAF"]
    assert mat.loc["999"].sum() == 0
    assert mat["mut_BRAF"].sum() == 0


def test_gene_list_skips_comments(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("# comment\nTP53\n\nKRAS\n")
    assert read_gene_list(path) == ["TP53", "KRAS"]
