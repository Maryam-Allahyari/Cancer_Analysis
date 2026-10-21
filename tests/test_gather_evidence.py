import pandas as pd

from gather_evidence import (current_symbol, frequencies, short_protein, summarise_civic,
                             summarise_clinvar, top_types_text, type_frequency_text)


#1. Database lookups use the short form (G12D); frameshifts and splice variants are skipped.
def test_short_protein():
    assert short_protein("p.Gly12Asp") == ("G", 12, "D", "G12D")
    assert short_protein("p.Gln153*") == ("Q", 153, "*", "Q153*")
    assert short_protein("p.Gly1122fs") is None
    assert short_protein("") is None


#2. Ensembl 69 names in the POG file must be translated before searching databases.
def test_renamed_genes():
    assert current_symbol("MLL3") == "KMT2C"
    assert current_symbol("KRAS") == "KRAS"


def _fake_tcga():
    """Two cohorts: paad with 4 samples, coadread with 5."""
    muts = pd.DataFrame({
        "study": ["paad", "paad", "paad", "coadread", "coadread"],
        "sample": ["p1", "p2", "p2", "c1", "c2"],
        "entrez": [3845] * 5,
        "protein_change": ["G12D", "G12R", "G13D", "G12D", "G12D"],
        "type": ["Missense_Mutation"] * 4 + ["Silent"],
    })
    totals = pd.Series({"paad": 4, "coadread": 5})
    return muts, totals


#3. Frequencies count samples (not mutations) and ignore silent changes.
def test_frequencies():
    muts, totals = _fake_tcga()
    freq = frequencies(muts, totals, entrez_id=3845, short="G12D")
    assert freq.loc["paad", "gene_pct"] == 50.0      # p1 and p2 out of 4
    assert freq.loc["paad", "variant_pct"] == 25.0   # only p1 has G12D
    assert freq.loc["coadread", "variant_pct"] == 20.0  # c2's G12D is silent, so excluded
    assert freq["variant_samples"].sum() == 2


#4. POG labels are mapped to TCGA cohorts; labels without a cohort say so.
def test_type_frequency_text():
    muts, totals = _fake_tcga()
    freq = frequencies(muts, totals, entrez_id=3845, short="G12D")
    assert type_frequency_text(freq, "PANC") == "PANC (paad gene 50.0%, variant 25.0%)"
    assert "no matching TCGA cohort" in type_frequency_text(freq, "MESO")
    assert top_types_text(freq, "variant_pct").startswith("paad 25.0%")


#5. Summaries handle empty results without failing.
def test_summaries():
    assert "no accepted CIViC evidence" in summarise_civic([])
    assert summarise_clinvar([]) == "not in ClinVar"
    item = {"id": 1, "evidenceType": "PROGNOSTIC", "evidenceLevel": "B",
            "significance": "POOR_OUTCOME", "disease": "Pancreatic Cancer"}
    assert summarise_civic([item]).startswith("1 accepted items: EID1 prognostic level B")
    rec = {"uid": "1", "germline": "Pathogenic", "oncogenicity": "", "somatic_impact": ""}
    assert summarise_clinvar([rec]) == "VCV1: germline Pathogenic; oncogenicity -; somatic impact -"
