"""Look up public evidence for the variants in the annotation worksheet (facts only).

Fills the factual columns of the worksheet from open databases, so the curator can focus
on interpretation. Tiers, variant class and "supports true / predicted type" are left
empty on purpose: those calls need a curator.

Sources (no login needed; only gene + protein change are sent, never patient data):
    UniProt     checks the reference amino acid at each position (POG uses Ensembl 69 transcripts)
    CIViC       accepted clinical evidence items for the exact variant
    cBioPortal  how often the gene and the exact variant are mutated in TCGA PanCancer Atlas types
    ClinVar     germline, oncogenicity and somatic clinical-impact classifications
COSMIC needs a login, so it is not queried.

Outputs (in the same git-ignored folder as the worksheet):
    annotation_worksheet_evidence.csv  the worksheet with the factual columns filled in
    evidence_frequencies.csv           gene and exact-variant frequency per TCGA cancer type
"""

import argparse
import logging
import re
import time
from pathlib import Path

import pandas as pd
import requests

log = logging.getLogger("gather_evidence")

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "*": "*",
}

# Ensembl 69 gene names in the POG file that have since been renamed (HGNC)
RENAMED = {"MLL3": "KMT2C", "MLL2": "KMT2D"}

# closest TCGA PanCancer Atlas cohorts for each POG cancer-type label
POG_TO_TCGA = {
    "BRCA": ["brca"], "COLO": ["coadread"], "LUNG": ["luad", "lusc"], "PANC": ["paad"],
    "SARC": ["sarc"], "OV": ["ov"], "CNS-PNS": ["gbm", "lgg"],
}

CBIO = "https://www.cbioportal.org/api"
CIVIC = "https://civicdb.org/api/graphql"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
TIMEOUT = 120

CIVIC_QUERY = """
query($genes: [String!], $variant: String) {
  genes(entrezSymbols: $genes) { nodes { variants(name: $variant) { nodes {
    name
    ... on GeneVariant { molecularProfiles { nodes { id name
      evidenceItems(first: 100) { nodes {
        id status evidenceLevel evidenceType significance disease { name } } } } } }
  } } } }
}"""


# ---------- pure helpers (tested offline) ----------

def short_protein(hgvs_p):
    """'p.Gly12Asp' -> ('G', 12, 'D', 'G12D'). Returns None for frameshifts, splice or missing changes."""
    m = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)", str(hgvs_p))
    if not m or m[1] not in AA3_TO_1 or m[3] not in AA3_TO_1:
        return None
    ref, pos, alt = AA3_TO_1[m[1]], int(m[2]), AA3_TO_1[m[3]]
    return ref, pos, alt, f"{ref}{pos}{alt}"


def current_symbol(gene):
    return RENAMED.get(gene, gene)


def frequencies(muts, totals, entrez_id, short):
    """Percent of samples per TCGA type with any non-silent mutation in the gene, and with the exact variant.

    muts: one row per mutation with columns study, sample, entrez, protein_change, type
    totals: Series of sequenced sample counts indexed by study
    """
    gene = muts[(muts["entrez"] == entrez_id) & (muts["type"] != "Silent")]
    gene_n = gene.groupby("study")["sample"].nunique()
    var_n = gene[gene["protein_change"] == short].groupby("study")["sample"].nunique()
    freq = pd.DataFrame({
        "gene_pct": (gene_n.reindex(totals.index, fill_value=0) / totals * 100).round(1),
        "variant_pct": (var_n.reindex(totals.index, fill_value=0) / totals * 100).round(2),
        "variant_samples": var_n.reindex(totals.index, fill_value=0),
    })
    freq.index.name = "study"
    return freq


def type_frequency_text(freq, pog_type):
    """e.g. 'PANC (paad): gene 65.4%, variant 27.37%'."""
    studies = [s for s in POG_TO_TCGA.get(pog_type, []) if s in freq.index]
    if not studies:
        return f"{pog_type}: no matching TCGA cohort"
    parts = [f"{s} gene {freq.loc[s, 'gene_pct']}%, variant {freq.loc[s, 'variant_pct']}%" for s in studies]
    return f"{pog_type} ({'; '.join(parts)})"


def top_types_text(freq, column, n=5):
    top = freq[freq[column] > 0].sort_values(column, ascending=False).head(n)
    return ", ".join(f"{s} {v}%" for s, v in top[column].items()) or "none"


def summarise_civic(items):
    """One line per worksheet cell from CIViC accepted evidence items."""
    if not items:
        return "no accepted CIViC evidence for this exact variant"
    shown = [f"EID{e['id']} {e['evidenceType'].lower()} level {e['evidenceLevel']} "
             f"{e['significance'].lower()} ({e['disease']})" for e in items[:4]]
    return f"{len(items)} accepted items: " + "; ".join(shown)


def summarise_clinvar(records):
    if not records:
        return "not in ClinVar"
    r = records[0]
    return (f"VCV{r['uid']}: germline {r['germline'] or '-'}; oncogenicity {r['oncogenicity'] or '-'}; "
            f"somatic impact {r['somatic_impact'] or '-'}")


# ---------- network lookups ----------

def uniprot_reference(session, gene):
    """Accession and canonical sequence of the reviewed human UniProt entry."""
    r = session.get(UNIPROT, timeout=TIMEOUT, params={
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,sequence", "format": "json", "size": 1})
    r.raise_for_status()
    hit = r.json()["results"][0]
    return hit["primaryAccession"], hit["sequence"]["value"]


def civic_evidence(session, gene, short):
    """Accepted evidence items for the single-variant profile 'GENE SHORT' (combination profiles skipped)."""
    r = session.post(CIVIC, timeout=TIMEOUT,
                     json={"query": CIVIC_QUERY, "variables": {"genes": [gene], "variant": short}})
    r.raise_for_status()
    items = []
    for g in r.json()["data"]["genes"]["nodes"]:
        for v in g["variants"]["nodes"]:
            for mp in v.get("molecularProfiles", {}).get("nodes", []):
                if mp["name"] != f"{gene} {v['name']}":
                    continue
                for e in mp["evidenceItems"]["nodes"]:
                    if e["status"] == "ACCEPTED":
                        items.append({**e, "disease": e["disease"]["name"] if e["disease"] else ""})
    return items


def clinvar_records(session, gene, hgvs_p, short):
    term = f'{gene}[gene] AND "{hgvs_p.replace("*", "Ter")}"'
    r = session.get(f"{EUTILS}/esearch.fcgi", timeout=TIMEOUT,
                    params={"db": "clinvar", "term": term, "retmode": "json"})
    r.raise_for_status()
    ids = r.json()["esearchresult"]["idlist"][:20]
    time.sleep(0.4)  # NCBI allows 3 requests per second without an API key
    if not ids:
        return []
    r = session.get(f"{EUTILS}/esummary.fcgi", timeout=TIMEOUT,
                    params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"})
    r.raise_for_status()
    time.sleep(0.4)
    result = r.json()["result"]
    records = []
    for uid in result["uids"]:
        rec = result[uid]
        changes = [c.strip() for c in rec.get("protein_change", "").split(",")]
        if short not in changes:
            continue
        records.append({
            "uid": uid,
            "germline": rec.get("germline_classification", {}).get("description", ""),
            "oncogenicity": rec.get("oncogenicity_classification", {}).get("description", ""),
            "somatic_impact": rec.get("clinical_impact_classification", {}).get("description", ""),
        })
    return records


def cbioportal_mutations(session, entrez_ids):
    """All mutations in the given genes across TCGA PanCancer Atlas studies, plus sequenced sample counts."""
    r = session.get(f"{CBIO}/studies", params={"keyword": "pan_can_atlas_2018"}, timeout=TIMEOUT)
    r.raise_for_status()
    studies = sorted(s["studyId"] for s in r.json())
    rows, totals = [], {}
    for study in studies:
        short_name = study.split("_")[0]
        samples = session.get(f"{CBIO}/sample-lists/{study}_sequenced", timeout=TIMEOUT)
        samples.raise_for_status()
        totals[short_name] = len(samples.json()["sampleIds"])
        muts = session.post(f"{CBIO}/molecular-profiles/{study}_mutations/mutations/fetch",
                            params={"projection": "SUMMARY"}, timeout=TIMEOUT,
                            json={"entrezGeneIds": entrez_ids, "sampleListId": f"{study}_sequenced"})
        muts.raise_for_status()
        rows += [{"study": short_name, "sample": m["sampleId"], "entrez": m["entrezGeneId"],
                  "protein_change": m.get("proteinChange", ""), "type": m.get("mutationType", "")}
                 for m in muts.json()]
    log.info("cBioPortal: %d mutations across %d TCGA studies", len(rows), len(studies))
    return pd.DataFrame(rows, columns=["study", "sample", "entrez", "protein_change", "type"]), pd.Series(totals)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--worksheet", default="results/annotation/annotation_worksheet.csv")
    parser.add_argument("--out-dir", default="results/annotation")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")
    out = Path(args.out_dir)
    sheet = pd.read_csv(args.worksheet, dtype=str).fillna("")
    session = requests.Session()

    genes = sorted({current_symbol(g) for g in sheet["gene_id"]})
    entrez = {}
    for g in genes:
        r = session.get(f"{CBIO}/genes/{g}", timeout=TIMEOUT)
        r.raise_for_status()
        entrez[g] = r.json()["entrezGeneId"]
    muts, totals = cbioportal_mutations(session, list(entrez.values()))
    references = {g: uniprot_reference(session, g) for g in genes}

    freq_tables = []
    for i, row in sheet.iterrows():
        gene = current_symbol(row["gene_id"])
        parsed = short_protein(row["hgvs_p"])
        if parsed is None:
            sheet.loc[i, "sources"] = "not looked up: no simple protein change (frameshift, splice or missing)"
            continue
        ref, pos, _, short = parsed
        accession, seq = references[gene]
        position_ok = pos <= len(seq) and seq[pos - 1] == ref

        freq = frequencies(muts, totals, entrez[gene], short)
        freq_tables.append(freq.reset_index().assign(gene=gene, hgvs_p=row["hgvs_p"]))
        civic = civic_evidence(session, gene, short)
        clinvar = clinvar_records(session, gene, row["hgvs_p"], short)

        sheet.loc[i, "civic_evidence"] = summarise_civic(civic)
        sheet.loc[i, "clinvar_somatic"] = summarise_clinvar(clinvar)
        sheet.loc[i, "cosmic_or_cbioportal"] = (
            f"TCGA: exact variant in {int(freq['variant_samples'].sum())} samples; "
            f"true type {type_frequency_text(freq, row['true_type'])}; "
            f"predicted type {type_frequency_text(freq, row['predicted_type'])}; COSMIC not checked (login)")
        sheet.loc[i, "typical_cancer_types"] = (
            f"variant top: {top_types_text(freq, 'variant_pct')} | gene top: {top_types_text(freq, 'gene_pct')}")
        sheet.loc[i, "sources"] = (
            f"UniProt {accession} ref {ref}{pos} {'matches' if position_ok else 'DOES NOT MATCH'}; "
            "CIViC; cBioPortal TCGA PanCancer Atlas 2018; ClinVar")
        if gene != row["gene_id"]:
            sheet.loc[i, "sources"] += f"; searched as {gene} (POG name {row['gene_id']})"
        log.info("%s %s: CIViC %d, ClinVar %d, TCGA samples %d, UniProt %s", gene, short, len(civic),
                 len(clinvar), int(freq["variant_samples"].sum()), "ok" if position_ok else "MISMATCH")

    sheet.to_csv(out / "annotation_worksheet_evidence.csv", index=False)
    if freq_tables:
        pd.concat(freq_tables).to_csv(out / "evidence_frequencies.csv", index=False)
    log.info("written to %s (tiers and support calls left for the curator)", out)


if __name__ == "__main__":
    main()
