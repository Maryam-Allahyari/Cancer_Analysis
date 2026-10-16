"""Build mutation features: a patient x driver gene table (1 = mutated) plus TMB.

Output: mutation_features.csv, one row per patient.
"""

import argparse
import logging

import numpy as np
import pandas as pd

log = logging.getLogger("features")


def read_gene_list(path):
    """Read one gene per line, skipping comments and blank lines."""
    with open(path) as f:
        genes = [line.strip() for line in f]
    return [g for g in genes if g and not g.startswith("#")]


def mutation_matrix(coding, patients, genes):
    """1 if the patient has at least one protein-changing mutation in the gene, else 0."""
    coding = coding[coding["gene_id"].isin(genes)]
    mat = pd.crosstab(coding["patient_id"], coding["gene_id"])
    mat = (mat > 0).astype(int)

    # patients with no driver mutations and genes never mutated still need a row/column
    mat = mat.reindex(index=patients, columns=genes, fill_value=0)
    mat.columns = ["mut_" + g for g in mat.columns]
    return mat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutations", required=True, help="mutations_coding.csv from clean.py")
    parser.add_argument("--tmb", required=True, help="tmb.csv from clean.py")
    parser.add_argument("--clinical", required=True, help="clinical.csv from clean.py")
    parser.add_argument("--genes", required=True, help="text file with driver genes")
    parser.add_argument("--out", default="Data/Processed/mutation_features.csv")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")

    patients = pd.read_csv(args.clinical, index_col=0, dtype={"patient_id": str}).index.astype(str)
    coding = pd.read_csv(args.mutations, dtype={"patient_id": str})
    tmb = pd.read_csv(args.tmb, index_col=0)
    tmb.index = tmb.index.astype(str)
    genes = read_gene_list(args.genes)

    feats = mutation_matrix(coding, patients, genes)
    # TMB is very skewed (a few hypermutated tumours), so use the log
    feats["log10_tmb"] = np.log10(tmb.loc[patients, "tmb_per_mb"] + 0.01)

    never = (feats.filter(like="mut_").sum() == 0).sum()
    log.info("%d patients x %d driver genes (%d never mutated) + TMB",
             feats.shape[0], len(genes), never)
    top = feats.filter(like="mut_").mean().sort_values(ascending=False).head(5)
    log.info("most often mutated:\n%s", top.round(2).to_string())

    feats.index.name = "patient_id"
    feats.to_csv(args.out)


if __name__ == "__main__":
    main()
