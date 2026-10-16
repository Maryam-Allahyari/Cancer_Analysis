"""Exploratory analysis: PCA of expression, liver signal by biopsy site, and simple stats.

Outputs go to --out-dir: figures (PNG) and stats.txt.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no screen needed (pipeline / Docker)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from markers import LIVER_GENES


def top_variable_genes(expr, n):
    return expr.loc[:, expr.var().sort_values(ascending=False).index[:n]]


def plot_pca(expr, clin, out_dir, n_genes=2000):
    x = StandardScaler().fit_transform(top_variable_genes(expr, n_genes))
    pca = PCA(n_components=10, random_state=0)
    pcs = pca.fit_transform(x)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, col in zip(axes, ["cancer_type", "biopsy_group"]):
        groups = clin[col]
        # biopsy site has many small groups, so lump the rare ones together
        if col == "biopsy_group":
            common = groups.value_counts().index[:5]
            groups = groups.where(groups.isin(common), "Other")
        for g in sorted(groups.unique()):
            mask = (groups == g).values
            ax.scatter(pcs[mask, 0], pcs[mask, 1], s=14, alpha=0.75, label=g)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
        ax.set_title(f"PCA coloured by {col.replace('_', ' ')}")
        ax.legend(fontsize=8, markerscale=1.5)
    fig.tight_layout()
    fig.savefig(out_dir / "pca.png", dpi=150)
    plt.close(fig)

    # scree plot
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.bar(range(1, 11), pca.explained_variance_ratio_ * 100)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("% variance explained")
    fig.tight_layout()
    fig.savefig(out_dir / "pca_scree.png", dpi=150)
    plt.close(fig)
    return pca.explained_variance_ratio_


def plot_liver_score(expr, clin, out_dir):
    """Average expression of liver genes, split by where the biopsy was taken."""
    score = expr[list(LIVER_GENES)].mean(axis=1)
    is_liver = clin["biopsy_group"].eq("Liver")

    fig, ax = plt.subplots(figsize=(8, 4))
    types = sorted(clin["cancer_type"].unique())
    data, labels = [], []
    for t in types:
        for site, mask in [("liver", is_liver), ("other", ~is_liver)]:
            vals = score[(clin["cancer_type"] == t) & mask]
            if len(vals) > 0:
                data.append(vals)
                labels.append(f"{t}\n{site}")
    ax.boxplot(data, tick_labels=labels)
    ax.set_ylabel("Liver gene score\n(mean log2 TPM+1)")
    ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "liver_score_by_biopsy.png", dpi=150)
    plt.close(fig)

    return stats.mannwhitneyu(score[is_liver], score[~is_liver])


def tmb_tests(muts, clin):
    """Does TMB differ between cancer types? Kruskal-Wallis, then pairwise tests with FDR."""
    tmb = muts["log10_tmb"]
    groups = {t: tmb[clin["cancer_type"] == t] for t in clin["cancer_type"].unique()}
    kw = stats.kruskal(*groups.values())

    pairs = []
    names = sorted(groups)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            p = stats.mannwhitneyu(groups[a], groups[b]).pvalue
            pairs.append((a, b, p))
    pairs = pd.DataFrame(pairs, columns=["type_a", "type_b", "p"])
    pairs["p_fdr"] = stats.false_discovery_control(pairs["p"])
    return kw, pairs.sort_values("p_fdr")


def tp53_test(muts, clin, a="BRCA", b="COLO"):
    """Is TP53 mutated more often in one cancer type than another? Fisher's exact test."""
    sub = clin["cancer_type"].isin([a, b])
    table = pd.crosstab(clin.loc[sub, "cancer_type"], muts.loc[sub, "mut_TP53"])
    return table, stats.fisher_exact(table.values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--mutation-features", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--out-dir", default="results/explore")
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    clin = pd.read_csv(args.clinical, index_col=0)
    clin.index = clin.index.astype(str)
    expr = pd.read_parquet(args.expression).loc[clin.index]
    muts = pd.read_csv(args.mutation_features, index_col=0)
    muts.index = muts.index.astype(str)
    muts = muts.loc[clin.index]

    var_ratio = plot_pca(expr, clin, out)
    liver_test = plot_liver_score(expr, clin, out)
    kw, pairs = tmb_tests(muts, clin)
    table, fisher = tp53_test(muts, clin)

    median_tmb = (10 ** muts["log10_tmb"]).groupby(clin["cancer_type"]).median()

    with open(out / "stats.txt", "w") as f:
        f.write("Patients per cancer type\n")
        f.write(clin["cancer_type"].value_counts().to_string() + "\n\n")
        f.write(f"PCA: PC1 explains {var_ratio[0]:.1%}, PC2 {var_ratio[1]:.1%} of variance "
                f"(top 2000 variable genes, scaled)\n\n")
        f.write("Liver gene score, liver biopsies vs other biopsies (Mann-Whitney U)\n")
        f.write(f"U = {liver_test.statistic:.0f}, p = {liver_test.pvalue:.2e}\n\n")
        f.write("Median TMB (mutations per Mb) by cancer type\n")
        f.write(median_tmb.sort_values().round(2).to_string() + "\n")
        f.write(f"Kruskal-Wallis on log10 TMB: H = {kw.statistic:.1f}, p = {kw.pvalue:.2e}\n")
        f.write("Pairwise Mann-Whitney with Benjamini-Hochberg FDR (top 10):\n")
        f.write(pairs.head(10).to_string(index=False) + "\n\n")
        f.write("TP53 mutation, BRCA vs COLO (Fisher's exact test)\n")
        f.write(table.to_string() + "\n")
        f.write(f"odds ratio = {fisher.statistic:.2f}, p = {fisher.pvalue:.2e}\n")
    print(open(out / "stats.txt").read())


if __name__ == "__main__":
    main()
