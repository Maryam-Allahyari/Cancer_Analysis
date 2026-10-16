"""Explain the final model with SHAP and check whether biopsy site affects predictions.

Outputs (in --out-dir): shap_top_features.png, shap_by_class.csv,
shap_waterfall_correct.png, shap_waterfall_wrong.png, biopsy_check.txt
"""

import argparse
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap

# train.py defines the custom transformer the saved model needs
from train import TopVarianceGenes  # noqa: F401
from markers import LIVER_GENES, pretty_name


def accuracy_by_group(pred):
    """Accuracy for liver biopsies vs all other biopsies."""
    group = np.where(pred["biopsy_group"] == "Liver", "liver biopsy", "other biopsy")
    correct = pred["true"] == pred["predicted"]
    return correct.groupby(group).agg(["mean", "size"]).rename(
        columns={"mean": "accuracy", "size": "n"})


def waterfall(expl, i, cls_idx, title, path):
    plt.figure()
    shap.plots.waterfall(expl[i, :, cls_idx], max_display=12, show=False)
    plt.title(title, fontsize=10)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="final_model.joblib from train.py")
    parser.add_argument("--expression", required=True)
    parser.add_argument("--mutation-features", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--test-predictions", required=True)
    parser.add_argument("--oof-predictions", required=True)
    parser.add_argument("--out-dir", default="results/explain")
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    saved = joblib.load(args.model)
    model, classes = saved["model"], saved["classes"]

    clin = pd.read_csv(args.clinical, index_col=0)
    clin.index = clin.index.astype(str)
    muts = pd.read_csv(args.mutation_features, index_col=0)
    muts.index = muts.index.astype(str)
    X = pd.read_parquet(args.expression).loc[clin.index].join(muts.loc[clin.index])

    # SHAP works on the tree model, so transform the data with the fitted preprocessing first
    prep, clf = model.named_steps["prep"], model.named_steps["clf"]
    names = [pretty_name(f) for f in prep.get_feature_names_out()]
    Xt = pd.DataFrame(prep.transform(X), index=X.index, columns=names)

    explainer = shap.TreeExplainer(clf)
    expl = explainer(Xt)  # values shape: patients x features x classes
    sv = expl.values

    # ---- global importance: mean |SHAP| summed over classes ----
    importance = pd.Series(np.abs(sv).mean(axis=0).sum(axis=1), index=names)
    top = importance.sort_values(ascending=False).head(20)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.barh(top.index[::-1], top.values[::-1])
    ax.set_xlabel("mean |SHAP value| (summed over classes)")
    ax.set_title("Top 20 features")
    fig.tight_layout()
    fig.savefig(out / "shap_top_features.png", dpi=150)
    plt.close(fig)

    # ---- top features pushing towards each cancer type ----
    rows = []
    for k, cls in enumerate(classes):
        in_class = (clin["cancer_type"] == cls).values
        mean_shap = pd.Series(sv[in_class, :, k].mean(axis=0), index=names)
        for feat, val in mean_shap.sort_values(ascending=False).head(10).items():
            rows.append({"cancer_type": cls, "feature": feat, "mean_shap": round(val, 4)})
    pd.DataFrame(rows).to_csv(out / "shap_by_class.csv", index=False)

    # ---- one correct and one wrong test patient ----
    test_pred = pd.read_csv(args.test_predictions, index_col=0)
    test_pred.index = test_pred.index.astype(str)
    pos = {pid: i for i, pid in enumerate(Xt.index)}
    for label, rows_ in [("correct", test_pred[test_pred.true == test_pred.predicted]),
                         ("wrong", test_pred[test_pred.true != test_pred.predicted])]:
        if len(rows_) == 0:
            continue
        pid = rows_.index[0]
        k = list(classes).index(rows_.loc[pid, "predicted"])
        title = (f"Patient {pid}: true {rows_.loc[pid, 'true']}, predicted "
                 f"{rows_.loc[pid, 'predicted']} (biopsy: {rows_.loc[pid, 'biopsy_group']})")
        waterfall(expl, pos[pid], k, title, out / f"shap_waterfall_{label}.png")

    # ---- biopsy site check ----
    oof = pd.read_csv(args.oof_predictions, index_col=0)
    liver_feats = [LIVER_GENES[g] for g in LIVER_GENES if LIVER_GENES[g] in names]
    is_liver = (clin["biopsy_group"] == "Liver").values

    lines = ["Accuracy by biopsy site, training set (out-of-fold predictions)",
             accuracy_by_group(oof).round(3).to_string(), "",
             "Accuracy by biopsy site, test set",
             accuracy_by_group(test_pred).round(3).to_string(), "",
             "What liver biopsies were wrongly predicted as (training, out-of-fold):",
             oof[(oof.biopsy_group == "Liver") & (oof.true != oof.predicted)]
             .groupby(["true", "predicted"]).size().to_string(), ""]

    if liver_feats:
        total = np.abs(sv).sum(axis=(1, 2))
        liver = np.abs(sv[:, [names.index(f) for f in liver_feats], :]).sum(axis=(1, 2))
        share = liver / total
        lines += [f"Liver genes used by the model: {', '.join(liver_feats)}",
                  "Share of each patient's total |SHAP| that comes from liver genes:",
                  f"  liver biopsies: median {np.median(share[is_liver]):.2%}",
                  f"  other biopsies: median {np.median(share[~is_liver]):.2%}"]
    else:
        lines.append("None of the liver marker genes were among the selected features.")

    text = "\n".join(lines)
    (out / "biopsy_check.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
