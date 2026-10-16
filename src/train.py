"""Train and compare models for predicting cancer type.

Two models (logistic regression, XGBoost) x three feature sets
(expression, mutations, both) are compared with stratified 5-fold CV on the
training set. Each model is then refitted on the whole training set and scored
once on the held-out test set.

Outputs (in --out-dir): cv_results.csv, test_results.csv, test_predictions.csv,
oof_predictions.csv, confusion_matrix.png, final_model.joblib
"""

import argparse
import json
import warnings
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import loguniform, randint, uniform
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (ConfusionMatrixDisplay, accuracy_score,
                             balanced_accuracy_score, f1_score)
from sklearn.model_selection import (RandomizedSearchCV, StratifiedKFold,
                                     cross_val_predict, cross_validate,
                                     train_test_split)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from xgboost import XGBClassifier


class TopVarianceGenes(BaseEstimator, TransformerMixin):
    """Keep the n columns with the highest variance.

    Written as a transformer so it is fitted on the training folds only (no leakage).
    """

    def __init__(self, n_genes=2000):
        self.n_genes = n_genes

    def fit(self, X, y=None):
        X = pd.DataFrame(X)
        self.columns_ = X.var().sort_values(ascending=False).index[:self.n_genes]
        return self

    def transform(self, X):
        return pd.DataFrame(X)[self.columns_]

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.columns_, dtype=object)


def make_preprocessor(feature_set, expr_cols, mut_cols, n_genes):
    expr_steps = Pipeline([("top_var", TopVarianceGenes(n_genes)),
                           ("scale", StandardScaler())])
    parts = []
    if feature_set in ("expression", "both"):
        parts.append(("expr", expr_steps, expr_cols))
    if feature_set in ("mutations", "both"):
        parts.append(("mut", StandardScaler(), mut_cols))
    # verbose_feature_names_out=False keeps plain gene names for SHAP later
    return ColumnTransformer(parts, verbose_feature_names_out=False)


def make_model(name, seed):
    if name == "logreg":
        # L2 penalty; class_weight balances the small classes
        return LogisticRegression(C=0.05, max_iter=5000, class_weight="balanced")
    if name == "xgboost":
        return XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05,
                             subsample=0.8, colsample_bytree=0.3,
                             tree_method="hist", n_jobs=-1, random_state=seed)
    raise ValueError(name)


def bootstrap_ci(y_true, y_pred, n_boot=1000, seed=0):
    """95% CI for balanced accuracy by resampling test patients."""
    rng = np.random.default_rng(seed)
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    scores = []
    # a resample can miss a small class completely, which makes sklearn warn every time
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        for _ in range(n_boot):
            idx = rng.integers(0, len(y_true), len(y_true))
            scores.append(balanced_accuracy_score(y_true[idx], y_pred[idx]))
    return np.percentile(scores, [2.5, 97.5])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--mutation-features", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--out-dir", default="results/models")
    parser.add_argument("--n-genes", type=int, default=2000)
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tune-iter", type=int, default=10,
                        help="random search iterations for the final XGBoost model (0 = no tuning)")
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- load and join everything on patient ID ----
    clin = pd.read_csv(args.clinical, index_col=0)
    clin.index = clin.index.astype(str)
    expr = pd.read_parquet(args.expression).loc[clin.index]
    muts = pd.read_csv(args.mutation_features, index_col=0)
    muts.index = muts.index.astype(str)
    X = expr.join(muts.loc[clin.index])
    expr_cols, mut_cols = list(expr.columns), list(muts.columns)

    # XGBoost wants labels as 0..k-1
    encoder = LabelEncoder()
    y = encoder.fit_transform(clin["cancer_type"])
    classes = encoder.classes_

    # ---- hold out a test set, used only once at the very end ----
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, stratify=y, random_state=args.seed)
    print(f"train: {len(X_train)} patients, test: {len(X_test)} patients, {len(classes)} classes")
    assert not set(X_train.index) & set(X_test.index), "patient in both train and test!"

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=args.seed)
    scoring = ["balanced_accuracy", "accuracy", "f1_macro"]

    # ---- compare models and feature sets with cross-validation ----
    cv_rows, test_rows = [], []
    fitted = {}
    for feature_set in ["expression", "mutations", "both"]:
        for model_name in ["logreg", "xgboost"]:
            pipe = Pipeline([
                ("prep", make_preprocessor(feature_set, expr_cols, mut_cols, args.n_genes)),
                ("clf", make_model(model_name, args.seed)),
            ])
            res = cross_validate(pipe, X_train, y_train, cv=cv, scoring=scoring, n_jobs=1)
            row = {"model": model_name, "features": feature_set}
            for s in scoring:
                row[f"cv_{s}_mean"] = res[f"test_{s}"].mean()
                row[f"cv_{s}_sd"] = res[f"test_{s}"].std()
            cv_rows.append(row)
            print(f"{model_name:8s} {feature_set:10s} CV balanced acc = "
                  f"{row['cv_balanced_accuracy_mean']:.3f} +/- {row['cv_balanced_accuracy_sd']:.3f}")

            pipe.fit(X_train, y_train)
            fitted[(model_name, feature_set)] = pipe

    cv_results = pd.DataFrame(cv_rows)
    cv_results.to_csv(out / "cv_results.csv", index=False)

    # ---- tune the XGBoost model on both feature sets (training data only) ----
    final_key = ("xgboost", "both")
    if args.tune_iter > 0:
        search = RandomizedSearchCV(
            fitted[final_key],
            param_distributions={
                "clf__n_estimators": randint(150, 600),
                "clf__max_depth": randint(2, 6),
                "clf__learning_rate": loguniform(0.01, 0.2),
                "clf__colsample_bytree": uniform(0.1, 0.5),
                "clf__min_child_weight": randint(1, 6),
            },
            n_iter=args.tune_iter, scoring="balanced_accuracy", cv=3,
            random_state=args.seed, n_jobs=1)
        search.fit(X_train, y_train)
        print(f"tuned XGBoost (3-fold) balanced acc = {search.best_score_:.3f}")
        print("best params:", search.best_params_)
        fitted[("xgboost_tuned", "both")] = search.best_estimator_
        final_key = ("xgboost_tuned", "both")
        with open(out / "best_params.json", "w") as f:
            json.dump({k: float(v) for k, v in search.best_params_.items()}, f, indent=2)

    # ---- score every model once on the test set ----
    for (model_name, feature_set), pipe in fitted.items():
        pred = pipe.predict(X_test)
        low, high = bootstrap_ci(y_test, pred, seed=args.seed)
        test_rows.append({
            "model": model_name, "features": feature_set,
            "test_balanced_accuracy": balanced_accuracy_score(y_test, pred),
            "ci_low": low, "ci_high": high,
            "test_accuracy": accuracy_score(y_test, pred),
            "test_f1_macro": f1_score(y_test, pred, average="macro"),
        })
    test_results = pd.DataFrame(test_rows)
    test_results.to_csv(out / "test_results.csv", index=False)
    print("\nTest set results:")
    print(test_results.round(3).to_string(index=False))

    # ---- save outputs of the final model for the explain step ----
    final = fitted[final_key]
    test_pred = pd.DataFrame({
        "true": classes[y_test], "predicted": classes[final.predict(X_test)],
        "biopsy_group": clin.loc[X_test.index, "biopsy_group"],
    }, index=X_test.index)
    test_pred.to_csv(out / "test_predictions.csv")

    # out-of-fold predictions on the training set give more patients
    # for checking accuracy by biopsy site than the small test set alone
    oof = cross_val_predict(final, X_train, y_train, cv=cv)
    pd.DataFrame({
        "true": classes[y_train], "predicted": classes[oof],
        "biopsy_group": clin.loc[X_train.index, "biopsy_group"],
    }, index=X_train.index).to_csv(out / "oof_predictions.csv")

    fig, ax = plt.subplots(figsize=(7, 6))
    ConfusionMatrixDisplay.from_predictions(test_pred["true"], test_pred["predicted"],
                                            labels=classes, ax=ax, colorbar=False)
    ax.set_title(f"Test set, {final_key[0]} ({final_key[1]})")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(out / "confusion_matrix.png", dpi=150)
    plt.close(fig)

    joblib.dump({"model": final, "classes": classes,
                 "train_ids": list(X_train.index), "test_ids": list(X_test.index)},
                out / "final_model.joblib")


if __name__ == "__main__":
    main()
