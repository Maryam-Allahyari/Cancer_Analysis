# To check there's no shared patients in both test and train groups
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from train import TopVarianceGenes, bootstrap_ci


def test_no_patient_in_both_sets():
    X = pd.DataFrame(np.zeros((40, 3)), index=[f"p{i}" for i in range(40)])
    y = np.repeat([0, 1, 2, 3], 10)
    X_train, X_test, _, _ = train_test_split(X, y, test_size=0.2, stratify=y, random_state=0)
    assert not set(X_train.index) & set(X_test.index)
    assert len(X_train) + len(X_test) == len(X)


def test_top_variance_picks_most_variable():
    X = pd.DataFrame({"flat": [1, 1, 1, 1], "some": [1, 2, 1, 2], "lots": [0, 10, 0, 10]})
    sel = TopVarianceGenes(n_genes=2).fit(X)
    assert list(sel.transform(X).columns) == ["lots", "some"]


def test_top_variance_only_uses_training_rows():
    # "a" varies only in the last two rows; if those are not used for fitting, "b" wins
    X = pd.DataFrame({"a": [1, 1, 1, 100, -100], "b": [0, 2, 0, 0, 0]})
    sel = TopVarianceGenes(n_genes=1).fit(X.iloc[:3])
    assert list(sel.columns_) == ["b"]


def test_bootstrap_ci_perfect_predictions():
    y = np.array([0, 1, 2] * 10)
    low, high = bootstrap_ci(y, y, n_boot=50)
    assert low == high == 1.0
