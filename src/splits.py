import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

LINEAGE_COLUMN = "cutoff_1168"
TRAIN_YEARS = (2001, 2008)
TEST_YEARS = (2009, 2012)


def count_time_split(year, labels, antibiotics, train_years, test_years):
    is_train = year.between(*train_years)
    is_test = year.between(*test_years)
    rows = []
    for antibiotic in antibiotics:
        train_labels = labels.loc[is_train, antibiotic].dropna()
        test_labels = labels.loc[is_test, antibiotic].dropna()
        rows.append({
            "antibiotic": antibiotic,
            "train_R": (train_labels == "R").sum(),
            "train_S": (train_labels == "S").sum(),
            "train_pct_R": round((train_labels == "R").mean() * 100, 1),
            "test_R": (test_labels == "R").sum(),
            "test_S": (test_labels == "S").sum(),
            "test_pct_R": round((test_labels == "R").mean() * 100, 1),
        })
    return pd.DataFrame(rows).set_index("antibiotic")


def random_split(labels, antibiotic):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    folds = []
    y = labels[antibiotic].dropna()
    for train_pos, test_pos in skf.split(np.zeros(len(y)), y):
        folds.append((y.index[train_pos], y.index[test_pos]))
    return folds


def lineage_split(labels, clusters, antibiotic):
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    folds = []
    y = labels[antibiotic].dropna()
    groups = clusters.loc[y.index, LINEAGE_COLUMN]
    for train_pos, test_pos in sgkf.split(np.zeros(len(y)), y, groups):
        folds.append((y.index[train_pos], y.index[test_pos]))
    return folds


def time_split(labels, year, antibiotic):
    y = labels[antibiotic].dropna()
    y_year = year.loc[y.index]
    in_train = y_year.between(*TRAIN_YEARS).fillna(False)
    in_test = y_year.between(*TEST_YEARS).fillna(False)
    return [(in_train[in_train].index, in_test[in_test].index)]