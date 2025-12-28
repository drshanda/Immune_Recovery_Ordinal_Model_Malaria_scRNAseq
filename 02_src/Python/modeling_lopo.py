# ============================================================
# HostScope — Modeling under LOPO (Ordinal + Pairwise)
# ============================================================
# This script is extracted verbatim (by cell) from:
#   LOPO_train_test_split.ipynb
#
# Purpose:
#   - Run LOPO splitting
#   - Fit pairwise logistic regression (Day0 vs Day28)
#   - Fit ordinal logistic regression (Day0 < Day7 < Day28)
#   - Compute and save core performance metrics/figures
#
# Notes:
#   - Code is included exactly as in the notebook cells referenced.
#   - Redundancies across scripts are intentional for reproducibility.
# ============================================================


# --------------------
# Notebook cell 1
# --------------------

# ============================================================
# HostScope Proof-of-Concept Modeling Script (UPDATED)
# ============================================================
# - LOPO splitting
# - Pairwise logistic regression (D0 vs D28)
# - Ordinal logistic regression (D0 < D7 < D28)
# - Metrics:
#     * Balanced Accuracy (both)
#     * AUC-ROC (pairwise)
#     * MAE (ordinal)
# - Saves:
#     * metrics CSV
#     * confusion matrix plots
#     * classification reports (CSV)
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    roc_auc_score,
    confusion_matrix,
    classification_report,
)

from statsmodels.miscmodels.ordinal_model import OrderedModel
from sklearn.preprocessing import label_binarize


# ============================================================
# Setup
# ============================================================

RESULTS_DIR = "/Users/lashandawilliams/Immune_Recovery_Ordinal_Model_Malaria_scRNAseq/03_results/"
os.makedirs(RESULTS_DIR, exist_ok=True)


# ============================================================
# Utility functions
# ============================================================

def select_feature_columns(df, exclude):
    return [
        c for c in df.columns
        if c not in exclude and pd.api.types.is_numeric_dtype(df[c])
    ]


def lopo_splits(df, donor_col="ID"):
    logo = LeaveOneGroupOut()
    groups = df[donor_col].astype(str).values
    X_dummy = np.zeros((len(df), 1))
    return list(logo.split(X_dummy, groups=groups))


def plot_confusion_matrix(
    cm,
    labels,
    title,
    out_path,
    normalize=False,
    figsize=(6, 6),
    cmap="Blues"
):
    """
    Plot a readable confusion matrix with large annotations.

    Parameters
    ----------
    cm : array-like (n_classes, n_classes)
        Confusion matrix
    labels : list of str
        Class labels in display order
    title : str
        Plot title
    out_path : str
        Where to save PNG
    normalize : bool
        If True, normalize rows to proportions
    figsize : tuple
        Figure size
    cmap : str
        Matplotlib colormap
    """

    cm = np.array(cm)

    if normalize:
        cm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    plt.figure(figsize=figsize)
    ax = sns.heatmap(
        cm,
        annot=True,
        fmt=".2f" if normalize else "d",
        cmap=cmap,
        cbar=True,
        square=True,
        linewidths=1.5,
        linecolor="white",
        annot_kws={
            "size": 18,
            "weight": "bold",
            "color": "black"
        }
    )

    ax.set_xticklabels(labels, fontsize=14)
    ax.set_yticklabels(labels, fontsize=14, rotation=0)

    ax.set_xlabel("Predicted", fontsize=16, labelpad=12)
    ax.set_ylabel("True", fontsize=16, labelpad=12)

    ax.set_title(title, fontsize=18, pad=16)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()



# ============================================================
# Pairwise Logistic Regression
# ============================================================


def run_pairwise_logistic(df):

    exclude_cols = {"donor_id", "timepoint", "y_binary", "y_stage"}
    feature_cols = select_feature_columns(df, exclude_cols)

    splits = lopo_splits(df)

    oof_rows = []
    fold_metrics = []

    for fold, (train_idx, test_idx) in enumerate(splits):

        train_df = df.iloc[train_idx]
        test_df  = df.iloc[test_idx]

        X_train = train_df[feature_cols].values
        y_train = train_df["y_binary"].values

        X_test  = test_df[feature_cols].values
        y_test  = test_df["y_binary"].values

        model = Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                ("clf", LogisticRegression(
                    penalty="l2",
                    C=1.0,
                    solver="liblinear",
                    max_iter=1000,
                )),
            ]
        )

        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)
        y_prob = model.predict_proba(X_test)[:, 1]

        # Store OOF predictions
        for yt, yp, yp_prob in zip(y_test, y_pred, y_prob):
            oof_rows.append({
                "fold": fold,
                "y_true": int(yt),
                "y_pred": int(yp),
                "y_prob": float(yp_prob),
            })

        fold_metrics.append({
            "model": "pairwise_logistic",
            "fold": fold,
            "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
            "auc": roc_auc_score(y_test, y_prob),
            "mae": float(np.mean(np.abs(y_pred - y_test))),
        })

    # =========================
    # OVERALL LOPO METRICS
    # =========================
    oof_df = pd.DataFrame(oof_rows)

    overall_metrics = {
        "model": "pairwise_logistic",
        "balanced_accuracy": float(
            balanced_accuracy_score(oof_df["y_true"], oof_df["y_pred"])
        ),
        "auc": float(
            roc_auc_score(oof_df["y_true"], oof_df["y_prob"])
        ),
        "mae": float(
            np.mean(np.abs(oof_df["y_pred"] - oof_df["y_true"]))
        ),
    }


    # =========================
    # CONFUSION MATRIX
    # =========================
    cm = confusion_matrix(oof_df["y_true"], oof_df["y_pred"])

    plot_confusion_matrix(
        cm,
        labels=["D28", "D0"],
        title="Pairwise Logistic (D0 vs D28)",
        out_path=f"{RESULTS_DIR}/figures/modeling/confusion_pairwise.png",
        normalize=False
    )

    plot_confusion_matrix(
        cm,
        labels=["D28", "D0"],
        title="Pairwise Logistic (D0 vs D28) — Normalized",
        out_path=f"{RESULTS_DIR}/figures/modeling/confusion_pairwise_normalized.png",
        normalize=True
    )

    # =========================
    # CLASSIFICATION REPORT
    # =========================
    report = classification_report(
        oof_df["y_true"],
        oof_df["y_pred"],
        target_names=["D28", "D0"],
        output_dict=True
    )

    pd.DataFrame(report).T.to_csv(
        f"{RESULTS_DIR}/tables/classification_report_pairwise.csv"
    )

    return fold_metrics, overall_metrics

# ============================================================
# Ordinal Logistic Regression
# ============================================================


def run_ordinal_logistic(df, return_artifacts=False):

    exclude_cols = {"ID", "timepoint", "y_binary", "y_stage"}
    feature_cols = select_feature_columns(df, exclude_cols)

    splits = lopo_splits(df)

    oof_rows = []
    fold_metrics = []

    for fold, (train_idx, test_idx) in enumerate(splits):

        train_df = df.iloc[train_idx].copy()
        test_df  = df.iloc[test_idx].copy()

        scaler = StandardScaler()
        X_train = scaler.fit_transform(train_df[feature_cols].values)
        X_test  = scaler.transform(test_df[feature_cols].values)

        y_train = train_df["y_stage"].values
        y_test  = test_df["y_stage"].values

        model = OrderedModel(y_train, X_train, distr="logit")
        result = model.fit(method="bfgs", disp=False)

        probs = result.predict(X_test)        # (n, 3)
        y_pred = np.argmax(probs, axis=1)

        for yt, yp, pr in zip(y_test, y_pred, probs):
            oof_rows.append({
                "model": "ordinal_logistic",
                "fold": fold,
                "y_true": int(yt),
                "y_pred": int(yp),
                "probs": pr,
            })

        fold_metrics.append({
            "model": "ordinal_logistic",
            "fold": fold,
            "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
            "auc": np.nan,  # computed globally
            "mae": float(np.mean(np.abs(y_pred - y_test))),
        })

    oof_df = pd.DataFrame(oof_rows)

    # --- Overall metrics ---
    y_true_bin = label_binarize(oof_df["y_true"], classes=[0, 1, 2])
    y_prob_mat = np.vstack(oof_df["probs"].values)

    overall_metrics = {
        "model": "ordinal_logistic",
        "balanced_accuracy": float(
            balanced_accuracy_score(oof_df["y_true"], oof_df["y_pred"])
        ),
        "auc": float(
            roc_auc_score(
                y_true_bin,
                y_prob_mat,
                average="macro",
                multi_class="ovr",
            )
        ),
        "mae": float(
            np.mean(np.abs(oof_df["y_pred"] - oof_df["y_true"]))
        ),
    }

    # --- Confusion matrix ---
    cm = confusion_matrix(
    oof_df["y_true"],
    oof_df["y_pred"],
    labels=[0, 1, 2]
)

    plot_confusion_matrix(
    cm,
    labels=["D0", "D7", "D28"],
    title="Ordinal Logistic (D0–D7–D28)",
    out_path=f"{RESULTS_DIR}/figures/modeling/confusion_ordinal.png",
    normalize=False,
    figsize=(7, 7)
)

    plot_confusion_matrix(
    cm,
    labels=["D0", "D7", "D28"],
    title="Ordinal Logistic (D0–D7–D28) — Normalized",
    out_path=f"{RESULTS_DIR}/figures/modeling/confusion_ordinal_normalized.png",
    normalize=True,
    figsize=(7, 7)
)

    # --- Classification report ---
    report = classification_report(
        oof_df["y_true"],
        oof_df["y_pred"],
        output_dict=True,
    )
    pd.DataFrame(report).T.to_csv(
        f"{RESULTS_DIR}/tables/classification_report_ordinal.csv"
    )

    if return_artifacts:
        return {
            "fold_metrics": fold_metrics,
            "overall_metrics": overall_metrics,
            "X": X_all,
            "y": y_all,
            "feature_names": feature_cols,
            "model": final_model,
            "result": final_result,
            "scaler": scaler_final,
        }

    return fold_metrics, overall_metrics

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    df_pairwise = pd.read_csv(
        "/Users/lashandawilliams/Immune_Recovery_Ordinal_Model_Malaria_scRNAseq/01_data/processed/hostscope_features_pairwise_D0_vs_D28.csv"
    )

    df_ordinal = pd.read_csv(
        "/Users/lashandawilliams/Immune_Recovery_Ordinal_Model_Malaria_scRNAseq/01_data/processed/hostscope_features_ordinal_D0_D7_D28.csv"
    )


    pairwise_folds, pairwise_overall = run_pairwise_logistic(df_pairwise)
    ordinal_folds, ordinal_overall = run_ordinal_logistic(df_ordinal)

    metrics_df = pd.DataFrame(pairwise_folds + ordinal_folds)
    metrics_df.to_csv(f"{RESULTS_DIR}/tables/hostscope_fold_metrics.csv", index=False)

    overall_df = pd.DataFrame([pairwise_overall, ordinal_overall])
    overall_df.to_csv(f"{RESULTS_DIR}/tables/hostscope_overall_metrics.csv", index=False)


