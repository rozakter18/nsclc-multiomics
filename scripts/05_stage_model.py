#!/usr/bin/env python3
"""Evaluate early-versus-late stage prediction with repeated cross-validation.

Compares models on IDENTICAL cross-validation splits (paired comparison):
    clinical baseline (age + sex)
    MOFA factors (all stable factors)        <- integrated representation
    MOFA factors + clinical
    RNA / methylation / CNV only             <- which layer carries signal?
    all omics concatenated                   <- early integration

Two CV schemes:
    patient : repeated stratified 5-fold, split by patient
    site    : repeated 5-fold holding out whole tissue source sites
              (tests whether performance is a hospital artefact)

Leakage controls:
  * imputation, scaling, feature selection and regularisation tuning are
    all inside the CV pipeline, fitted on training folds only;
  * the factor model uses ALL stable factors, not only the factors that
    correlated with stage in step 4 (choosing them on the full data and then
    testing would leak);
  * MOFA+ was fitted on all patients without labels (transductive, no label
    leakage) - the concatenated-omics model has no such dependence and serves
    as the clean comparison.

Outputs (results/model/):
    cv_results.tsv          per model x scheme: AUROC, balanced accuracy
    cv_per_repeat.tsv       per repeat AUROC (for paired comparisons)
    oof_predictions.tsv     out-of-fold probabilities (patient scheme, repeat 1)
    adjusted_factors.tsv    stage-linked factors adjusted for age and sex
    roc.png                 ROC curves
    model_report.md

Run from the repository root. Preprocessing, feature selection, and model
tuning remain inside each training fold.
"""
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

PROC = Path("data/processed")
MOFA = Path("results/mofa")
OUT = Path("results/model")
VIEWS = ["rna", "methylation", "cnv"]
K_BEST = 500          # features kept per omics model (selected inside folds)
N_FOLDS = 5
LOG = []


def log(msg=""):
    print(msg)
    LOG.append(msg)


# ---- data ----------------------------------------------------------------------
def load():
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0,
                       dtype={"tss": str})
    meta = meta[meta["label_late"].notna()]
    y = meta["label_late"].astype(int)

    clin = pd.DataFrame(index=meta.index)
    clin["age"] = pd.to_numeric(meta["age"], errors="coerce")
    clin["male"] = meta["gender"].astype(str).str.upper().map(
        {"MALE": 1.0, "FEMALE": 0.0})

    ann = pd.read_csv(MOFA / "factor_annotation.tsv", sep="\t")
    stable = ann.loc[ann["stable"].astype(str) == "True", "factor"].tolist()
    Z = pd.read_csv(MOFA / "factors.tsv", sep="\t", index_col=0).loc[meta.index]

    omics = {v: pd.read_parquet(PROC / f"{v}.parquet").loc[meta.index]
             for v in VIEWS}
    allomics = pd.concat({v: omics[v] for v in VIEWS}, axis=1)
    allomics.columns = [f"{v}:{c}" for v, c in allomics.columns]

    feats = {
        "clinical (age+sex)": clin,
        "MOFA factors": Z[stable],
        "MOFA factors + clinical": pd.concat([Z[stable], clin], axis=1),
        "RNA only": omics["rna"],
        "methylation only": omics["methylation"],
        "CNV only": omics["cnv"],
        "all omics concatenated": allomics,
    }
    return meta, y, feats, ann, stable, Z


def pipeline(n_features):
    if n_features <= 50:   # low-dimensional: plain L2 logistic regression
        return make_pipeline(
            SimpleImputer(strategy="median"), StandardScaler(),
            LogisticRegression(class_weight="balanced", max_iter=5000))
    # high-dimensional: select top-K inside the fold, tune L2 strength by
    # inner 3-fold CV (also inside the fold)
    return make_pipeline(
        SimpleImputer(strategy="median"),
        SelectKBest(f_classif, k=min(K_BEST, n_features)),
        StandardScaler(),
        LogisticRegressionCV(Cs=[1e-3, 1e-2, 1e-1, 1.0], cv=3,
                             scoring="roc_auc", class_weight="balanced",
                             max_iter=5000))


# ---- CV splits -----------------------------------------------------------------
def make_splits(meta, y, scheme, repeats, seed=0):
    strat = (y.astype(str) + "_" + meta["histology"]).values
    X0 = np.zeros(len(y))
    out = []
    for r in range(repeats):
        if scheme == "patient":
            cv = StratifiedKFold(N_FOLDS, shuffle=True, random_state=seed + r)
            out.append(list(cv.split(X0, strat)))
        else:
            cv = StratifiedGroupKFold(N_FOLDS, shuffle=True,
                                      random_state=seed + r)
            out.append(list(cv.split(X0, y.values, meta["tss"].values)))
    return out


def oof_predict(X, y, folds):
    p = np.full(len(y), np.nan)
    for tr, te in folds:
        m = pipeline(X.shape[1]).fit(X.iloc[tr], y.iloc[tr])
        p[te] = m.predict_proba(X.iloc[te])[:, 1]
    return p


def evaluate(X, y, meta, splits):
    aucs, baccs, by_hist, first = [], [], {"LUAD": [], "LUSC": []}, None
    for folds in splits:
        p = oof_predict(X, y, folds)
        if first is None:
            first = p
        aucs.append(roc_auc_score(y, p))
        baccs.append(balanced_accuracy_score(y, (p >= 0.5).astype(int)))
        for h in by_hist:
            m = (meta["histology"] == h).values
            by_hist[h].append(roc_auc_score(y[m], p[m]))
    return np.array(aucs), np.array(baccs), by_hist, first


def permutation_p(X, y, folds, n_perm, observed, rng):
    null = []
    for _ in range(n_perm):
        yp = pd.Series(rng.permutation(y.values), index=y.index)
        null.append(roc_auc_score(yp, oof_predict(X, yp, folds)))
    null = np.array(null)
    return (1 + (null >= observed).sum()) / (1 + n_perm), null


# ---- adjusted association of stage-linked factors ----------------------------
def adjusted_factors(y, Z, meta, ann, stable):
    import statsmodels.api as sm

    cand = ann[(ann["stage_fdr"] < 0.05) & ann["factor"].isin(stable)]["factor"]
    if cand.empty:
        return pd.DataFrame()
    age = pd.to_numeric(meta["age"], errors="coerce")
    male = meta["gender"].astype(str).str.upper().map({"MALE": 1.0, "FEMALE": 0.0})
    rows = []
    for f in cand:
        z = (Z[f] - Z[f].mean()) / Z[f].std()
        for label, X in [("unadjusted", pd.DataFrame({f: z})),
                         ("adjusted for age+sex+histology",
                          pd.DataFrame({f: z, "age": (age - age.mean()) / age.std(),
                                        "male": male,
                                        "lusc": (meta["histology"] == "LUSC").astype(float)}))]:
            d = pd.concat([X, y.rename("y")], axis=1).dropna()
            fit = sm.Logit(d["y"], sm.add_constant(d.drop(columns="y"))).fit(disp=0)
            ci = fit.conf_int().loc[f]
            rows.append({"factor": f, "model": label, "n": len(d),
                         "OR_per_SD": np.exp(fit.params[f]),
                         "CI_low": np.exp(ci[0]), "CI_high": np.exp(ci[1]),
                         "p": fit.pvalues[f]})
    return pd.DataFrame(rows)


# ---- main ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Compare stage-prediction models with repeated cross-validation."
        )
    )
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--perms", type=int, default=200,
                    help="label permutations for the low-dimensional models")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    meta, y, feats, ann, stable, Z = load()
    log("# Stage prediction report (step 5a)\n")
    log(f"- patients: {len(y)} ({int(y.sum())} late, {int((y == 0).sum())} early); "
        f"LUAD {int((meta.histology == 'LUAD').sum())}, "
        f"LUSC {int((meta.histology == 'LUSC').sum())}")
    log(f"- tissue source sites: {meta['tss'].nunique()}")
    log(f"- CV: {N_FOLDS}-fold x {args.repeats} repeats; identical splits "
        "for every model")
    log(f"- stable MOFA factors used: {', '.join(stable)}\n")

    schemes = {s: make_splits(meta, y, s, args.repeats) for s in ["patient", "site"]}
    rows, per_rep, oof, roc_data = [], [], pd.DataFrame(index=y.index), {}
    for name, X in feats.items():
        for scheme, splits in schemes.items():
            print(f"  {name:26s} {scheme:8s} ...", flush=True)
            aucs, baccs, by_hist, first = evaluate(X, y, meta, splits)
            rows.append({"model": name, "cv": scheme, "n_features": X.shape[1],
                         "auroc_mean": aucs.mean(), "auroc_sd": aucs.std(),
                         "bal_acc_mean": baccs.mean(),
                         "auroc_LUAD": np.mean(by_hist["LUAD"]),
                         "auroc_LUSC": np.mean(by_hist["LUSC"])})
            for i, a in enumerate(aucs):
                per_rep.append({"model": name, "cv": scheme, "repeat": i + 1,
                                "auroc": a})
            if scheme == "patient":
                oof[name] = first
                roc_data[name] = first

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "cv_results.tsv", sep="\t", index=False)
    pr = pd.DataFrame(per_rep)
    pr.to_csv(OUT / "cv_per_repeat.tsv", sep="\t", index=False)
    oof.assign(label_late=y, histology=meta["histology"]).to_csv(
        OUT / "oof_predictions.tsv", sep="\t")

    log("## Cross-validated performance\n")
    log("AUROC 0.5 = chance. Balanced accuracy 0.5 = chance (corrects for the "
        "2.7:1 class imbalance). Mean over repeats; sd = spread across repeats.\n")
    show = res.copy()
    for c in ["auroc_mean", "auroc_sd", "bal_acc_mean", "auroc_LUAD", "auroc_LUSC"]:
        show[c] = show[c].round(3)
    log(show.to_markdown(index=False) + "\n")

    # paired comparison against the clinical baseline
    log("## Paired comparison with the clinical baseline (patient CV)\n")
    log("Same splits for every model, so per-repeat differences are paired. "
        "`wins` = repeats in which the model beat the baseline.\n")
    base = pr.query("cv == 'patient' and model == 'clinical (age+sex)'")["auroc"].values
    comp = []
    for name in feats:
        if name == "clinical (age+sex)":
            continue
        a = pr.query("cv == 'patient' and model == @name")["auroc"].values
        d = a - base
        comp.append({"model": name, "delta_auroc_mean": round(d.mean(), 3),
                     "delta_sd": round(d.std(), 3),
                     "wins": f"{int((d > 0).sum())}/{len(d)}"})
    log(pd.DataFrame(comp).to_markdown(index=False) + "\n")

    # permutation tests for the low-dimensional models
    if args.perms > 0:
        log("## Label-permutation test (patient CV, repeat 1)\n")
        log(f"Labels shuffled {args.perms} times; p = share of shuffles that "
            "reach the observed AUROC. Low-dimensional models only (fast).\n")
        prow = []
        for name in ["clinical (age+sex)", "MOFA factors", "MOFA factors + clinical"]:
            obs = roc_auc_score(y, oof[name])
            p, null = permutation_p(feats[name], y, schemes["patient"][0],
                                    args.perms, obs, rng)
            prow.append({"model": name, "observed_auroc": round(obs, 3),
                         "null_mean": round(null.mean(), 3),
                         "null_95th": round(np.quantile(null, 0.95), 3),
                         "p": round(p, 4)})
        log(pd.DataFrame(prow).to_markdown(index=False) + "\n")

    adj = adjusted_factors(y, Z, meta, ann, stable)
    if not adj.empty:
        adj.to_csv(OUT / "adjusted_factors.tsv", sep="\t", index=False)
        log("## Stage-linked factors adjusted for age, sex and histology\n")
        log("Logistic regression of late stage on each factor (per SD). These "
            "factors were SELECTED for stage association in step 4 on the same "
            "patients, so p-values are optimistic; the question here is only "
            "whether each association survives adjustment.\n")
        a = adj.copy()
        for c in ["OR_per_SD", "CI_low", "CI_high"]:
            a[c] = a[c].round(2)
        a["p"] = a["p"].map(lambda v: f"{v:.2g}")
        log(a.to_markdown(index=False) + "\n")

    # ROC figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.5, 5))
    for name, p in roc_data.items():
        fpr, tpr, _ = roc_curve(y, p)
        ax.plot(fpr, tpr, lw=1.4, label=f"{name} ({roc_auc_score(y, p):.2f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Late vs early stage (patient CV, repeat 1)")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "roc.png", dpi=130)
    plt.close(fig)
    log("Figure: results/model/roc.png")

    (OUT / "model_report.md").write_text("\n".join(LOG), encoding="utf-8")
    print(f"\nReport: {OUT / 'model_report.md'}")


if __name__ == "__main__":
    main()
