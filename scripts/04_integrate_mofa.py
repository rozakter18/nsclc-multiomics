#!/usr/bin/env python3
"""Fit grouped MOFA+ models and annotate their latent factors.

LUAD and LUSC are separate groups. MOFA+ is trained without stage labels;
associations with stage are descriptive, while predictive evaluation is in
``05_stage_model.py``. Outputs are written to ``results/mofa/``.

Run from the repository root.
"""
import argparse
import contextlib
import io
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

PROC = Path("data/processed")
OUT = Path("results/mofa")
VIEWS = ["rna", "methylation", "cnv"]
GROUPS = ["LUAD", "LUSC"]
N_PERM = 200
STABLE_R = 0.8     # min cross-seed |r| for a factor to count as reproducible
LOG = []


def log(msg=""):
    print(msg)
    LOG.append(msg)


# ---- data ----------------------------------------------------------------------
def load():
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0,
                       dtype={"tss": str})
    views = {v: pd.read_parquet(PROC / f"{v}.parquet").loc[meta.index]
             for v in VIEWS}
    return meta, views


# ---- training ------------------------------------------------------------------
def train(meta, views, factors, seed, path):
    from mofapy2.run.entry_point import entry_point

    data = []
    for v in VIEWS:
        data.append([views[v].loc[meta["histology"] == g].values.astype(np.float32)
                     for g in GROUPS])
    samples = [list(meta.index[meta["histology"] == g]) for g in GROUPS]

    buf = io.StringIO()                # MOFA prints a lot; capture it
    try:
        with contextlib.redirect_stdout(buf):
            ent = entry_point()
            # scale_views: each view gets unit total variance, so the view
            #   with the most features or largest values does not dominate.
            # center_groups: each histology is centred separately, so factors
            #   capture variation WITHIN histology, not LUAD-vs-LUSC.
            ent.set_data_options(scale_views=True, center_groups=True,
                                 use_float32=True)
            ent.set_data_matrix(
                data, likelihoods=["gaussian"] * len(VIEWS),
                views_names=VIEWS, groups_names=GROUPS,
                samples_names=samples,
                # RNA and CNV both use gene symbols (e.g. TP53), but MOFA+
                # needs names unique across views -> prefix with the view.
                features_names=[[f"{v}:{c}" for c in views[v].columns]
                                for v in VIEWS])
            ent.set_model_options(factors=factors, spikeslab_weights=True,
                                  ard_factors=True, ard_weights=True)
            # dropR2: factors explaining < 1% variance everywhere are pruned
            ent.set_train_options(iter=1000, convergence_mode="fast",
                                  seed=seed, dropR2=0.01, quiet=True)
            ent.build()
            ent.run()
            ent.save(str(path), save_data=False)
    except Exception:
        print(buf.getvalue())          # show MOFA's output if training fails
        raise


def read_model(path, meta):
    def s(x):
        return x.decode() if isinstance(x, bytes) else str(x)

    with h5py.File(path, "r") as f:
        Z = pd.concat([
            pd.DataFrame(f[f"expectations/Z/{g}"][()].T,
                         index=[s(x) for x in f[f"samples/{g}"][()]])
            for g in GROUPS])
        K = Z.shape[1]
        Z.columns = [f"Factor{k + 1}" for k in range(K)]
        W = {v: pd.DataFrame(f[f"expectations/W/{v}"][()].T,
                             index=[s(x).removeprefix(f"{v}:")
                                    for x in f[f"features/{v}"][()]],
                             columns=Z.columns) for v in VIEWS}
        r2 = []
        for g in GROUPS:
            arr = f[f"variance_explained/r2_per_factor/{g}"][()]   # views x K
            for i, v in enumerate(VIEWS):
                for k in range(K):
                    r2.append({"group": g, "view": v, "factor": Z.columns[k],
                               "r2": float(arr[i, k])})
        e = f["training_stats/elbo"][()]
        e = e[np.isfinite(e)]                 # ELBO is not computed every iter
        elbo = float(e[-1]) if len(e) else np.nan
    return Z.loc[meta.index], W, pd.DataFrame(r2), elbo


# ---- annotation helpers -------------------------------------------------------
def auc_mwu(x, y):
    """AUC of x for y=1 vs y=0 and the Mann-Whitney p-value."""
    x1, x0 = x[y == 1], x[y == 0]
    if len(x1) < 5 or len(x0) < 5:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(x1, x0, alternative="two-sided")
    return u / (len(x1) * len(x0)), p


def eta_squared(values, groups):
    df = pd.DataFrame({"v": values, "g": groups}).dropna()
    grand = df["v"].mean()
    between = df.groupby("g")["v"].apply(
        lambda x: len(x) * (x.mean() - grand) ** 2).sum()
    total = ((df["v"] - grand) ** 2).sum()
    return between / total if total else np.nan


def site_test(z, meta, rng):
    """Site effect beyond chance, within histology (as in the step 3 QC)."""
    hist = meta["histology"].values
    site = meta["tss"].astype(str).values
    zs = pd.Series(z)
    resid = (zs - zs.groupby(hist).transform("mean")).values
    obs = eta_squared(resid, site)
    null = []
    for _ in range(N_PERM):
        perm = site.copy()
        for h in GROUPS:
            idx = np.where(hist == h)[0]
            perm[idx] = rng.permutation(perm[idx])
        null.append(eta_squared(resid, perm))
    null = np.array(null)
    return obs - null.mean(), (1 + (null >= obs).sum()) / (1 + N_PERM)


def cox(z, meta):
    """Hazard ratio per SD of the factor, stratified by histology."""
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        return np.nan, np.nan
    df = pd.DataFrame({"z": (z - z.mean()) / z.std(),
                       "T": pd.to_numeric(meta["os_days"], errors="coerce"),
                       "E": pd.to_numeric(meta["os_event"], errors="coerce"),
                       "h": meta["histology"]}).dropna()
    df = df[df["T"] > 0]
    if len(df) < 50 or df["E"].sum() < 10:
        return np.nan, np.nan
    try:
        cph = CoxPHFitter().fit(df, "T", "E", strata=["h"])
        return float(np.exp(cph.params_["z"])), float(cph.summary.loc["z", "p"])
    except Exception:
        return np.nan, np.nan


def bh(p):
    """Benjamini-Hochberg FDR."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty_like(pv)
    res[order] = np.minimum(ranked, 1)
    out[ok] = res
    return out


def annotate(Z, r2, meta):
    rng = np.random.default_rng(0)
    y = meta["label_late"]
    sup = y.notna()
    rows = []
    for f in Z.columns:
        z = Z[f]
        row = {"factor": f}
        for g in GROUPS:
            for v in VIEWS:
                row[f"r2_{v}_{g}"] = round(float(r2.query(
                    "factor == @f and view == @v and group == @g")["r2"].iloc[0]), 2)
        auc, p = auc_mwu(z[sup].values, y[sup].values)
        row["stage_auc"], row["stage_p"] = auc, p
        for g in GROUPS:
            m = sup & (meta["histology"] == g)
            row[f"stage_auc_{g}"] = auc_mwu(z[m].values, y[m].values)[0]
        row["site_excess"], row["site_p"] = site_test(z.values, meta, rng)
        row["os_hr_per_sd"], row["os_p"] = cox(z, meta)
        age = pd.to_numeric(meta["age"], errors="coerce")
        ok = age.notna()
        row["age_rho"] = stats.spearmanr(z[ok], age[ok])[0] if ok.sum() > 20 else np.nan
        sex = meta["gender"].astype(str).str.upper().map({"MALE": 1, "FEMALE": 0})
        row["sex_auc"] = auc_mwu(z[sex.notna()].values, sex.dropna().values)[0]
        rows.append(row)
    ann = pd.DataFrame(rows)
    for col in ["stage_p", "site_p", "os_p"]:
        ann[col.replace("_p", "_fdr")] = bh(ann[col])
    return ann


def robustness(best_Z, others):
    """For each factor of the best model: best |r| with any factor of each
    other run. Values near 1 = the factor is reproducible across seeds."""
    out = pd.DataFrame(index=best_Z.columns)
    for name, Z2 in others.items():
        C = np.corrcoef(best_Z.values.T, Z2.values.T)[:best_Z.shape[1], best_Z.shape[1]:]
        out[f"max_abs_r_{name}"] = np.abs(C).max(axis=1).round(3)
    return out


# ---- figures -------------------------------------------------------------------
def plot_r2(r2, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(GROUPS), figsize=(5 * len(GROUPS), 4.5),
                             sharey=True)
    vmax = r2["r2"].max()
    for ax, g in zip(axes, GROUPS):
        m = r2[r2["group"] == g].pivot(index="factor", columns="view",
                                       values="r2")
        m = m.loc[sorted(m.index, key=lambda s: int(s[6:])), VIEWS]
        im = ax.imshow(m.values, cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(VIEWS)), VIEWS)
        ax.set_yticks(range(len(m)), m.index)
        ax.set_title(f"{g}: % variance explained")
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                ax.text(j, i, f"{m.values[i, j]:.1f}", ha="center",
                        va="center", fontsize=7,
                        color="white" if m.values[i, j] > vmax / 2 else "black")
    fig.colorbar(im, ax=axes, shrink=0.8)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_stage(Z, meta, ann, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top = ann.assign(d=(ann["stage_auc"] - 0.5).abs()).nlargest(4, "d")["factor"]
    sup = meta["label_late"].notna()
    fig, axes = plt.subplots(1, len(top), figsize=(3.5 * len(top), 3.8))
    for ax, f in zip(np.atleast_1d(axes), top):
        data, labels = [], []
        for g in GROUPS:
            for lab, name in [(0, "early"), (1, "late")]:
                m = sup & (meta["histology"] == g) & (meta["label_late"] == lab)
                data.append(Z.loc[m, f].values)
                labels.append(f"{g}\n{name}")
        ax.boxplot(data, showfliers=False)
        ax.set_xticks(range(1, len(labels) + 1), labels, fontsize=8)
        a = ann.set_index("factor").loc[f]
        ax.set_title(f"{f}  AUC={a['stage_auc']:.2f}  FDR={a['stage_fdr']:.2g}",
                     fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---- main ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Fit grouped MOFA+ models and annotate latent factors."
    )
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--factors", type=int, default=15)
    ap.add_argument("--reuse", action="store_true",
                    help="reuse trained models in results/mofa if present")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    meta, views = load()
    log("# MOFA+ integration report (step 4)\n")
    log(f"- patients: {len(meta)} "
        f"({', '.join(f'{g} {int((meta.histology == g).sum())}' for g in GROUPS)})")
    log("- views: " + ", ".join(f"{v} {views[v].shape[1]}" for v in VIEWS))
    log(f"- settings: up to {args.factors} factors, scale_views, groups "
        f"centred, factors < 1% R2 dropped, seeds {args.seeds}\n")

    runs = {}
    for seed in args.seeds:
        path = OUT / f"model_seed{seed}.hdf5"
        if not (args.reuse and path.exists()):
            print(f"training seed {seed} ...")
            train(meta, views, args.factors, seed, path)
        runs[seed] = read_model(path, meta)
    # highest final ELBO (the model's fit score); nan-safe
    best = max(runs, key=lambda s: -np.inf if np.isnan(runs[s][3]) else runs[s][3])
    Z, W, r2, _ = runs[best]

    log("## Training and robustness\n")
    log("| seed | factors | final ELBO |\n|---|---|---|")
    for s, (z, _, _, e) in runs.items():
        log(f"| {s}{' (best)' if s == best else ''} | {z.shape[1]} | {e:,.0f} |")
    log("")
    rob = robustness(Z, {f"seed{s}": runs[s][0] for s in runs if s != best})
    if not rob.empty:
        rob["stable"] = rob.min(axis=1) >= STABLE_R
        log("Best |correlation| of each factor with any factor of the other "
            f"runs (near 1 = reproducible). Factors below {STABLE_R} in any "
            "run are marked unstable and excluded from conclusions:\n")
        log(rob.to_markdown() + "\n")
        stable = rob["stable"]
    else:
        stable = pd.Series(True, index=Z.columns)

    # save tables
    Z.assign(histology=meta["histology"]).to_csv(OUT / "factors.tsv", sep="\t")
    for v in VIEWS:
        W[v].to_csv(OUT / f"weights_{v}.tsv", sep="\t")
    r2.to_csv(OUT / "variance_explained.tsv", sep="\t", index=False)

    log("## Variance explained (%)\n")
    tot = r2.groupby(["group", "view"])["r2"].sum().unstack()[VIEWS]
    log("Total per view and histology:\n")
    log(tot.round(1).to_markdown() + "\n")

    ann = annotate(Z, r2, meta)
    ann.insert(1, "stable", ann["factor"].map(stable).astype(bool))
    ann.to_csv(OUT / "factor_annotation.tsv", sep="\t", index=False)
    show = ann[["factor"] + [f"r2_{v}_{g}" for g in GROUPS for v in VIEWS]]
    log("Per factor (columns = view_histology):\n")
    log(show.to_markdown(index=False) + "\n")

    log("## Factor annotation\n")
    log("Stage AUC: late vs early (0.5 = no association; <0.5 means lower in "
        "late stage). The SIGN of a MOFA factor is arbitrary, so AUC 0.30 is "
        "as strong as 0.70 - read the distance from 0.5. Site: excess eta2 over a within-histology permutation "
        "null. OS: Cox hazard ratio per SD, stratified by histology. FDR = "
        "Benjamini-Hochberg across factors. All associations are descriptive "
        "(no labels were used to train MOFA+).\n")
    cols = ["factor", "stable", "stage_auc", "stage_auc_LUAD", "stage_auc_LUSC",
            "stage_fdr", "site_excess", "site_fdr", "os_hr_per_sd", "os_fdr",
            "age_rho", "sex_auc"]
    log(ann[cols].round(3).to_markdown(index=False) + "\n")

    def names(df):
        return ", ".join(df["factor"]) if len(df) else "none"

    batch = (ann["site_fdr"] < 0.05) & (ann["site_excess"] > 0.05)
    stage = ann["stage_fdr"] < 0.05
    log(f"- unstable across seeds (not interpreted): {names(ann[~ann['stable']])}")
    log(f"- factors with a site effect (FDR < 0.05, excess > 0.05): "
        f"{names(ann[batch])}")
    log(f"- stage-associated, stable, no strong site effect (FDR < 0.05): "
        f"**{names(ann[stage & ann['stable'] & ~batch])}**")
    log(f"- stage-associated but stable with a site effect (interpret with "
        f"caution): {names(ann[stage & ann['stable'] & batch])}")
    os_sig = ann[(ann["os_fdr"] < 0.05) & ann["stable"]]
    log(f"- associated with overall survival (FDR < 0.05, stable): "
        f"{names(os_sig)}")

    plot_r2(r2, OUT / "variance_explained.png")
    plot_stage(Z, meta, ann, OUT / "factors_vs_stage.png")
    log("\nFigures: results/mofa/variance_explained.png, "
        "results/mofa/factors_vs_stage.png")
    (OUT / "mofa_report.md").write_text("\n".join(LOG), encoding="utf-8")
    print(f"\nReport: {OUT / 'mofa_report.md'}")


if __name__ == "__main__":
    main()