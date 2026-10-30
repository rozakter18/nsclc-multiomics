#!/usr/bin/env python3
"""Preprocess the three omics layers and run cohort-level quality checks.

Inputs  (from step 2):  data/raw/*
Outputs:
    data/processed/rna.parquet           patients x genes   (log2 expression)
    data/processed/methylation.parquet   patients x probes  (M-values, NA kept)
    data/processed/cnv.parquet           patients x genes   (GISTIC2 -2..2)
    data/processed/metadata.tsv          patient labels and covariates
    data/processed/samples/*.csv         small representative extracts (tracked)
    results/qc/qc_report.md              feature counts, missingness, batch check
    results/qc/pca_<layer>.png           PCA coloured by histology/stage/site

Design decisions (see ``docs/objective.md``):
  * One primary tumour (sample type 01) per patient; patients must have all
    three layers and clinical data.
  * Features are ranked by the MEAN OF WITHIN-HISTOLOGY VARIANCES, so the
    selection is not dominated by LUAD-vs-LUSC differences.
  * No filter uses stage labels -> no label leakage. Supervised feature
    selection happens later, inside cross-validation folds.

Run from the repository root. ``--reuse`` can load selected existing
processed layers instead of recalculating them.
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path("data/raw")
PROC = Path("data/processed")
QC = Path("results/qc")
COHORTS = ["LUAD", "LUSC"]

# ---- parameters (documented in the QC report) --------------------------------
N_TOP = {"rna": 5000, "methylation": 5000, "cnv": 3000}
RNA_MIN_EXPR = 1.0        # log2(norm count + 1)
RNA_MIN_FRAC = 0.20       # gene must exceed RNA_MIN_EXPR in >= 20% of patients
METH_MAX_NA = 0.10        # drop probes missing in > 10% of patients
METH_CHUNK = 20000        # probes per chunk when streaming the methylation file
CNV_MIN_ALT_FRAC = 0.05   # gene altered (|value| >= 1) in >= 5% of patients
CNV_MAX_CORR = 0.9        # prune CNV features correlated above this
CNV_CANDIDATES = 8000     # how far down the variance ranking pruning scans
N_PERM = 200              # permutations for the site (batch) null distribution
PROBEMAP_GLOB = "probeMap_illuminaMethyl450*"   # .tsv or .gz from step 2

LOG = []  # lines collected for the QC report


class _Args:            # defaults when the module is imported, not run
    reuse = []


ARGS = _Args()


def log(msg=""):
    print(msg)
    LOG.append(msg)


# ---- helpers ------------------------------------------------------------------
def raw_file(cohort, stem):
    """Find a downloaded file by cohort prefix and dataset stem."""
    hits = sorted(RAW.glob(f"{cohort}_{stem}*"))
    if not hits:
        sys.exit(f"Missing raw file for {cohort} {stem}. Run step 2 first.")
    return hits[0]


def header(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return f.readline().rstrip("\n").split("\t")


def primary_tumour_columns(cols):
    """Map patient ID -> column name for primary tumours (type 01)."""
    out = {}
    for c in sorted(cols):
        if len(c) >= 15 and c.startswith("TCGA-") and c[13:15] == "01":
            out.setdefault(c[:12], c)  # keep the first sample per patient
    return out


def within_group_variance(df, groups):
    """Mean over groups of per-feature variance. df: patients x features."""
    variances = [df.loc[groups == group].var(axis=0) for group in COHORTS]
    return pd.concat(variances, axis=1).mean(axis=1)


def stage_group(s):
    if not isinstance(s, str):
        return "missing"
    s = s.upper().replace("STAGE", "").strip()
    if s.startswith("IV") or s.startswith("III"):
        return "late"
    if s.startswith("II"):
        return "intermediate"
    if s.startswith("I"):
        return "early"
    return "missing"


def overall_survival(m):
    """Overall survival (days, event) from a TCGA clinical matrix.
    Uses Xena's _OS/_OS_IND if present; otherwise derives it:
    time = days_to_death if deceased, else days to last follow-up."""
    def num(col):
        if col in m:
            return pd.to_numeric(m[col], errors="coerce")
        return pd.Series(np.nan, index=m.index)

    if "_OS" in m and "_OS_IND" in m and num("_OS").notna().any():
        return num("_OS"), num("_OS_IND")
    death = num("days_to_death")
    follow = num("days_to_last_followup")
    if follow.isna().all():
        follow = num("days_to_last_follow_up")
    vs = m["vital_status"].astype(str).str.upper() if "vital_status" in m \
        else pd.Series("", index=m.index)
    event = pd.Series(np.nan, index=m.index)
    event[vs.isin(["DECEASED", "DEAD"])] = 1.0
    event[vs.isin(["LIVING", "ALIVE"])] = 0.0
    time = death.where(event == 1, follow)
    time = time.where((time > 0) & event.notna())   # unknown status -> missing
    event = event.where(time.notna())
    return time, event


# ---- 1. cohort definition ------------------------------------------------------
def define_cohort():
    log("## 1. Cohort\n")
    layer_stems = {"rna": "HiSeqV2", "methylation": "HumanMethylation450",
                   "cnv": "Gistic2"}
    colmaps = {l: {} for l in layer_stems}   # layer -> patient -> (cohort, col)
    for c in COHORTS:
        for layer, stem in layer_stems.items():
            for pid, col in primary_tumour_columns(header(raw_file(c, stem))[1:]).items():
                colmaps[layer][pid] = (c, col)

    clin = []
    for c in COHORTS:
        df = pd.read_csv(raw_file(c, "clinicalMatrix"), sep="\t", low_memory=False)
        df = df[df["sampleID"].str[13:15] == "01"].copy()
        df["patient"] = df["sampleID"].str[:12]
        df["histology"] = c
        clin.append(df)
    clin = pd.concat(clin).drop_duplicates("patient").set_index("patient")

    patients = sorted(
        set(colmaps["rna"])
        & set(colmaps["methylation"])
        & set(colmaps["cnv"])
        & set(clin.index)
    )
    for layer in layer_stems:
        log(f"- {layer}: {len(colmaps[layer])} primary-tumour patients")
    log(f"- clinical: {len(clin)} patients")
    log(f"- **all three layers + clinical: {len(patients)} patients**\n")

    # metadata table
    m = clin.loc[patients]
    stage_col = ("pathologic_stage" if "pathologic_stage" in m else
                 next(c for c in m.columns if "stage" in c.lower()))
    meta = pd.DataFrame(index=pd.Index(patients, name="patient"))
    meta["histology"] = m["histology"]
    meta["stage_raw"] = m[stage_col]
    meta["stage_group"] = meta["stage_raw"].map(stage_group)
    meta["label_late"] = meta["stage_group"].map({"early": 0, "late": 1})
    meta["tss"] = [p[5:7] for p in patients]      # tissue source site (batch)
    for src, dst in [("gender", "gender"),
                     ("age_at_initial_pathologic_diagnosis", "age"),
                     ("tobacco_smoking_history", "smoking")]:
        meta[dst] = m[src] if src in m else np.nan
    meta["os_days"], meta["os_event"] = overall_survival(m)
    n_os = int(meta["os_days"].notna().sum())
    log(f"Overall survival available for {n_os} patients "
        f"({int(meta['os_event'].sum())} deaths)\n")

    counts = meta.groupby(["histology", "stage_group"]).size().unstack(fill_value=0)
    log(counts.to_markdown() + "\n")
    log(f"Supervised cohort (early vs late): "
        f"{int(meta['label_late'].notna().sum())} patients "
        f"({int((meta['label_late'] == 0).sum())} early, "
        f"{int((meta['label_late'] == 1).sum())} late)\n")
    return meta, colmaps


# ---- 2. RNA --------------------------------------------------------------------
def process_rna(meta, colmaps):
    log("## 2. RNA (HiSeqV2, log2 normalised counts)\n")
    parts = []
    for c in COHORTS:
        pids = meta.index[meta["histology"] == c]
        cols = [colmaps["rna"][p][1] for p in pids]
        path = raw_file(c, "HiSeqV2")
        df = pd.read_csv(path, sep="\t", index_col=0,
                         usecols=[header(path)[0]] + cols)
        df = df[cols].astype(np.float32)
        df.columns = list(pids)
        parts.append(df.T)
    X = pd.concat(parts, axis=0).loc[meta.index]
    X = X.loc[:, ~X.columns.duplicated()]
    log(f"- start: {X.shape[1]} genes")

    expressed = (X > RNA_MIN_EXPR).mean(axis=0) >= RNA_MIN_FRAC
    X = X.loc[:, expressed]
    log(
        f"- expressed (> {RNA_MIN_EXPR} in >= {RNA_MIN_FRAC:.0%}): "
        f"{X.shape[1]}"
    )

    v = within_group_variance(X, meta["histology"])
    X = X[v.nlargest(N_TOP["rna"]).index]
    log(f"- top {N_TOP['rna']} by within-histology variance: {X.shape[1]}\n")
    return X


# ---- 3. Methylation (two streaming passes) -------------------------------------
def load_probemap():
    hits = sorted(RAW.glob(PROBEMAP_GLOB))
    if not hits:
        log("- WARNING: probe annotation missing; chrX/Y probes NOT removed")
        return None
    pm = pd.read_csv(hits[0], sep="\t")
    pm = pm.rename(columns={pm.columns[0]: "probe"}).set_index("probe")
    return pm


def beta_to_m(b):
    b = np.clip(b, 0.001, 0.999)
    return np.log2(b / (1 - b))


def meth_chunks(cohort, pids, colmaps):
    path = raw_file(cohort, "HumanMethylation450")
    cols = [colmaps["methylation"][p][1] for p in pids]
    first = header(path)[0]
    reader = pd.read_csv(path, sep="\t", index_col=0, chunksize=METH_CHUNK,
                         usecols=[first] + cols,
                         dtype={c: np.float32 for c in cols})
    for chunk in reader:
        chunk = chunk[cols]
        chunk.columns = list(pids)
        yield beta_to_m(chunk.astype(np.float32))   # vectorised over the chunk


def process_methylation(meta, colmaps):
    log("## 3. DNA methylation (450K, beta -> M-values)\n")
    pm = load_probemap()

    # Pass 1: per-probe statistics within each cohort
    stats = {}
    for c in COHORTS:
        pids = list(meta.index[meta["histology"] == c])
        n_valid, s, ss = [], [], []
        for i, chunk in enumerate(meth_chunks(c, pids, colmaps)):
            n_valid.append(chunk.notna().sum(axis=1))
            s.append(chunk.sum(axis=1).astype(np.float64))
            ss.append((chunk.astype(np.float64) ** 2).sum(axis=1))
            print(f"\r    pass 1 {c}: {(i + 1) * METH_CHUNK:,} probes", end="")
        print()
        n = pd.concat(n_valid)
        mean = pd.concat(s) / n
        var = pd.concat(ss) / n - mean ** 2
        stats[c] = pd.DataFrame({"n": n, "var": var, "N": len(pids)})

    common = stats["LUAD"].index.intersection(stats["LUSC"].index)
    log(f"- start: {len(common):,} probes in both cohorts")
    na_frac = 1 - (stats["LUAD"].loc[common, "n"] + stats["LUSC"].loc[common, "n"]) \
        / (stats["LUAD"]["N"].iloc[0] + stats["LUSC"]["N"].iloc[0])
    keep = na_frac[na_frac <= METH_MAX_NA].index
    log(f"- missing <= {METH_MAX_NA:.0%}: {len(keep):,}")

    if pm is not None and "chrom" in pm:
        sex = pm.index[pm["chrom"].isin(["chrX", "chrY", "X", "Y"])]
        keep = keep.difference(sex)
        log(f"- autosomal (chrX/Y removed): {len(keep):,}")

    v = (stats["LUAD"].loc[keep, "var"] + stats["LUSC"].loc[keep, "var"]) / 2
    top = v.nlargest(N_TOP["methylation"]).index
    log(f"- top {N_TOP['methylation']} by within-histology variance: {len(top)}")

    # Pass 2: keep only the selected probes
    parts = []
    topset = set(top)
    for c in COHORTS:
        pids = list(meta.index[meta["histology"] == c])
        kept = [ch.loc[ch.index.intersection(topset)]
                for ch in meth_chunks(c, pids, colmaps)]
        parts.append(pd.concat(kept).T)
        print(f"    pass 2 {c}: done")
    X = pd.concat(parts, axis=0).loc[meta.index, list(top)]
    log(f"- remaining missing values: {X.isna().mean().mean():.2%} "
        "(kept as NA; MOFA+ handles them, ML imputes inside CV folds)\n")
    return X


# ---- 4. CNV --------------------------------------------------------------------
def process_cnv(meta, colmaps):
    log("## 4. Copy number (GISTIC2 thresholded, -2..2)\n")
    parts = []
    for c in COHORTS:
        pids = meta.index[meta["histology"] == c]
        cols = [colmaps["cnv"][p][1] for p in pids]
        path = raw_file(c, "Gistic2")
        first = header(path)[0]
        df = pd.read_csv(path, sep="\t", index_col=0, usecols=[first] + cols)
        df = df[cols].astype(np.float32)
        df.columns = list(pids)
        parts.append(df.T)
    X = pd.concat(parts, axis=0).loc[meta.index]
    X = X.loc[:, ~X.columns.duplicated()].dropna(axis=1)
    log(f"- start: {X.shape[1]} genes")

    altered = (X.abs() >= 1).mean(axis=0) >= CNV_MIN_ALT_FRAC
    X = X.loc[:, altered]
    log(f"- altered in >= {CNV_MIN_ALT_FRAC:.0%} of patients: {X.shape[1]}")

    # Neighbouring genes on the same DNA segment get IDENTICAL thresholded
    # values. Collapse identical profiles into one feature so CNV is not
    # over-counted in integration; keep the gene membership of each feature.
    key = pd.util.hash_pandas_object(X.T, index=False)
    groups = pd.Series(X.columns, index=X.columns).groupby(key.values).agg(list)
    rep = {g[0]: g for g in groups}                 # representative -> members
    X = X[list(rep)]
    log(f"- unique copy-number profiles (identical genes collapsed): {X.shape[1]}")

    # Neighbouring genes are often NEARLY identical (breakpoints differ a
    # little between patients), so also prune by correlation: walk down the
    # variance ranking and keep a feature only if |r| <= CNV_MAX_CORR with
    # every feature already kept. Dropped genes join the kept feature they
    # correlate with most, so no gene is lost from the membership table.
    v = within_group_variance(X, meta["histology"])
    cand = list(v.nlargest(min(len(v), CNV_CANDIDATES)).index)
    Xc = X[cand] - X[cand].groupby(meta["histology"].values).transform("mean")
    Zs = ((Xc - Xc.mean()) / Xc.std().replace(0, 1)).values.astype(np.float32)
    n = Zs.shape[0]
    kept = []                                   # column indices into cand
    for i in range(len(cand)):
        if kept:
            r = np.abs(Zs[:, i] @ Zs[:, kept]) / (n - 1)
            j = int(r.argmax())
            if r[j] > CNV_MAX_CORR:
                rep[cand[kept[j]]] += rep.pop(cand[i])
                continue
        kept.append(i)
        if len(kept) == N_TOP["cnv"]:
            break
    X = X[[cand[i] for i in kept]]
    log(f"- correlation-pruned (|r| <= {CNV_MAX_CORR}, scanning the top "
        f"{len(cand)} by within-histology variance): {X.shape[1]} features")
    members = pd.DataFrame({
        "feature": X.columns,
        "n_genes": [len(rep[f]) for f in X.columns],
        "genes": [",".join(rep[f]) for f in X.columns],
    })
    members.to_csv(PROC / "cnv_feature_groups.tsv", sep="\t", index=False)
    log(f"- each CNV feature stands for a median of "
        f"{int(members['n_genes'].median())} genes "
        f"(max {members['n_genes'].max()}); membership -> "
        "data/processed/cnv_feature_groups.tsv\n")
    return X


# ---- 5. QC: PCA and batch check -------------------------------------------------
def eta_squared(values, groups):
    """Fraction of variance in `values` explained by categorical `groups`."""
    df = pd.DataFrame({"v": values, "g": groups}).dropna()
    grand = df["v"].mean()
    between = df.groupby("g")["v"].apply(lambda x: len(x) * (x.mean() - grand) ** 2).sum()
    total = ((df["v"] - grand) ** 2).sum()
    return between / total if total else np.nan


def qc_layer(name, X, meta):
    from sklearn.decomposition import PCA
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Z = X.fillna(X.mean())
    Z = (Z - Z.mean()) / Z.std().replace(0, 1)
    pca = PCA(n_components=5, random_state=0)
    pcs = pca.fit_transform(Z.values)

    # Site (batch) effect, done properly:
    #  1. remove the histology effect from each PC (subtract histology means),
    #     because sites contribute mostly one histology;
    #  2. measure how much of the remainder site explains (eta2);
    #  3. compare with a permutation null: shuffle site labels WITHIN each
    #     histology. Many small groups explain some variance by chance alone,
    #     so only the EXCESS over the null is evidence of a batch effect.
    hist = meta["histology"].values
    site = meta["tss"].astype(str).values
    stage = meta["stage_group"].replace("missing", np.nan).values
    rng = np.random.default_rng(0)
    rows = []
    for i in range(5):
        pc = pd.Series(pcs[:, i])
        resid = (pc - pc.groupby(hist).transform("mean")).values
        obs = eta_squared(resid, site)
        null = []
        for _ in range(N_PERM):
            perm = site.copy()
            for h in COHORTS:
                idx = np.where(hist == h)[0]
                perm[idx] = rng.permutation(perm[idx])
            null.append(eta_squared(resid, perm))
        null = np.array(null)
        rows.append({
            "PC": f"PC{i + 1}",
            "var_explained": f"{pca.explained_variance_ratio_[i]:.1%}",
            "eta2_histology": round(eta_squared(pcs[:, i], hist), 3),
            "eta2_stage": round(eta_squared(pcs[:, i], stage), 3),
            "eta2_site_resid": round(obs, 3),
            "site_null_mean": round(null.mean(), 3),
            "site_excess": round(obs - null.mean(), 3),
            "site_p": f"{(1 + (null >= obs).sum()) / (1 + N_PERM):.3f}",
        })
    tab = pd.DataFrame(rows)
    log(f"### {name}\n")
    log(tab.to_markdown(index=False) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    top_sites = meta["tss"].value_counts().index[:8]
    site = meta["tss"].where(meta["tss"].isin(top_sites), "other")
    for ax, (title, lab) in zip(axes, [("Histology", meta["histology"]),
                                       ("Stage", meta["stage_group"]),
                                       ("Tissue source site (top 8)", site)]):
        for g in sorted(lab.unique()):
            mask = (lab == g).values
            ax.scatter(pcs[mask, 0], pcs[mask, 1], s=8, alpha=0.6, label=g)
        ax.set_title(f"{name}: {title}")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
        ax.legend(fontsize=7, markerscale=2, frameon=False)
    fig.tight_layout()
    fig.savefig(QC / f"pca_{name}.png", dpi=130)
    plt.close(fig)


def sample_qc(layers, meta):
    log("### Per-patient checks\n")
    flags = pd.DataFrame(index=meta.index)
    flags["meth_na_frac"] = layers["methylation"].isna().mean(axis=1)
    flags["rna_median"] = layers["rna"].median(axis=1)
    z = (flags["rna_median"] - flags["rna_median"].mean()) / flags["rna_median"].std()
    out = flags[(flags["meth_na_frac"] > 0.2) | (z.abs() > 4)]
    log(f"- patients with > 20% missing methylation or extreme RNA median "
        f"(|z| > 4): {len(out)}")
    if len(out):
        log(out.round(3).to_markdown() + "\n")
    log("- flagged patients are reported, not removed; decide after review\n")


# ---- main ----------------------------------------------------------------------
def main():
    PROC.mkdir(parents=True, exist_ok=True)
    (PROC / "samples").mkdir(exist_ok=True)
    QC.mkdir(parents=True, exist_ok=True)

    log("# QC report - preprocessing (step 3)\n")
    meta, colmaps = define_cohort()
    steps = {"rna": process_rna, "methylation": process_methylation,
             "cnv": process_cnv}
    layers = {}
    for name, fn in steps.items():
        if name in ARGS.reuse and (PROC / f"{name}.parquet").exists():
            layers[name] = pd.read_parquet(PROC / f"{name}.parquet").loc[meta.index]
            log(f"## {name}: reused data/processed/{name}.parquet "
                f"({layers[name].shape[1]} features; --reuse)\n")
        else:
            layers[name] = fn(meta, colmaps)

    log("## 5. QC\n")
    log("PCA on each processed layer (standardised features, NA -> feature "
        "mean for PCA only). eta2 = share of a PC's variance explained by a "
        "grouping.\n")
    log("Batch check: `eta2_site_resid` is the share of a PC's variance "
        "explained by tissue source site AFTER removing histology means. "
        f"`site_null_mean` is the same quantity with site labels shuffled "
        f"within histology ({N_PERM} permutations); many small sites explain "
        "some variance by chance. `site_excess` = observed - null is the "
        "batch effect beyond chance; `site_p` is the permutation p-value.\n")
    for name, X in layers.items():
        qc_layer(name, X, meta)
    sample_qc(layers, meta)

    log("## 6. Outputs\n")
    meta.to_csv(PROC / "metadata.tsv", sep="\t")
    for name, X in layers.items():
        X.to_parquet(PROC / f"{name}.parquet")
        # small tracked extract for the submission ("representative sample")
        X.iloc[:50, :200].round(4).to_csv(PROC / "samples" / f"{name}_sample.csv")
        log(f"- {name}: {X.shape[0]} patients x {X.shape[1]} features "
            f"-> data/processed/{name}.parquet")
    log("- metadata -> data/processed/metadata.tsv")

    log("\n## Parameters\n")
    for k, v in [("N_TOP", N_TOP), ("RNA_MIN_EXPR", RNA_MIN_EXPR),
                 ("RNA_MIN_FRAC", RNA_MIN_FRAC), ("METH_MAX_NA", METH_MAX_NA),
                 ("CNV_MIN_ALT_FRAC", CNV_MIN_ALT_FRAC),
                 ("CNV_MAX_CORR", CNV_MAX_CORR),
                 ("CNV_CANDIDATES", CNV_CANDIDATES), ("N_PERM", N_PERM)]:
        log(f"- {k} = {v}")

    (QC / "qc_report.md").write_text("\n".join(LOG), encoding="utf-8")
    print(f"\nQC report: {QC / 'qc_report.md'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Preprocess the TCGA omics layers and run QC."
    )
    ap.add_argument("--reuse", nargs="*", default=[],
                    choices=["rna", "methylation", "cnv"],
                    help="load these layers from existing data/processed "
                         "files instead of reprocessing (e.g. --reuse "
                         "methylation to skip the slow two-pass step)")
    ARGS = ap.parse_args()
    main()