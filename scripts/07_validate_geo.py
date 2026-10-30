#!/usr/bin/env python3
"""Validate transferable MOFA RNA-weight scores in an independent GEO cohort.

MOFA+ cannot be re-fitted on a new cohort, but each factor's RNA WEIGHTS can
be transferred: in the new cohort every gene is z-scored within histology
(mirroring MOFA's group centring) and a patient's factor score is the
weighted sum of gene z-scores. If a factor captures real biology, its score
should show the SAME direction of association with stage and survival as in
TCGA. Expected directions are written to the report BEFORE testing.

Workflow:
  1. python scripts/07_validate_geo.py --gse GSE37745 --inspect
       downloads the series matrix + platform annotation, prints all sample
       characteristic fields and the fields auto-detected for stage,
       histology, survival time and survival event.
  2. python scripts/07_validate_geo.py --gse GSE37745
       runs the validation. Override any auto-detected field if needed:
       --stage-field "tumor stage" --time-field "days to ..." --event-field "dead"

Outputs (``results/validation/<GSE>/``):
    sample_characteristics.tsv   parsed GEO annotations
    scores.tsv                   factor scores + clinical variables
    validation_report.md         hypotheses, TCGA self-check, results
    km_<factor>.png              Kaplan-Meier by factor-score median
    stage_<factor>.png           factor score by stage
Run from the repository root. Use ``--inspect`` to review parsed clinical
fields before running validation.
"""
import argparse
import datetime as dt
import gzip
import hashlib
import io
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROC = Path("data/processed")
MOFA = Path("results/mofa")
GEO_RAW = Path("data/raw/geo")
UA = {"User-Agent": "nsclc-multiomics-hackathon/1.0"}
LOG = []


def log(msg=""):
    print(msg)
    LOG.append(msg)


# ---- download --------------------------------------------------------------------
def http(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception as e:                       # noqa: BLE001
            if i == tries - 1:
                raise
            print(f"    retry ({e})")
            time.sleep(2 ** i)


def stub(acc):
    """GSE37745 -> GSE37nnn ; GPL570 -> GPLnnn ; GPL6884 -> GPL6nnn"""
    m = re.match(r"([A-Z]+)(\d+)", acc)
    pre, num = m.group(1), m.group(2)
    return f"{pre}{num[:-3]}nnn"


def fetch(url, dest, sources):
    if not dest.exists():
        print(f"  downloading {url}")
        raw = http(url)
        dest.write_bytes(raw)
    raw = dest.read_bytes()
    sources.append({"url": url, "file": str(dest), "bytes": len(raw),
                    "sha256": hashlib.sha256(raw).hexdigest(),
                    "accessed_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")})
    return raw


def read_series_matrix(gse, sources):
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub(gse)}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    raw = fetch(url, GEO_RAW / f"{gse}_series_matrix.txt.gz", sources)
    text = gzip.decompress(raw).decode("utf-8", errors="replace").splitlines()
    meta, table, in_table, platform = {}, [], False, None
    for line in text:
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            table.append(line)
            continue
        if line.startswith("!Series_platform_id"):
            platform = line.split("\t")[1].strip('"')
        if line.startswith("!Sample_"):
            parts = [p.strip('"') for p in line.split("\t")]
            meta.setdefault(parts[0], []).append(parts[1:])
    expr = pd.read_csv(io.StringIO("\n".join(table)), sep="\t", index_col=0,
                       na_values=["null", "NA", ""])
    expr.index = expr.index.astype(str).str.strip('"')
    expr.columns = [c.strip('"') for c in expr.columns]

    # characteristics: each "!Sample_characteristics_ch1" row holds "key: value"
    samples = expr.columns
    ch = pd.DataFrame(index=samples)
    for row in meta.get("!Sample_characteristics_ch1", []):
        for s, cell in zip(samples, row):
            if ":" in cell:
                k, v = cell.split(":", 1)
                ch.loc[s, k.strip().lower()] = v.strip()
    for key in ["!Sample_title", "!Sample_source_name_ch1"]:
        if key in meta:
            ch[key.replace("!Sample_", "").lower()] = meta[key][0]
    return expr, ch, platform


def platform_symbols(gpl, sources):
    """probe -> gene symbol from the GEO platform annotation."""
    candidates = [
        (f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub(gpl)}/{gpl}/annot/{gpl}.annot.gz",
         GEO_RAW / f"{gpl}.annot.gz"),
        (f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub(gpl)}/{gpl}/soft/{gpl}_family.soft.gz",
         GEO_RAW / f"{gpl}_family.soft.gz"),
    ]
    for url, dest in candidates:
        try:
            raw = fetch(url, dest, sources)
        except Exception as e:                       # noqa: BLE001
            print(f"  {url}: {e}")
            continue
        lines = gzip.decompress(raw).decode("utf-8", errors="replace").splitlines()
        begin = next((i for i, l in enumerate(lines)
                      if l.startswith("!platform_table_begin")), None)
        end = next((i for i, l in enumerate(lines)
                    if l.startswith("!platform_table_end")), len(lines))
        if begin is None:
            continue
        tab = pd.read_csv(io.StringIO("\n".join(lines[begin + 1:end])), sep="\t",
                          dtype=str, low_memory=False)
        col = next((c for c in tab.columns if re.fullmatch(
            r"(gene[ _]?symbol|symbol|ilmn_gene)", c.strip(), re.I)), None)
        if col is None:
            col = next((c for c in tab.columns if "symbol" in c.lower()), None)
        if col is None:
            continue
        sym = tab[col].fillna("").str.split(r"\s*///\s*|\s*//\s*|,").str[0].str.strip()
        out = pd.Series(sym.values, index=tab.iloc[:, 0].astype(str)).replace("", np.nan)
        print(f"  platform {gpl}: {out.notna().sum():,} probes with gene symbols "
              f"(column '{col}')")
        return out.dropna()
    sys.exit(f"Could not get gene symbols for platform {gpl}.")


# ---- characteristic parsing --------------------------------------------------------
def detect(ch, patterns, numeric=False):
    for p in patterns:
        for c in ch.columns:
            if re.search(p, c, re.I):
                if numeric and pd.to_numeric(ch[c], errors="coerce").notna().mean() < 0.5:
                    continue
                return c
    return None


ROMAN = {"IV": 4, "III": 3, "II": 2, "I": 1}


def parse_stage(v):
    if not isinstance(v, str):
        return np.nan
    s = v.upper().replace("STAGE", "").strip()
    m = re.match(r"^\s*([1-4])", s)
    if m:
        return int(m.group(1))
    for r in ["IV", "III", "II", "I"]:
        if s.startswith(r):
            return ROMAN[r]
    return np.nan


def parse_event(v):
    if not isinstance(v, str):
        return np.nan
    s = v.strip().lower()
    if s in {"1", "yes", "dead", "deceased", "died", "death", "true"} or s.startswith("dead"):
        return 1.0
    if s in {"0", "no", "alive", "living", "censored", "false"} or s.startswith("alive"):
        return 0.0
    return np.nan


def parse_hist(v):
    if not isinstance(v, str):
        return np.nan
    s = v.lower()
    if "adeno" in s:
        return "LUAD"
    if "squam" in s or s.strip() in {"scc", "sqcc"}:
        return "LUSC"
    return "other"


def clinical(ch, args):
    f = {
        "stage": args.stage_field or detect(ch, [r"^tumou?r stage", r"stage"]),
        "histology": args.hist_field or detect(ch, [r"histolog", r"histotype", r"cell type"]),
        "time": args.time_field or detect(ch, [r"death", r"surviv.*(day|month|time)",
                                               r"os.*(day|month|time)", r"follow"],
                                          numeric=True),
        "event": args.event_field or detect(ch, [r"^dead$", r"vital", r"status",
                                                 r"^death$", r"os.*(event|status)"]),
        "age": detect(ch, [r"^age"], numeric=True),
        "sex": detect(ch, [r"gender", r"^sex"]),
    }
    cl = pd.DataFrame(index=ch.index)
    cl["stage_num"] = ch[f["stage"]].map(parse_stage) if f["stage"] else np.nan
    cl["histology"] = ch[f["histology"]].map(parse_hist) if f["histology"] else "unknown"
    cl["os_time"] = pd.to_numeric(ch[f["time"]], errors="coerce") if f["time"] else np.nan
    cl["os_event"] = ch[f["event"]].map(parse_event) if f["event"] else np.nan
    cl["age"] = pd.to_numeric(ch[f["age"]], errors="coerce") if f["age"] else np.nan
    cl["label_late"] = cl["stage_num"].map({1: 0.0, 3: 1.0, 4: 1.0})
    return cl, f


# ---- scoring -------------------------------------------------------------------------
def gene_matrix(expr, symbols):
    X = expr.apply(pd.to_numeric, errors="coerce")
    if np.nanmax(X.values) > 50:          # not yet log-scaled
        X = np.log2(X.clip(lower=0) + 1)
    X = X.loc[X.index.intersection(symbols.index)]
    X["gene"] = symbols.reindex(X.index).values
    X["mean"] = X.drop(columns="gene").mean(axis=1)
    # one probe per gene: the most highly expressed (a common, stable choice)
    X = X.sort_values("mean").groupby("gene").last().drop(columns="mean")
    return X.T                             # samples x genes


def zscore_within(X, groups):
    Z = X.copy()
    for g in pd.Series(groups).dropna().unique():
        m = (groups == g)
        sub = X.loc[m]
        Z.loc[m] = (sub - sub.mean()) / sub.std().replace(0, np.nan)
    return Z.fillna(0)


def factor_scores(Z, W, factors):
    genes = Z.columns.intersection(W.index)
    S = pd.DataFrame(Z[genes].values @ W.loc[genes, factors].values,
                     index=Z.index, columns=factors)
    return (S - S.mean()) / S.std(), genes


# ---- statistics ----------------------------------------------------------------------
def auc_test(score, y):
    ok = y.notna() & score.notna()
    x1, x0 = score[ok & (y == 1)], score[ok & (y == 0)]
    if len(x1) < 5 or len(x0) < 5:
        return np.nan, np.nan, len(x1), len(x0)
    u, p = stats.mannwhitneyu(x1, x0, alternative="two-sided")
    return u / (len(x1) * len(x0)), p, len(x1), len(x0)


def cox(df, cols, strata=None):
    from lifelines import CoxPHFitter
    d = df[cols + ["os_time", "os_event"] + ([strata] if strata else [])].dropna()
    d = d[d["os_time"] > 0]
    if len(d) < 30 or d["os_event"].sum() < 10:
        return None
    if strata and d[strata].nunique() < 2:
        strata = None
    try:
        return CoxPHFitter().fit(d, "os_time", "os_event",
                                 strata=[strata] if strata else None), len(d), int(d["os_event"].sum())
    except Exception:                                # noqa: BLE001
        return None


def one_sided(p, observed_dir, expected_dir):
    if np.isnan(p):
        return np.nan
    return p / 2 if observed_dir == expected_dir else 1 - p / 2


# ---- main ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Validate MOFA factor scores in an independent GEO cohort."
        )
    )
    ap.add_argument("--gse", default="GSE37745")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--factors", nargs="*", default=None)
    ap.add_argument("--stage-field")
    ap.add_argument("--hist-field")
    ap.add_argument("--time-field")
    ap.add_argument("--event-field")
    args = ap.parse_args()

    GEO_RAW.mkdir(parents=True, exist_ok=True)
    out = Path("results/validation") / args.gse
    out.mkdir(parents=True, exist_ok=True)
    sources = []

    print(f"loading {args.gse} ...")
    expr, ch, gpl = read_series_matrix(args.gse, sources)
    ch.to_csv(out / "sample_characteristics.tsv", sep="\t")
    cl, fields = clinical(ch, args)

    if args.inspect:
        print(f"\n{args.gse}: {expr.shape[1]} samples, {expr.shape[0]:,} probes, platform {gpl}")
        print("\nSample characteristic fields (first 3 distinct values):")
        for c in ch.columns:
            vals = ch[c].dropna().unique()[:3]
            print(f"  - {c!r}: {', '.join(map(str, vals))}")
        print("\nAuto-detected fields:")
        for k, v in fields.items():
            print(f"  {k:9s} -> {v!r}")
        print("\nParsed summary:")
        print(f"  histology: {cl['histology'].value_counts(dropna=False).to_dict()}")
        print(f"  stage:     {cl['stage_num'].value_counts(dropna=False).sort_index().to_dict()}")
        print(f"  survival:  time for {int(cl['os_time'].notna().sum())}, "
              f"events {int(np.nansum(cl['os_event']))}")
        print("\nIf a field is wrong, rerun with e.g. --stage-field \"<name>\".")
        return

    # ---- factors and expected directions (from TCGA) ---------------------------------
    ann = pd.read_csv(MOFA / "factor_annotation.tsv", sep="\t")
    ann["stable"] = ann["stable"].astype(str) == "True"
    batch = (ann["site_fdr"] < 0.05) & (ann["site_excess"] > 0.05)
    factors = args.factors or ann.loc[ann["stable"] & (ann["stage_fdr"] < 0.05) & ~batch,
                                      "factor"].tolist()
    ann = ann.set_index("factor")
    W = pd.read_csv(MOFA / "weights_rna.tsv", sep="\t", index_col=0)

    log(f"# External validation report: {args.gse} (step 6)\n")
    log(f"- run date: {dt.date.today()}; platform {gpl}; {expr.shape[1]} samples")
    log("- fields used: " + ", ".join(f"{k} = '{v}'" for k, v in fields.items()))
    log("\n## Pre-specified hypotheses (directions taken from TCGA)\n")
    log("Expected directions were specified from TCGA before validation "
        "testing. One-sided p-values test the expected direction and are "
        "nominal.\n")
    expect = {}
    for f in factors:
        a = ann.loc[f]
        sd = "lower" if a["stage_auc"] < 0.5 else "higher"
        hd = "< 1" if a["os_hr_per_sd"] < 1 else "> 1"
        expect[f] = (np.sign(a["stage_auc"] - 0.5), np.sign(np.log(a["os_hr_per_sd"])))
        log(f"- **{f}**: score {sd} in late (III-IV) than early (I) stage "
            f"(TCGA AUC {a['stage_auc']:.3f}); survival HR per SD {hd} "
            f"(TCGA HR {a['os_hr_per_sd']:.2f}, FDR {a['os_fdr']:.2g})")

    # ---- TCGA self-check ----------------------------------------------------------------
    symbols = platform_symbols(gpl, sources)
    G = gene_matrix(expr, symbols)
    rna = pd.read_parquet(PROC / "rna.parquet")
    meta = pd.read_csv(PROC / "metadata.tsv", sep="\t", index_col=0, dtype={"tss": str})
    Ztcga = zscore_within(rna.loc[meta.index], meta["histology"].values)
    Zmofa = pd.read_csv(MOFA / "factors.tsv", sep="\t", index_col=0).loc[meta.index]
    overlap = G.columns.intersection(W.index)
    log("\n## Transfer check in TCGA\n")
    log(f"- genes with RNA weights: {len(W)}; measured in {args.gse}: "
        f"{len(overlap)} ({len(overlap) / len(W):.0%})")
    log("- correlation in TCGA between the MOFA factor and the transferred "
        "score (all genes / only genes measured in the validation cohort):\n")
    rows = []
    for f in factors:
        s_all, _ = factor_scores(Ztcga, W, [f])
        s_ov, _ = factor_scores(Ztcga[overlap], W.loc[overlap], [f])
        rows.append({"factor": f,
                     "r_all_genes": round(np.corrcoef(s_all[f], Zmofa[f])[0, 1], 3),
                     "r_overlap_genes": round(np.corrcoef(s_ov[f], Zmofa[f])[0, 1], 3)})
    log(pd.DataFrame(rows).to_markdown(index=False))
    log("\nValues near 1 mean the transferred score reproduces the factor; a "
        "low value means the validation result for that factor is not "
        "interpretable.\n")

    # ---- validation cohort ----------------------------------------------------------------
    keep = cl["histology"].isin(["LUAD", "LUSC"])
    cl = cl.loc[keep]
    G = G.loc[cl.index]
    Z = zscore_within(G, cl["histology"].values)
    S, genes = factor_scores(Z, W, factors)
    df = pd.concat([cl, S], axis=1)
    df.to_csv(out / "scores.tsv", sep="\t")

    log("## Validation cohort\n")
    log(f"- LUAD/LUSC samples used: {len(df)} "
        f"({df['histology'].value_counts().to_dict()})")
    log(f"- stage: {df['stage_num'].value_counts().sort_index().to_dict()} "
        f"-> early (I) {int((df.label_late == 0).sum())}, "
        f"late (III-IV) {int((df.label_late == 1).sum())}")
    log(f"- survival: {int(df['os_time'].notna().sum())} with follow-up, "
        f"{int(np.nansum(df['os_event']))} deaths\n")

    log("## Results\n")
    res = []
    for f in factors:
        auc, p, n1, n0 = auc_test(df[f], df["label_late"])
        exp_s, exp_h = expect[f]
        row = {"factor": f, "stage_auc": auc, "n_late": n1, "n_early": n0,
               "stage_p_two_sided": p,
               "stage_p_one_sided": one_sided(p, np.sign(auc - 0.5), exp_s),
               "stage_concordant": bool(np.sign(auc - 0.5) == exp_s) if pd.notna(auc) else None}
        for h in ["LUAD", "LUSC"]:
            row[f"stage_auc_{h}"] = auc_test(df.loc[df.histology == h, f],
                                             df.loc[df.histology == h, "label_late"])[0]
        c = cox(df, [f], strata="histology")
        if c:
            fit, n, ev = c
            hr = float(np.exp(fit.params_[f]))
            ph = float(fit.summary.loc[f, "p"])
            row.update({"os_hr_per_sd": hr, "os_p_two_sided": ph,
                        "os_p_one_sided": one_sided(ph, np.sign(np.log(hr)), exp_h),
                        "os_concordant": bool(np.sign(np.log(hr)) == exp_h),
                        "os_n": n, "os_events": ev})
            c2 = cox(df.assign(**{f"st{k}": (df.stage_num == k).astype(float)
                                  for k in [2, 3, 4]}).dropna(subset=["stage_num"]),
                     [f, "st2", "st3", "st4"], strata="histology")
            if c2:
                row["os_hr_adj_stage"] = float(np.exp(c2[0].params_[f]))
                row["os_p_adj_stage"] = float(c2[0].summary.loc[f, "p"])
        res.append(row)
    R = pd.DataFrame(res)
    R.to_csv(out / "validation_results.tsv", sep="\t", index=False)
    show = R.copy()
    for c in show.columns:
        if show[c].dtype.kind == "f":
            show[c] = show[c].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
    log(show.to_markdown(index=False) + "\n")
    log("Stage AUC compares late (III-IV) with early (I); stage II is "
        "excluded, as in TCGA. Survival HRs are per SD of the factor score "
        "from Cox models stratified by histology. `os_hr_adj_stage` and "
        "`os_p_adj_stage` are from a model additionally adjusted for stage; "
        "this is covariate adjustment, not multiple-testing correction. "
        "One-sided p-values are derived from the two-sided Cox p-value using "
        "the expected TCGA direction.\n")

    rep = R[(R["stage_concordant"] == True) & (R["stage_p_one_sided"] < 0.05)]["factor"].tolist()  # noqa: E712
    rep_os = R[(R.get("os_concordant") == True) & (R.get("os_p_one_sided") < 0.05)]["factor"].tolist() \
        if "os_concordant" in R else []
    log(f"- stage association meeting the nominal directional threshold "
        f"(one-sided p < 0.05): "
        f"{', '.join(rep) or 'none'}")
    log(f"- survival association meeting the nominal directional threshold "
        f"(one-sided p < 0.05): {', '.join(rep_os) or 'none'}")
    log("- Validation p-values across factors and outcomes (stage, "
        "unadjusted survival, and stage-adjusted survival) are not "
        "multiplicity-adjusted; no multiplicity-adjusted validation "
        "p-value is reported.")

    plots(df, factors, out)
    pd.DataFrame(sources).to_csv(out / "geo_sources.tsv", sep="\t", index=False)
    log(f"\nFigures: {out}/km_<factor>.png, {out}/stage_<factor>.png")
    log("\n## Limitations\n")
    log("- Only the RNA part of each factor can be transferred; methylation and "
        "CNV contributions are not tested here.")
    log("- Microarray vs RNA-seq: within-cohort z-scoring removes scale "
        "differences, but platform-specific probe behaviour remains.")
    log("- One probe per gene (highest mean expression); genes absent from the "
        "platform are dropped (see transfer check).")
    (out / "validation_report.md").write_text("\n".join(LOG), encoding="utf-8")
    print(f"\nReport: {out / 'validation_report.md'}")


def plots(df, factors, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from lifelines import KaplanMeierFitter
    except ImportError:
        KaplanMeierFitter = None
    for f in factors:
        if KaplanMeierFitter is not None and df["os_time"].notna().sum() > 30:
            fig, ax = plt.subplots(figsize=(5, 4))
            d = df.dropna(subset=["os_time", "os_event", f])
            hi = d[f] >= d[f].median()
            for lab, m in [("score above median", hi), ("score below median", ~hi)]:
                KaplanMeierFitter().fit(d.loc[m, "os_time"], d.loc[m, "os_event"],
                                        label=f"{lab} (n={int(m.sum())})").plot_survival_function(ax=ax, ci_show=False)
            ax.set_title(f"{f}: overall survival")
            ax.set_xlabel("time (units as in GEO)")
            fig.tight_layout()
            fig.savefig(out / f"km_{f}.png", dpi=130)
            plt.close(fig)
        fig, ax = plt.subplots(figsize=(5, 4))
        data = [df.loc[df.stage_num == s, f].dropna() for s in [1, 2, 3, 4]]
        ax.boxplot(data, showfliers=False)
        ax.set_xticks([1, 2, 3, 4], ["I", "II", "III", "IV"])
        ax.set_xlabel("stage")
        ax.set_ylabel(f"{f} score")
        ax.set_title(f"{f} score by stage")
        fig.tight_layout()
        fig.savefig(out / f"stage_{f}.png", dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    main()
