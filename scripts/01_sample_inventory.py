#!/usr/bin/env python3
"""Inventory TCGA LUAD/LUSC samples available across the Xena omics layers.

Writes sample-level layer presence, stage counts by layer set, and a summary
under ``results/inventory/``. Run from the repository root.
"""
import io
import sys
import urllib.request
from pathlib import Path

import pandas as pd
import xenaPython as xena

HUB = "https://tcga.xenahubs.net"
COHORTS = {"LUAD": "TCGA Lung Adenocarcinoma (LUAD)",
           "LUSC": "TCGA Lung Squamous Cell Carcinoma (LUSC)"}

# Keyword rules for picking one dataset per layer from each cohort's list.
# The first dataset whose name contains ALL keywords (case-insensitive) wins.
LAYER_RULES = {
    "rna":         [["/hiseqv2$"], ["hiseqv2"]],
    "methylation": [["humanmethylation450"]],
    "cnv":         [["gistic2", "thresholded"], ["gistic2"]],
    "rppa":        [["rppa_rbn"], ["rppa"]],
}

OUT = Path("results/inventory")

EXCLUDE = ["percentile", "pancan", "exon", "_paired"]


def _match(key, name):
    # A key ending in "$" must match the end of the name exactly.
    if key.endswith("$"):
        return name.endswith(key[:-1])
    return key in name


def pick_dataset(names, rules):
    lower = {n: n.lower() for n in names
             if not any(x in n.lower() for x in EXCLUDE)}
    for keys in rules:
        for n, ln in lower.items():
            if all(_match(k, ln) for k in keys):
                return n
    return None


def clinical_matrix(cohort_code):
    url = f"{HUB}/download/TCGA.{cohort_code}.sampleMap%2F{cohort_code}_clinicalMatrix"
    with urllib.request.urlopen(url, timeout=120) as r:
        return pd.read_csv(io.BytesIO(r.read()), sep="\t", low_memory=False)


def stage_group(s):
    if not isinstance(s, str):
        return "missing"
    s = s.upper().replace("STAGE", "").strip()
    if s.startswith("IV"):
        return "late"
    if s.startswith("III"):
        return "late"
    if s.startswith("II"):
        return "intermediate"
    if s.startswith("I") or s == "0":
        return "early"
    return "missing"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    presence, clin_rows, chosen = [], [], []

    for code, cohort in COHORTS.items():
        names = [d["name"] for d in xena.dataset_list(HUB, [cohort])]
        for layer, rules in LAYER_RULES.items():
            ds = pick_dataset(names, rules)
            chosen.append({"cohort": code, "layer": layer, "dataset": ds})
            if ds is None:
                print(f"[warn] {code}: no dataset matched layer '{layer}'",
                      file=sys.stderr)
                continue
            samples = xena.dataset_samples(HUB, ds, None)
            presence += [{"sample": s, "cohort": code, "layer": layer}
                         for s in samples]
            print(f"{code:4s} {layer:11s} {len(samples):5d}  {ds}")

        cm = clinical_matrix(code)
        stage_col = next((c for c in cm.columns
                          if c.lower() == "pathologic_stage"), None)
        if stage_col is None:
            stage_col = next(c for c in cm.columns if "stage" in c.lower())
        clin_rows.append(pd.DataFrame({
            "sample": cm["sampleID"],
            "cohort": code,
            "stage_raw": cm[stage_col],
        }))

    pres = (pd.DataFrame(presence)
            .assign(v=1)
            .pivot_table(index=["sample", "cohort"], columns="layer",
                         values="v", fill_value=0)
            .reset_index())
    clin = pd.concat(clin_rows).drop_duplicates("sample")
    df = pres.merge(clin.drop(columns="cohort"), on="sample", how="left")

    df["sample_type"] = df["sample"].str[13:15]
    df["patient"] = df["sample"].str[:12]
    df["stage_group"] = df["stage_raw"].map(stage_group)
    df.to_csv(OUT / "layer_samples.tsv", sep="\t", index=False)
    pd.DataFrame(chosen).to_csv(OUT / "datasets_used.tsv", sep="\t", index=False)

    # Primary tumours only, one sample per patient
    tum = (df[df["sample_type"] == "01"]
           .sort_values("sample")
           .drop_duplicates("patient"))

    layer_sets = {
        "RNA only":                    ["rna"],
        "RNA+meth":                    ["rna", "methylation"],
        "RNA+meth+CNV":                ["rna", "methylation", "cnv"],
        "RNA+meth+CNV+RPPA":           ["rna", "methylation", "cnv", "rppa"],
    }
    rows = []
    for label, layers in layer_sets.items():
        have = [l for l in layers if l in tum.columns]
        if len(have) < len(layers):
            continue
        sub = tum[(tum[layers] == 1).all(axis=1)]
        counts = sub.groupby(["cohort", "stage_group"]).size().unstack(fill_value=0)
        for cohort, r in counts.iterrows():
            rows.append({"layers": label, "cohort": cohort, **r.to_dict(),
                         "total": int(r.sum())})
        tot = sub["stage_group"].value_counts().to_dict()
        rows.append({"layers": label, "cohort": "LUAD+LUSC", **tot,
                     "total": len(sub)})
    summary = pd.DataFrame(rows).fillna(0)
    for c in ["early", "intermediate", "late", "missing", "total"]:
        if c in summary:
            summary[c] = summary[c].astype(int)
    summary.to_csv(OUT / "stage_by_layers.tsv", sep="\t", index=False)

    normals = df[df["sample_type"] == "11"].groupby("cohort").size()
    with open(OUT / "inventory_report.md", "w") as f:
        f.write("# Sample inventory — TCGA LUAD/LUSC (Xena TCGA hub)\n\n")
        f.write("Primary tumours (type 01), one per patient.\n\n")
        f.write(summary.to_markdown(index=False))
        f.write("\n\nSolid-tissue normals (type 11): "
                + ", ".join(f"{k} {v}" for k, v in normals.items()) + "\n")
    print("\n" + summary.to_string(index=False))
    print(f"\nWritten to {OUT}/")


if __name__ == "__main__":
    main()