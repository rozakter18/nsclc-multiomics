#!/usr/bin/env python3
"""Mine biological evidence for stage-linked MOFA factors.

For each selected MOFA factor:
  1. Gene-level multi-omics support
       RNA weights per gene; methylation probes -> genes (probe annotation);
       CNV features -> all member genes. Per view, |weight| is scaled to
       [0, 1]; the DATA-SUPPORT score is the sum across views, so genes
       supported in several layers rank higher.
  2. Pathway enrichment (over-representation, hypergeometric test) of the
       top genes against Reactome and MSigDB Hallmark, using the MEASURED
       genes as background (not the whole genome), BH-FDR per library.
  3. Literature evidence (PubMed E-utilities): papers mentioning the gene
       with lung cancer, and with lung cancer progression/stage/prognosis.
  4. Curated evidence: Open Targets NSCLC association score; STRING
       interactions with established NSCLC genes.
  5. Evidence tags (known / plausible / limited evidence) and a GAP score:
       percentile(data support) - percentile(progression literature).
       High gap = strong multi-omics support, little prior study.
  6. A small knowledge graph (factor - gene - pathway - NSCLC) as TSV and
       GraphML.

All network results are cached in results/knowledge/cache/, so re-runs are
fast and reproducible; delete the cache to refresh. Every external resource
is logged with URL, access date and checksum.

Inputs:  results/mofa/{weights_*.tsv, factor_annotation.tsv}
         data/processed/cnv_feature_groups.tsv
         data/raw/probeMap_illuminaMethyl450*  (probe -> gene)
Outputs: results/knowledge/ (see report)

Run from the repository root. External resources are cached under
``results/knowledge/cache/`` and provenance is written to the output report.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import math
import os
import re
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

MOFA = Path("results/mofa")
PROC = Path("data/processed")
RAW = Path("data/raw")
OUT = Path("results/knowledge")
CACHE = OUT / "cache"
VIEWS = ["rna", "methylation", "cnv"]
UA = {"User-Agent": "nsclc-multiomics-hackathon/1.0 (research; python urllib)"}
LOG = []

# ---- parameters (reported) -----------------------------------------------------
TOP_ORA = 200            # top genes per factor used for enrichment
LAYER_SUPPORT = 0.2      # scaled |w| >= this counts as support in a layer
GS_MIN, GS_MAX = 10, 500 # gene-set size limits (within the measured universe)
KNOWN_PROG_PAPERS = 20   # >= this many progression papers -> "known"
KNOWN_OT_SCORE = 0.5     # Open Targets NSCLC score >= this -> "known"
PLAUSIBLE_OT_SCORE = 0.1
N_KNOWN_DRIVERS = 50     # top Open Targets NSCLC genes used as network anchors
STRING_MIN_SCORE = 700   # STRING combined score (0-1000), "high confidence"

LUNG = ('("lung neoplasms"[mh] OR NSCLC[tiab] OR "lung cancer"[tiab] OR '
        '"lung adenocarcinoma"[tiab] OR "lung squamous"[tiab] OR '
        '"non-small cell lung"[tiab])')
PROG = ('(stage[tiab] OR staging[tiab] OR progression[tiab] OR metasta*[tiab] '
        'OR invasion[tiab] OR prognos*[tiab] OR survival[tiab])')

# Gene symbols that are also common words/abbreviations give noisy PubMed counts
AMBIGUOUS = {"SET", "MET", "CAT", "ACE", "ARC", "CAP", "CAR", "MAX", "REST",
             "SHE", "WAS", "IMPACT", "CLOCK", "SPP1", "PC", "PI", "GC", "CD",
             "HR", "AR", "ER", "PR", "TF", "IMP", "MAL", "LARGE", "NOV", "BAD",
             "GAS", "AGO", "PIGS", "FAT", "MARCH", "SEPT", "DEC", "APEX",
             "CASP", "SON", "TANK", "WARS", "CARS", "MAP", "STAR", "BID",
             "ATM", "RAN", "PAH", "COPE", "MICE"}


def log(msg=""):
    print(msg)
    LOG.append(msg)


# ---- HTTP with retries and provenance -------------------------------------------
SOURCES = []   # provenance records


def http(url, data=None, headers=None, timeout=120, tries=4):
    h = dict(UA, **(headers or {}))
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=h)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:                       # noqa: BLE001
            if i == tries - 1:
                raise
            wait = 2 ** i
            print(f"    retry in {wait}s ({type(e).__name__}: {e})")
            time.sleep(wait)


def record(name, url, payload=None, note=""):
    SOURCES.append({
        "resource": name, "url": url,
        "accessed_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        "sha256": hashlib.sha256(payload).hexdigest() if payload else "",
        "note": note})


def cached_json(name, fn):
    """Return cached JSON result, or compute with fn() and cache it."""
    path = CACHE / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    val = fn()
    empty = not val or (isinstance(val, dict) and not val.get("scores", True))
    if not empty:                     # never cache a failed fetch
        path.write_text(json.dumps(val), encoding="utf-8")
    return val


# ---- 1. gene-level multi-omics support ------------------------------------------
def probe_to_genes():
    """probe -> genes, and gene -> (chrom, approx. position) from probe
    coordinates (used for the positional-cluster check in enrichment)."""
    hits = sorted(RAW.glob("probeMap_illuminaMethyl450*"))
    if not hits:
        log("- WARNING: probe annotation missing; methylation not mapped to genes")
        return {}, {}
    pm = pd.read_csv(hits[0], sep="\t", dtype=str)
    pm = pm.rename(columns={pm.columns[0]: "probe"})
    out, pos = {}, {}
    start = pd.to_numeric(pm.get("chromStart"), errors="coerce") \
        if "chromStart" in pm else pd.Series(np.nan, index=pm.index)
    chrom = pm["chrom"] if "chrom" in pm else pd.Series("", index=pm.index)
    from collections import Counter
    for p, g, c, s in zip(pm["probe"], pm["gene"].fillna("."), chrom, start):
        genes = [x for x in re.split(r"[,;]", g) if x and x != "."]
        if genes:
            out[p] = genes
            # positions only from probes with a real chromosome and coordinate
            if pd.notna(s) and isinstance(c, str) and c.strip() not in ("", "."):
                for x in genes:
                    pos.setdefault(x, []).append((c, s))
    gene_pos = {}
    for g, v in pos.items():
        c = Counter(ch for ch, _ in v).most_common(1)[0][0]
        gene_pos[g] = (c, float(np.median([s for ch, s in v if ch == c])))
    return out, gene_pos


def gene_support(factor, W, p2g, cnv_groups):
    """DataFrame indexed by gene: scaled |w| per view, signed weights, score."""
    per_view = {}
    # RNA
    w = W["rna"][factor]
    per_view["rna"] = pd.DataFrame({"w": w, "absw": w.abs()})
    # methylation: probe -> genes, keep the probe with the largest |w| per gene
    w = W["methylation"][factor]
    rows = [(g, val) for probe, val in w.items() for g in p2g.get(probe, [])]
    if rows:
        df = pd.DataFrame(rows, columns=["gene", "w"])
        df["absw"] = df["w"].abs()
        per_view["methylation"] = df.sort_values("absw").groupby("gene").last()
    # CNV: feature -> member genes, each member inherits the feature weight
    w = W["cnv"][factor]
    rows = [(g, w[f], len(cnv_groups.get(f, [f]))) for f in w.index
            for g in cnv_groups.get(f, [f])]
    df = pd.DataFrame(rows, columns=["gene", "w", "region_genes"])
    df["absw"] = df["w"].abs()
    per_view["cnv"] = df.sort_values("absw").groupby("gene").last()

    genes = sorted(set().union(*[set(d.index) for d in per_view.values()]))
    out = pd.DataFrame(index=pd.Index(genes, name="gene"))
    for v, d in per_view.items():
        scaled = d["absw"] / d["absw"].max() if d["absw"].max() > 0 else d["absw"]
        out[f"s_{v}"] = scaled.reindex(genes).fillna(0.0)
        out[f"w_{v}"] = d["w"].reindex(genes)
    for v in VIEWS:
        if f"s_{v}" not in out:
            out[f"s_{v}"], out[f"w_{v}"] = 0.0, np.nan
    out["cnv_region_genes"] = per_view["cnv"]["region_genes"].reindex(genes)
    # CNV support is region-level: every gene in a region shares one weight.
    # Discount it by region size so large regions do not flood the ranking.
    s_cnv_adj = out["s_cnv"] / np.sqrt(out["cnv_region_genes"].fillna(1))
    # RNA-ANCHORED support: a gene is a candidate only if its own expression
    # carries the factor (s_rna >= LAYER_SUPPORT). Methylation and CNV then
    # add support. Without this, passenger genes that merely sit in a CNV
    # region or near a methylation probe rank highly.
    out["rna_supported"] = out["s_rna"] >= LAYER_SUPPORT
    out["data_support"] = out["s_rna"] + out["s_methylation"] + s_cnv_adj
    out["n_layers"] = (out[[f"s_{v}" for v in VIEWS]] >= LAYER_SUPPORT).sum(axis=1)
    out["layers"] = out.apply(lambda r: "+".join(
        v for v in VIEWS if r[f"s_{v}"] >= LAYER_SUPPORT), axis=1)
    return out.sort_values("data_support", ascending=False)


# ---- 2. gene sets and enrichment ------------------------------------------------
def load_gmt_text(text):
    sets = {}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        genes = {g.split(",")[0].strip() for g in parts[1:] if g.strip()}
        genes.discard("")
        # Reactome lines: name \t id \t genes...; Enrichr: name \t \t genes...
        genes = {g for g in genes if not g.startswith("R-HSA-")}
        sets[parts[0]] = genes
    return sets


def fetch_reactome():
    def get():
        for url in ["https://reactome.org/download/current/ReactomePathways.gmt.zip",
                    "https://download.reactome.org/current/ReactomePathways.gmt.zip"]:
            try:
                raw = http(url)
                with zipfile.ZipFile(io.BytesIO(raw)) as z:
                    text = z.read(z.namelist()[0]).decode("utf-8")
                record("Reactome pathways (GMT)", url, raw)
                return text
            except Exception as e:                   # noqa: BLE001
                print(f"    Reactome {url}: {e}")
        return ""
    return load_gmt_text(cached_json("reactome_gmt", get))


def fetch_hallmark():
    url = ("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text"
           "&libraryName=MSigDB_Hallmark_2020")

    def get():
        try:
            raw = http(url)
            record("MSigDB Hallmark 2020 (via Enrichr library)", url, raw)
            return raw.decode("utf-8")
        except Exception as e:                       # noqa: BLE001
            print(f"    Hallmark: {e}")
            return ""
    return load_gmt_text(cached_json("hallmark_gmt", get))


def locus_share(genes, gene_pos, bin_mb=2):
    """Share of genes (with known position) that fall in the single most
    common 2-Mb genomic bin. Near 1 = the 'pathway' is one gene cluster."""
    loci = [f"{gene_pos[g][0]}:{int(gene_pos[g][1] // (bin_mb * 1e6))}"
            for g in genes if g in gene_pos]
    if len(loci) < 2:
        return np.nan, len(loci)
    return pd.Series(loci).value_counts().iloc[0] / len(loci), len(set(loci))


def ora(top_genes, universe, gene_sets):
    """Hypergeometric over-representation with BH-FDR."""
    top = set(top_genes) & universe
    N, n = len(universe), len(top)
    rows = []
    for name, gs in gene_sets.items():
        gs = gs & universe
        K = len(gs)
        if K < GS_MIN or K > GS_MAX:
            continue
        k = len(gs & top)
        # Zero-overlap eligible sets have p=1 and remain in the BH family.
        p = stats.hypergeom.sf(k - 1, N, K, n)
        rows.append({"term": name, "set_size": K, "overlap": k,
                     "expected": round(n * K / N, 2), "p": p,
                     "genes": ",".join(sorted(gs & top))})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.sort_values("p")
    m = len(df)
    q = df["p"].values * m / np.arange(1, m + 1)
    df["fdr"] = np.minimum(np.minimum.accumulate(q[::-1])[::-1], 1)
    return df


# ---- 3. literature (PubMed) -----------------------------------------------------
def _norm(s):
    return re.sub(r"[-_.\s]", "", str(s)).upper()


def hgnc_aliases(genes, other_symbols=()):
    """Previous symbols and aliases from HGNC (e.g. SERPINB3 -> SCCA1, SCCA).
    Dropped: short or ambiguous aliases, and aliases that equal ANOTHER
    gene's official symbol (e.g. IRX1 alias 'IRX-5' = gene IRX5), which
    would pull in that other gene's papers."""
    taken = {_norm(s) for s in other_symbols}
    path = CACHE / "hgnc_aliases.json"
    cache = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    todo = [g for g in genes if g not in cache]
    if todo:
        print(f"    HGNC aliases: {len(todo)} genes to query")
    ok = 0
    for g in todo:
        try:
            raw = http(f"https://rest.genenames.org/fetch/symbol/{urllib.parse.quote(g)}",
                       headers={"Accept": "application/json"}, tries=2)
            docs = json.loads(raw)["response"]["docs"]
            al = []
            if docs:
                al = list(docs[0].get("alias_symbol", [])) + \
                     list(docs[0].get("prev_symbol", []))
            cache[g] = al
            ok += 1
            time.sleep(0.1)
        except Exception as e:                       # noqa: BLE001
            print(f"    HGNC failed for {g}: {e}")
    path.write_text(json.dumps(cache), encoding="utf-8")
    if ok:
        record("HGNC REST (aliases, previous symbols)",
               "https://rest.genenames.org/fetch/symbol/")
    out = {}
    for g in genes:
        keep = []
        for a in cache.get(g, []):
            a = str(a).strip()
            if (len(a) >= 4 or (len(a) == 3 and any(ch.isdigit() for ch in a))) \
                    and a.upper() not in AMBIGUOUS and a.upper() != g.upper() \
                    and _norm(a) not in taken - {_norm(g)} \
                    and re.fullmatch(r"[A-Za-z0-9\-\.]+", a):
                keep.append(a)
        out[g] = keep[:6]
    return out


def pubmed_count(term):
    key = os.environ.get("NCBI_API_KEY", "")
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "rettype": "count",
                                "retmode": "json", **({"api_key": key} if key else {})})
    raw = http(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{q}")
    time.sleep(0.11 if key else 0.36)      # NCBI: 3 req/s without key, 10 with
    return int(json.loads(raw)["esearchresult"]["count"])


def literature(genes, aliases):
    path = CACHE / "pubmed_counts_v2.json"     # v2: alias-aware queries
    cache = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    # re-query a gene if its search names changed (e.g. an alias was dropped)
    todo = [g for g in genes if g not in cache
            or cache[g].get("names") != [g] + aliases.get(g, [])]
    if todo:
        print(f"    PubMed: {len(todo)} genes to query (cached: {len(cache)})")
    for i, g in enumerate(todo):
        try:
            names = [g] + aliases.get(g, [])
            gq = "(" + " OR ".join(f'"{n}"[tiab]' for n in names) + ")"
            cache[g] = {"lung": pubmed_count(f"{gq} AND {LUNG}"),
                        "prog": pubmed_count(f"{gq} AND {LUNG} AND {PROG}"),
                        "names": names}
        except Exception as e:                       # noqa: BLE001
            print(f"    PubMed failed for {g}: {e}")
            continue
        if (i + 1) % 25 == 0:
            path.write_text(json.dumps(cache), encoding="utf-8")
            print(f"    PubMed: {i + 1}/{len(todo)}")
    path.write_text(json.dumps(cache), encoding="utf-8")
    if todo:
        record("PubMed E-utilities esearch (counts)",
               "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
               note=f"lung query: {LUNG}; progression query adds: {PROG}")
    return cache


# ---- 4. curated evidence: Open Targets, STRING -----------------------------------
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
OT_QUERY = """
query assoc($id: String!, $index: Int!, $size: Int!) {
  disease(efoId: $id) {
    id name
    associatedTargets(page: {index: $index, size: $size}) {
      count
      rows { score target { approvedSymbol } }
    }
  }
}"""


def fetch_open_targets():
    def get():
        for efo in ["EFO_0003060", "MONDO_0005233"]:   # NSCLC (EFO, MONDO)
            for size in [3000, 500]:
                try:
                    scores, index, total, name = {}, 0, None, ""
                    while total is None or index * size < total:
                        body = json.dumps({"query": OT_QUERY, "variables":
                                           {"id": efo, "index": index, "size": size}})
                        raw = http(OT_URL, data=body.encode(),
                                   headers={"Content-Type": "application/json"})
                        d = json.loads(raw)
                        dis = (d.get("data") or {}).get("disease")
                        if not dis:
                            raise ValueError(d.get("errors", "no disease"))
                        at = dis["associatedTargets"]
                        total, name = at["count"], dis["name"]
                        for r in at["rows"]:
                            scores[r["target"]["approvedSymbol"]] = r["score"]
                        index += 1
                        if index > 100:
                            break
                    if scores:
                        record("Open Targets Platform GraphQL", OT_URL,
                               note=f"{efo} ({name}): {len(scores)} targets")
                        return {"disease": f"{efo} {name}", "scores": scores}
                except Exception as e:               # noqa: BLE001
                    print(f"    Open Targets {efo} size {size}: {e}")
        return {"disease": "", "scores": {}}
    return cached_json("open_targets_nsclc", get)


def fetch_string(genes, tag):
    def get():
        url = "https://string-db.org/api/tsv/network"
        data = urllib.parse.urlencode({
            "identifiers": "\r".join(genes), "species": 9606,
            "required_score": STRING_MIN_SCORE,
            "caller_identity": "nsclc-multiomics-hackathon"}).encode()
        try:
            raw = http(url, data=data)
            record("STRING v12 network API", url, note=f"{tag}: {len(genes)} genes")
            df = pd.read_csv(io.BytesIO(raw), sep="\t")
            return df[["preferredName_A", "preferredName_B", "score"]].values.tolist()
        except Exception as e:                       # noqa: BLE001
            print(f"    STRING ({tag}): {e}")
            return []
    return cached_json(f"string_{tag}", get)


# ---- helpers ----------------------------------------------------------------------
def pct_rank(s):
    return s.rank(pct=True, method="average")


def is_ambiguous(g):
    return g.upper() in AMBIGUOUS or len(g) <= 2 or g.startswith(("LOC", "ORF"))


def tag(row):
    if row["ambiguous_symbol"]:
        return "ambiguous symbol"
    if pd.isna(row["pubmed_progression"]):
        # never call a gene "novel" just because the literature lookup failed
        return "no literature data"
    if row["pubmed_progression"] >= KNOWN_PROG_PAPERS or row["ot_score"] >= KNOWN_OT_SCORE:
        return "known"
    if row["in_enriched_pathway"] or row["string_link_known"] or \
            row["ot_score"] >= PLAUSIBLE_OT_SCORE:
        return "plausible"
    # below every curated-evidence threshold. NOT "novel": a gene can still
    # have some literature (see pubmed columns) - verify before claiming.
    return "limited evidence"


# ---- main --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Rank genes and evidence for selected MOFA factors."
    )
    ap.add_argument("--factors", nargs="*", default=None,
                    help="factors to mine (default: stable, stage FDR < 0.05, "
                         "no strong site effect)")
    ap.add_argument("--top", type=int, default=100,
                    help="top genes per factor for literature/curation")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(exist_ok=True)

    ann = pd.read_csv(MOFA / "factor_annotation.tsv", sep="\t")
    ann["stable"] = ann["stable"].astype(str) == "True"
    batch = (ann["site_fdr"] < 0.05) & (ann["site_excess"] > 0.05)
    if args.factors:
        factors = args.factors
    else:
        factors = ann.loc[ann["stable"] & (ann["stage_fdr"] < 0.05) & ~batch,
                          "factor"].tolist()
    ann = ann.set_index("factor")

    W = {v: pd.read_csv(MOFA / f"weights_{v}.tsv", sep="\t", index_col=0)
         for v in VIEWS}
    p2g, gene_pos = probe_to_genes()
    cg = pd.read_csv(PROC / "cnv_feature_groups.tsv", sep="\t")
    cnv_groups = {f: g.split(",") for f, g in zip(cg["feature"], cg["genes"])}

    log("# Knowledge-mining report (step 5b)\n")
    log(f"- factors mined: {', '.join(factors)}")
    log(f"- top {TOP_ORA} genes per factor for enrichment; top {args.top} for "
        "literature and curated evidence")
    log("- candidates are RNA-anchored: a gene must carry the factor in its own "
        f"expression (scaled RNA weight >= {LAYER_SUPPORT}); methylation and "
        "CNV add support (CNV discounted by region size)")
    log("- background (universe) for enrichment: genes measured in the RNA "
        "layer (the layer candidates are drawn from), not the whole genome")
    log("- literature queries include HGNC aliases and previous symbols\n")

    print("fetching knowledge resources (cached after first run) ...")
    libs = {"Reactome": fetch_reactome(), "Hallmark": fetch_hallmark()}
    ot = fetch_open_targets()
    ot_scores = ot["scores"]
    for k, v in libs.items():
        log(f"- {k}: {len(v)} gene sets" + ("" if v else "  (UNAVAILABLE)"))
    log(f"- Open Targets: {len(ot_scores)} targets for {ot['disease'] or 'NSCLC'}"
        + ("" if ot_scores else "  (UNAVAILABLE)"))
    known_drivers = [g for g, _ in sorted(ot_scores.items(),
                                          key=lambda x: -x[1])[:N_KNOWN_DRIVERS]]
    log(f"- network anchors: top {len(known_drivers)} Open Targets NSCLC genes"
        + (f" (e.g. {', '.join(known_drivers[:8])})" if known_drivers else ""))
    log("")

    all_cands, all_enrich, edges = [], [], []
    for f in factors:
        log(f"## {f}\n")
        a = ann.loc[f]
        late_dir = -1 if a["stage_auc"] < 0.5 else 1   # +1: higher factor = later stage
        r2 = {c: a[c] for c in ann.columns if c.startswith("r2_")}
        log(f"- stage AUC {a['stage_auc']:.3f} (LUAD {a['stage_auc_LUAD']:.3f}, "
            f"LUSC {a['stage_auc_LUSC']:.3f}); survival HR/SD "
            f"{a['os_hr_per_sd']:.2f} (FDR {a['os_fdr']:.2g})")
        log("- variance explained: " + ", ".join(
            f"{k[3:]} {v:.1f}%" for k, v in r2.items()))

        sup_all = gene_support(f, W, p2g, cnv_groups)
        sup_all.assign(factor=f).to_csv(OUT / f"gene_support_{f}.tsv", sep="\t")
        sup = sup_all[sup_all["rna_supported"]].copy()
        universe = set(W["rna"].index)
        log(f"- RNA-supported genes: {len(sup)} of {len(universe)} measured")
        # direction relative to LATE stage (RNA, methylation)
        for v in ["rna", "methylation"]:
            sup[f"{v}_higher_in"] = np.where(
                sup[f"w_{v}"].isna(), "",
                np.where(np.sign(sup[f"w_{v}"]) * late_dir > 0, "late", "early"))

        # enrichment
        top_ora = list(sup.index[:TOP_ORA])
        enr = []
        for lib, gs in libs.items():
            if not gs:
                continue
            e = ora(top_ora, universe, gs)
            if not e.empty:
                e.insert(0, "library", lib)
                e.insert(0, "factor", f)
                enr.append(e)
        enr = pd.concat(enr) if enr else pd.DataFrame()
        if not enr.empty:
            ls = [locus_share(g.split(","), gene_pos) for g in enr["genes"]]
            enr["max_locus_share"] = [round(x[0], 2) if pd.notna(x[0]) else np.nan
                                      for x in ls]
            enr["n_loci"] = [x[1] for x in ls]
            # positional cluster: most overlapping genes sit in one 2-Mb bin
            enr["positional_cluster"] = enr["max_locus_share"] >= 0.5
        sig = enr[enr["fdr"] < 0.05] if not enr.empty else enr
        if not enr.empty:
            all_enrich.append(enr)
            n_pc = int(sig["positional_cluster"].fillna(False).sum())
            log(f"\n**Enriched gene sets** (FDR < 0.05: {len(sig)}, of which "
                f"{n_pc} driven mainly by one genomic locus; top 10 shown)\n")
            show = enr.head(10)[["library", "term", "overlap", "set_size",
                                 "expected", "fdr", "max_locus_share",
                                 "positional_cluster"]].copy()
            show["fdr"] = show["fdr"].map(lambda x: f"{x:.2g}")
            log(show.to_markdown(index=False) + "\n")
            sig = sig[~sig["positional_cluster"].fillna(False)]  # for tags
        enriched_genes = set()
        for gl in sig.get("genes", []):
            enriched_genes |= set(gl.split(","))

        # candidates: literature + curated evidence
        cand = sup.head(args.top).copy()
        aliases = hgnc_aliases(list(cand.index),
                               other_symbols=set(W["rna"].index) | set(ot_scores))
        lit = literature(list(cand.index), aliases)
        cand["aliases_searched"] = [",".join(aliases.get(g, [])) for g in cand.index]
        cand["pubmed_lung"] = [lit.get(g, {}).get("lung", np.nan) for g in cand.index]
        cand["pubmed_progression"] = [lit.get(g, {}).get("prog", np.nan)
                                      for g in cand.index]
        cand["ot_score"] = [ot_scores.get(g, 0.0) for g in cand.index]
        net = fetch_string(sorted(set(cand.index) | set(known_drivers)), f)
        kd = set(known_drivers)
        link = {g: set() for g in cand.index}
        for a_, b_, s_ in net:
            if a_ in link and b_ in kd and a_ != b_:
                link[a_].add(b_)
            if b_ in link and a_ in kd and a_ != b_:
                link[b_].add(a_)
            if a_ in link and b_ in link:
                edges.append(("gene", a_, "gene", b_, "STRING", s_))
        cand["string_known_partners"] = [",".join(sorted(link[g])) for g in cand.index]
        cand["string_link_known"] = cand["string_known_partners"] != ""
        cand["in_enriched_pathway"] = cand.index.isin(enriched_genes)
        cand["ambiguous_symbol"] = [is_ambiguous(g) for g in cand.index]
        cand["evidence"] = cand.apply(tag, axis=1)
        # GAP: data-support percentile minus literature percentile, computed
        # only over genes with literature data and unambiguous symbols
        ok = cand["pubmed_progression"].notna() & ~cand["ambiguous_symbol"]
        cand["gap_score"] = np.nan
        if ok.any():
            cand.loc[ok, "gap_score"] = (
                pct_rank(cand.loc[ok, "data_support"])
                - pct_rank(np.log1p(cand.loc[ok, "pubmed_progression"]))).round(3)
        n_missing = int(cand["pubmed_progression"].isna().sum())
        if n_missing:
            log(f"- WARNING: PubMed counts missing for {n_missing} genes; they "
                "are tagged 'no literature data' and excluded from gap scores. "
                "Re-run to retry (successful lookups are cached).")
        if not ot_scores:
            log("- WARNING: Open Targets unavailable; 'known' and 'plausible' "
                "tags are incomplete for this run.")
        if not net:
            log("- note: STRING returned no high-confidence interactions for "
                "these genes (or was unreachable); network-based 'plausible' "
                "tags could not be assigned.")
        cand.insert(0, "factor", f)
        all_cands.append(cand.reset_index())

        counts = cand["evidence"].value_counts().to_dict()
        log(f"**Evidence tags** (top {args.top} genes): " + ", ".join(
            f"{k} {v}" for k, v in counts.items()))
        cols = ["gene", "layers", "data_support", "rna_higher_in", "pubmed_lung",
                "pubmed_progression", "ot_score", "string_known_partners",
                "evidence", "gap_score", "aliases_searched"]
        top_known = cand.reset_index().head(10)[cols]
        log("\n**Top 10 genes by multi-omics support**\n")
        log(top_known.round(3).to_markdown(index=False) + "\n")
        gap = cand.reset_index().dropna(subset=["gap_score"]) \
            .sort_values("gap_score", ascending=False).head(10)[cols]
        log("**Top 10 under-studied candidates** (highest gap score: strong "
            "data support, little progression literature)\n")
        log(gap.round(3).to_markdown(index=False) + "\n")

        # knowledge-graph edges
        for g, r in cand.iterrows():
            edges.append(("factor", f, "gene", g, "weight", round(r["data_support"], 4)))
            if r["ot_score"] > 0:
                edges.append(("gene", g, "disease", "NSCLC", "open_targets",
                              round(r["ot_score"], 4)))
        for _, e in sig.iterrows():
            for g in e["genes"].split(","):
                if g in cand.index:
                    edges.append(("gene", g, "pathway", e["term"], "member", 1))

    # ---- save outputs --------------------------------------------------------------
    cands = pd.concat(all_cands) if all_cands else pd.DataFrame()
    cands.to_csv(OUT / "candidates.tsv", sep="\t", index=False)
    if all_enrich:
        pd.concat(all_enrich).to_csv(OUT / "enrichment.tsv", sep="\t", index=False)
    E = pd.DataFrame(edges, columns=["source_type", "source", "target_type",
                                     "target", "relation", "weight"]).drop_duplicates()
    E.to_csv(OUT / "knowledge_graph_edges.tsv", sep="\t", index=False)
    write_graphml(E, OUT / "knowledge_graph.graphml")
    sources_path = OUT / "knowledge_sources.tsv"
    if SOURCES or not sources_path.exists():
        pd.DataFrame(SOURCES).to_csv(sources_path, sep="\t", index=False)
    plot_gap(cands, OUT / "support_vs_literature.png")

    log("## Method notes and limitations\n")
    log(f"- Evidence tags: known = >= {KNOWN_PROG_PAPERS} lung-cancer progression "
        f"papers or Open Targets score >= {KNOWN_OT_SCORE}; plausible = in an "
        f"enriched gene set (FDR < 0.05), a STRING partner (score >= "
        f"{STRING_MIN_SCORE}) of a top-{N_KNOWN_DRIVERS} Open Targets NSCLC gene, "
        f"or Open Targets score >= {PLAUSIBLE_OT_SCORE}; limited evidence = none "
        "of these (not the same as novel - such genes may still have some "
        "literature and should be checked manually before any novelty claim).")
    log("- Aliases that equal another gene's official symbol (e.g. 'IRX-5' "
        "listed for IRX1) are excluded from literature queries.")
    log("- PubMed counts use title/abstract mentions of the gene symbol and its "
        "HGNC aliases/previous symbols (aliases shorter than 4 characters or "
        "ambiguous are dropped); symbols that are common words or very short "
        "are flagged 'ambiguous symbol' and excluded from novelty claims. "
        "Counts measure attention, not truth.")
    log("- Candidates are RNA-anchored; CNV features stand for several "
        "neighbouring genes, so CNV support is region-level and discounted by "
        "the square root of region size.")
    log("- Enriched gene sets whose overlapping genes mostly fall in one 2-Mb "
        "genomic bin are flagged `positional_cluster` (co-located gene families "
        "such as CYP4F or SERPINB can look like pathway enrichment) and are not "
        "used for 'plausible' tags.")
    log("- Methylation probes are assigned to their annotated genes; promoter vs "
        "gene-body effects are not distinguished.")
    log("- Candidates are hypotheses for follow-up, not validated biomarkers.\n")
    log("Outputs: results/knowledge/candidates.tsv, enrichment.tsv, "
        "knowledge_graph_edges.tsv, knowledge_graph.graphml, "
        "knowledge_sources.tsv, support_vs_literature.png")
    (OUT / "knowledge_report.md").write_text("\n".join(LOG), encoding="utf-8")
    print(f"\nReport: {OUT / 'knowledge_report.md'}")


def write_graphml(E, path):
    def esc(s):
        return (str(s).replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;").replace('"', "&quot;"))
    nodes = {}
    for _, r in E.iterrows():
        nodes[f"{r.source_type}:{r.source}"] = r.source_type
        nodes[f"{r.target_type}:{r.target}"] = r.target_type
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
             '<key id="type" for="node" attr.name="type" attr.type="string"/>',
             '<key id="rel" for="edge" attr.name="relation" attr.type="string"/>',
             '<key id="w" for="edge" attr.name="weight" attr.type="double"/>',
             '<graph id="nsclc_knowledge" edgedefault="undirected">']
    for n, t in nodes.items():
        lines.append(f'<node id="{esc(n)}"><data key="type">{esc(t)}</data></node>')
    for i, r in E.iterrows():
        lines.append(f'<edge id="e{i}" source="{esc(r.source_type)}:{esc(r.source)}" '
                     f'target="{esc(r.target_type)}:{esc(r.target)}">'
                     f'<data key="rel">{esc(r.relation)}</data>'
                     f'<data key="w">{float(r.weight)}</data></edge>')
    lines += ["</graph>", "</graphml>"]
    path.write_text("\n".join(lines), encoding="utf-8")


def plot_gap(cands, path):
    if cands.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    facs = cands["factor"].unique()
    fig, axes = plt.subplots(1, len(facs), figsize=(5 * len(facs), 4.5), squeeze=False)
    colors = {"known": "#4C72B0", "plausible": "#DD8452",
              "limited evidence": "#55A868", "ambiguous symbol": "#BBBBBB"}
    for ax, f in zip(axes[0], facs):
        d = cands[(cands["factor"] == f) & cands["pubmed_progression"].notna()]
        for t, c in colors.items():
            m = d["evidence"] == t
            ax.scatter(d.loc[m, "data_support"], np.log1p(d.loc[m, "pubmed_progression"]),
                       s=14, c=c, label=t, alpha=0.8)
        for _, r in d.dropna(subset=["gap_score"]).nlargest(8, "gap_score").iterrows():
            ax.annotate(r["gene"], (r["data_support"], math.log1p(r["pubmed_progression"])),
                        fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("multi-omics data support")
        ax.set_ylabel("log(1 + progression papers)")
        ax.set_title(f"{f}: data support vs literature")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()