#!/usr/bin/env python3
"""Download selected Xena inputs and maintain their checksum manifest.

Run from the repository root. Raw files are written to ``data/raw/`` and
their URLs, sizes, timestamps, and SHA-256 digests to
``data/raw_manifest.tsv``.
"""
import argparse
import csv
import datetime as dt
import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

HUB = "https://tcga.xenahubs.net"
SOURCES = Path("data/sources.tsv")
RAW = Path("data/raw")
MANIFEST = Path("data/raw_manifest.tsv")
FIELDS = ["source_id", "dataset", "url", "file", "bytes", "sha256",
          "downloaded_utc"]


def selected_xena_sources():
    with open(SOURCES, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return [r for r in rows
            if r["status"] == "selected" and "Xena" in r["repository"]]


def candidate_urls(dataset, hub):
    # Xena serves most matrices gzipped; clinical matrices are often plain.
    # The "/" in dataset names is usually URL-encoded, but try it raw too.
    enc = dataset.replace("/", "%2F")
    return [f"{hub}/download/{enc}.gz", f"{hub}/download/{enc}",
            f"{hub}/download/{dataset}.gz", f"{hub}/download/{dataset}"]


def local_name(dataset, url):
    """Return a local filename, prefixing cohort-specific datasets."""
    ext = ".gz" if url.endswith(".gz") else ".tsv"
    prefix, base = dataset.split("/", 1)
    if ".sampleMap" in prefix:
        cohort = prefix.split(".")[1]                 # "LUAD" / "LUSC"
        return f"{cohort}_{base.removeprefix(cohort + '_')}{ext}"
    return f"{prefix}_{base}{ext}"                    # e.g. probeMap_...


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def download(url, dest):
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as out:
        total = int(r.headers.get("Content-Length") or 0)
        done = 0
        while block := r.read(1 << 20):
            out.write(block)
            done += len(block)
            if total:
                print(f"\r    {done / 1e6:7.1f} / {total / 1e6:.1f} MB",
                      end="", flush=True)
            else:
                print(f"\r    {done / 1e6:7.1f} MB", end="", flush=True)
    print()
    tmp.replace(dest)  # only a complete download gets the real name


def load_manifest():
    if not MANIFEST.exists():
        return {}
    with open(MANIFEST, newline="") as f:
        return {r["source_id"]: r for r in csv.DictReader(f, delimiter="\t")}


def write_manifest(entries):
    with open(MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for sid in sorted(entries):
            w.writerow(entries[sid])


def verify(manifest):
    bad = 0
    for sid, e in sorted(manifest.items()):
        p = RAW / e["file"]
        if not p.exists():
            print(f"MISSING   {sid:16s} {e['file']}")
            bad += 1
        elif sha256(p) != e["sha256"]:
            print(f"MISMATCH  {sid:16s} {e['file']}")
            bad += 1
        else:
            print(f"OK        {sid:16s} {e['file']}")
    return bad


def main():
    ap = argparse.ArgumentParser(
        description="Download selected Xena inputs or verify their checksums."
    )
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--hub", default=HUB)
    args = ap.parse_args()

    manifest = load_manifest()
    if args.verify:
        sys.exit(1 if verify(manifest) else 0)

    RAW.mkdir(parents=True, exist_ok=True)
    sources = selected_xena_sources()
    print(f"{len(sources)} selected Xena sources in {SOURCES}\n")

    failed = []
    for s in sources:
        sid, ds = s["source_id"], s["dataset_or_accession"]
        prev = manifest.get(sid)
        if prev and not args.force and (RAW / prev["file"]).exists():
            print(f"skip  {sid:16s} already downloaded ({prev['file']})")
            continue

        print(f"get   {sid:16s} {ds}")
        for url in candidate_urls(ds, args.hub):
            fname = local_name(ds, url)
            dest = RAW / fname
            try:
                download(url, dest)
            except urllib.error.HTTPError as e:
                # Xena redirects to Amazon S3, which answers 403 (not 404)
                # for files that don't exist. Treat both as "try next form".
                if e.code in (403, 404):
                    print(f"      {e.code} for {url.rsplit('/', 1)[-1]}, "
                          "trying next URL form")
                    continue
                raise
            manifest[sid] = {
                "source_id": sid, "dataset": ds, "url": url, "file": fname,
                "bytes": dest.stat().st_size, "sha256": sha256(dest),
                "downloaded_utc": dt.datetime.now(dt.timezone.utc)
                                    .strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
            write_manifest(manifest)  # save after each file
            break
        else:
            print(f"      FAILED: no URL worked for {ds}", file=sys.stderr)
            failed.append(sid)

    print(f"\nManifest: {MANIFEST}")
    if failed:
        print("Failed:", ", ".join(failed), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()