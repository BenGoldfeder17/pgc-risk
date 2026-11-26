"""Download the PGC summary-statistic shards listed in the catalog from Hugging Face."""
from __future__ import annotations

import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .catalog import SOURCES, Source

DATA = Path(__file__).resolve().parent.parent / "data"
RAW = DATA / "raw"


def raw_dir(src: Source) -> Path:
    return RAW / src.repo.replace("/", "__")


def wanted_files(src: Source, repo_files: list[str]) -> list[str]:
    """Shard paths of this source that actually exist in the repo (a few are missing upstream)."""
    if src.file:
        return [src.file]
    pat = re.compile(rf"^data/{re.escape(src.config)}/train-(\d{{5}})-of-{src.series:05d}\.parquet$")
    lo = 0 if src.first is None else src.first
    hi = src.series - 1 if src.last is None else src.last
    out = []
    for f in repo_files:
        m = pat.match(f)
        if m and lo <= int(m.group(1)) <= hi:
            out.append(f)
    return sorted(out)


def local_files(src: Source) -> list[Path]:
    """Shards of this source already on disk."""
    base = raw_dir(src)
    if src.file:
        p = base / src.file
        return [p] if p.exists() else []
    lo = 0 if src.first is None else src.first
    hi = src.series - 1 if src.last is None else src.last
    pat = re.compile(rf"^train-(\d{{5}})-of-{src.series:05d}\.parquet$")
    out = []
    d = base / "data" / src.config
    if d.is_dir():
        for p in d.iterdir():
            m = pat.match(p.name)
            if m and lo <= int(m.group(1)) <= hi:
                out.append(p)
    return sorted(out)


def _download_one(repo: str, revision: str, filename: str, dest: Path) -> str:
    from huggingface_hub import hf_hub_download

    target = dest / filename
    if target.exists() and target.stat().st_size > 0:
        return "cached"
    last = None
    for attempt in range(8):
        try:
            hf_hub_download(repo, filename, repo_type="dataset", revision=revision, local_dir=dest)
            return "downloaded"
        except Exception as e:  # rate limits and transient network errors: back off and retry
            last = e
            time.sleep(min(90, 5 * 2 ** attempt))
    raise RuntimeError(f"{repo}/{filename}: {last}")


def fetch(keys: list[str] | None = None, workers: int = 6) -> None:
    fetch_sources([s for s in SOURCES if not keys or s.key in keys], workers)


def fetch_sources(todo: list[Source], workers: int = 6) -> None:
    from huggingface_hub import HfApi

    api = HfApi()
    listings: dict[tuple[str, str], list[str]] = {}
    for src in todo:
        k = (src.repo, src.revision)
        if k not in listings:
            listings[k] = api.list_repo_files(src.repo, repo_type="dataset", revision=src.revision)
        files = wanted_files(src, listings[k])
        dest = raw_dir(src)
        dest.mkdir(parents=True, exist_ok=True)
        done = fresh = 0
        t0 = time.time()
        print(f"[{src.key}] {src.name}: {len(files)} file(s) from {src.repo}", flush=True)
        with ThreadPoolExecutor(workers) as pool:
            futs = [pool.submit(_download_one, src.repo, src.revision, f, dest) for f in files]
            for fut in as_completed(futs):
                fresh += fut.result() == "downloaded"
                done += 1
                if done % 100 == 0 or done == len(files):
                    print(f"[{src.key}]   {done}/{len(files)}  ({fresh} new, {time.time() - t0:.0f}s)",
                          flush=True)
    print("fetch complete", flush=True)


if __name__ == "__main__":
    fetch(sys.argv[1:] or None)
