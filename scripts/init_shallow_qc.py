#!/usr/bin/env python3
"""Create a separate shallow Micro-C pilot workspace from explicit FASTQ pairs.

Never guesses genotype, pairs multiple sequencing lanes silently, or overwrites
an existing configuration. Only the generated YAML file is created.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
import yaml

FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)(?:_S\d+)?(?:_L\d{3})?_R(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$",
    re.IGNORECASE,
)
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def discover_pairs(folder):
    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise ValueError(f"FASTQ directory does not exist: {folder}")
    samples = {}
    filenames = sorted(x for x in folder.rglob("*") if x.is_file() and
                       (x.name.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))))
    for path in filenames:
        match = FASTQ_RE.fullmatch(path.name)
        if not match:
            raise ValueError(f"Unrecognized FASTQ naming: {path.name}; expected *_R1_001.fastq.gz / *_R2_001.fastq.gz")
        sample = match.group("sample")
        if not ID_RE.fullmatch(sample):
            raise ValueError(f"Invalid sample ID derived from {path.name}: {sample}")
        mate = match.group("mate")
        slots = samples.setdefault(sample, {})
        if mate in slots:
            raise ValueError(f"Multiple R{mate} FASTQs for {sample}; merge lanes explicitly before generating a config")
        slots[mate] = path
    if not samples:
        raise ValueError(f"No FASTQ files found under {folder}")
    for sample, mates in samples.items():
        if set(mates) != {"1", "2"}:
            raise ValueError(f"Sample {sample} has incomplete mates: {list(mates)}")
    return {s: {"r1": str(mates["1"]), "r2": str(mates["2"])} for s, mates in sorted(samples.items())}


def make_config(base_config, fastq_dir, workspace):
    base = Path(base_config)
    cfg = yaml.safe_load(base.read_text(encoding="utf8"))
    if not isinstance(cfg, dict):
        raise ValueError("Base config is not a YAML mapping")
    cfg["tmpdir"] = str((Path(workspace).resolve() / "tmp"))
    cfg["results_dir"] = "results"
    pairs = discover_pairs(fastq_dir)
    cfg["samples"] = {
        sample: {"condition": "BATCH2", "clone": sample, **paths}
        for sample, paths in pairs.items()
    }
    cfg["analysis"] = {"mode": "shallow_qc", "differential": False, "require_differential": False}
    cfg["comparison"] = {}
    cfg["contact_maps"]["reuse_verified_existing_mcools"] = False
    cfg["contact_maps"].pop("existing_mcool_equivalence_tsv", None)
    cfg["contacts"]["prefilter_min_nonzero_samples"] = min(
        len(pairs), int(cfg["contacts"]["prefilter_min_nonzero_samples"])
    )
    return cfg


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", required=True,
                        help="Existing validated Micro-C YAML from the same reference genome")
    parser.add_argument("--fastq-dir", required=True)
    parser.add_argument("--workspace", required=True)
    args = parser.parse_args(argv)
    workspace = Path(args.workspace).resolve()
    target = workspace / "config_shallow_qc.yaml"
    if target.exists():
        parser.error(f"Refusing to overwrite existing configuration: {target}")
    try:
        cfg = make_config(args.base_config, args.fastq_dir, workspace)
        from validate_config import load_and_validate
        workspace.mkdir(parents=True, exist_ok=True)
        target.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf8")
        try:
            load_and_validate(target, check_paths=True)
        except Exception:
            target.unlink()
            raise
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"Shallow QC config: {target}")
    print(f"{len(cfg['samples'])} paired-end libraries, no deep-level structural analysis")
    for name, files in cfg["samples"].items():
        print(f"  {name}: {files['r1']} | {files['r2']}")


if __name__ == "__main__":
    main()
