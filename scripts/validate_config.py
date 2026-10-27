#!/usr/bin/env python3
"""Preflight for SnakeMicroC experiment configurations.

Validation does not create, modify, or delete project data. It intentionally
requires descriptive sample IDs without shell metacharacters because the
legacy Snakemake shell rules pass them as space/comma-separated arguments.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import yaml

SAMPLE_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class ConfigurationError(ValueError):
    pass


def nested(cfg, dotted):
    value = cfg
    for part in dotted.split("."):
        if not isinstance(value, dict) or part not in value:
            raise ConfigurationError(f"Missing required configuration key: {dotted}")
        value = value[part]
    return value


def positive_int(cfg, dotted, allow_zero=False):
    val = nested(cfg, dotted)
    if isinstance(val, bool) or not isinstance(val, int) or val < (0 if allow_zero else 1):
        raise ConfigurationError(f"{dotted} must be a {'nonnegative' if allow_zero else 'positive'} integer")
    return val


def load_and_validate(config_path, check_paths=False):
    path = Path(config_path).resolve()
    if not path.is_file():
        raise ConfigurationError(f"Cannot read config file: {path}")
    cfg = yaml.safe_load(path.read_text(encoding="utf8"))
    if not isinstance(cfg, dict):
        raise ConfigurationError("Configuration must be a YAML mapping")

    samples = nested(cfg, "samples")
    if not isinstance(samples, dict) or not samples:
        raise ConfigurationError("samples must be a nonempty mapping")
    conditions = {}
    for sample, data in samples.items():
        if not isinstance(sample, str) or not SAMPLE_PATTERN.fullmatch(sample):
            raise ConfigurationError(f"Invalid sample ID {sample!r}: use letters, digits, dot, underscore or hyphen")
        if not isinstance(data, dict):
            raise ConfigurationError(f"Sample {sample}: expected a mapping")
        group = data.get("condition")
        if not isinstance(group, str) or not SAMPLE_PATTERN.fullmatch(group):
            raise ConfigurationError(f"Sample {sample}: invalid condition {group!r}")
        conditions[group] = conditions.get(group, 0) + 1
        for field in ("r1", "r2"):
            value = data.get(field)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Sample {sample}: missing {field}")
            if check_paths and not Path(value).expanduser().is_file():
                raise ConfigurationError(f"Sample {sample}: {field} path not found: {value}")
        if data["r1"] == data["r2"]:
            raise ConfigurationError(f"Sample {sample}: R1 and R2 must be different files")

    tmpdir = nested(cfg, "tmpdir")
    if not isinstance(tmpdir, str) or not Path(tmpdir).expanduser().is_absolute() or Path(tmpdir).resolve() == Path("/"):
        raise ConfigurationError("tmpdir must be an absolute non-root path (processing cleans its sample subdirectories)")
    resultdir = cfg.get("results_dir", "results")
    if not isinstance(resultdir, str) or not resultdir.strip() or resultdir in ("/", ".", ".."):
        raise ConfigurationError("results_dir must designate a results directory, not the workspace root")

    for field in ("reference.fasta", "reference.fai", "reference.assembly"):
        value = nested(cfg, field)
        if not isinstance(value, str) or not value.strip():
            raise ConfigurationError(f"{field} must be non-empty")
        if check_paths and field != "reference.assembly" and not Path(value).expanduser().is_file():
            raise ConfigurationError(f"Missing reference file {field}: {value}")

    for field in (
        "processing.threads", "processing.io_threads", "processing.mapq", "qc.qc_map_resolution",
        "contact_maps.base_resolution", "reproducibility.resolution",
        "compartments.primary_resolution", "boundaries.resolution", "boundaries.primary_window",
        "contacts.resolution", "contacts.min_distance_bp", "contacts.max_distance_bp",
        "contacts.prefilter_min_total_count", "contacts.prefilter_min_nonzero_samples",
        "loops.primary_resolution", "loops.min_distance_bp", "loops.max_distance_bp",
    ):
        positive_int(cfg, field)
    positive_int(cfg, "loops.quant_radius_bins", allow_zero=True)
    cmaps = nested(cfg, "contact_maps.resolutions")
    if not isinstance(cmaps, list) or not cmaps or any(type(x) is not int or x <= 0 for x in cmaps):
        raise ConfigurationError("contact_maps.resolutions must be a nonempty list of positive integers")
    if len(set(cmaps)) != len(cmaps):
        raise ConfigurationError("contact_maps.resolutions contains duplicates")
    base_res = cfg["contact_maps"]["base_resolution"]
    if any(x % base_res != 0 for x in cmaps):
        raise ConfigurationError("All contact-map resolutions must be integer multiples of base_resolution")
    required_cool_res = [
        cfg["reproducibility"]["resolution"], cfg["compartments"]["primary_resolution"],
        cfg["contacts"]["resolution"], cfg["qc"]["qc_map_resolution"],
    ]
    if any(r not in cmaps for r in required_cool_res):
        raise ConfigurationError("HiCRep, compartment, contact and QC resolutions must occur in contact_maps.resolutions")
    windows = nested(cfg, "boundaries.windows")
    if not isinstance(windows, list) or not windows or any(type(w) is not int or w <= 0 for w in windows):
        raise ConfigurationError("boundaries.windows must list positive integers")
    if cfg["boundaries"]["primary_window"] not in windows:
        raise ConfigurationError("boundaries.primary_window must be in boundaries.windows")
    if any(w % cfg["boundaries"]["resolution"] != 0 for w in windows):
        raise ConfigurationError("All boundaries.windows must be multiples of boundaries.resolution")
    for section in ("contacts", "loops"):
        a, b = cfg[section]["min_distance_bp"], cfg[section]["max_distance_bp"]
        if a >= b:
            raise ConfigurationError(f"{section}: min_distance_bp must be less than max_distance_bp")
    if cfg["contacts"]["prefilter_min_nonzero_samples"] > len(samples):
        raise ConfigurationError("contacts.prefilter_min_nonzero_samples exceeds number of samples")

    cm = cfg["contact_maps"]
    if not isinstance(cm.get("reuse_verified_existing_mcools", False), bool):
        raise ConfigurationError("contact_maps.reuse_verified_existing_mcools must be true/false")
    if cm.get("reuse_verified_existing_mcools", False):
        if int(cfg["loops"]["primary_resolution"]) != 10000 or 50000 not in cmaps:
            raise ConfigurationError("Legacy 50-kb matrix adoption requires 10-kb loop coolers and 50-kb mcool resolution")
        evidence = cm.get("existing_mcool_equivalence_tsv")
        if not isinstance(evidence, str) or not Path(evidence).is_absolute():
            raise ConfigurationError("contact_maps.existing_mcool_equivalence_tsv must be an absolute path")
    comparison = cfg.get("comparison") or {}
    if not isinstance(comparison, dict):
        raise ConfigurationError("comparison must be a mapping")
    ctrl = comparison.get("control_condition")
    test = comparison.get("test_condition")
    if ctrl and test and ctrl == test:
        raise ConfigurationError("Comparison groups must differ")
    if bool(ctrl) != bool(test):
        raise ConfigurationError("Specify both comparison groups, or omit both")
    ncontrol, ntest = conditions.get(ctrl, 0), conditions.get(test, 0)
    enabled = bool(cfg.get("analysis", {}).get("differential", True))
    ready = enabled and ncontrol >= 2 and ntest >= 2
    if cfg.get("analysis", {}).get("require_differential", False) and not ready:
        raise ConfigurationError("analysis.require_differential is true but fewer than two samples per comparison group exist")
    return {
        "config": cfg,
        "samples": len(samples),
        "conditions": conditions,
        "comparison": f"{test} versus {ctrl}" if test and ctrl else "not configured",
        "differential_enabled": ready,
        "message": "Differential analyses enabled" if ready else "QC and descriptive structure only (differential analyses skipped)",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--check-paths", action="store_true")
    args = parser.parse_args()
    try:
        result = load_and_validate(args.config, check_paths=args.check_paths)
    except ConfigurationError as exc:
        parser.error(str(exc))
    print(f"VALID: {result['samples']} samples; conditions={result['conditions']}")
    print(f"Comparison: {result['comparison']}")
    print(result["message"])


if __name__ == "__main__":
    main()
