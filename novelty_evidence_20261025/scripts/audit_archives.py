"""Arithmetic audit of archived aggregate CSVs; no new experiment or inference.

Inputs are read through common.gitbytes at the pinned archival commit.  The
factorial table is a comparison target only, never an input to the recomputed
contrasts.  rhoB uses the archived squared bias (bias2), not its square root.
"""
import csv
import hashlib
import io
import math

from common import ROOT, gitbytes, writejson


COMMIT = "198f04b76d54c5ca77cd7b59df3a5653d6001a3e"
DECOMPOSITION = "combo_bias_variance_20261024/analysis/decomposition.csv"
FACTORIAL = "combo_bias_variance_20261024/analysis/factorial.csv"
METHODS = {
    "A": "FullCalRefit",
    "B": "Standard",
    "C": "QScanRefit",
    "D": "QScan_Standard",
}
METRICS = {"mse": "dmse", "bias2": "dbias2", "var": "dvar"}
ATOL = 1e-12
RTOL = 1e-10


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite {key}: {row.get('scene')} / {row.get('method')}")
    return value


def table(path):
    raw = gitbytes(path)
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if not rows:
        raise ValueError(f"Empty archived CSV: {path}")
    return rows, {
        "path": path,
        "commit": COMMIT,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
        "total_rows": len(rows),
        "region_all_rows": sum(row["region"] == "all" for row in rows),
    }


def ratio(numerator, denominator):
    # A zero contrast has no defined ratio. Do not emit NaN/Infinity or silently
    # substitute a claimed effect. Preserve nonzero, even small denominators.
    return numerator / denominator if denominator != 0.0 else None


def main():
    decomposition, source_d = table(DECOMPOSITION)
    factorial, source_f = table(FACTORIAL)
    all_rows = [row for row in decomposition if row["region"] == "all"]
    observed_methods = sorted({row["method"] for row in all_rows})
    missing_methods = set(METHODS.values()) - set(observed_methods)
    if missing_methods:
        raise ValueError(f"Required methods absent: {sorted(missing_methods)}")

    scenes = {}
    for row in all_rows:
        scene = row["scene"]
        by_method = scenes.setdefault(scene, {})
        if row["method"] in by_method:
            raise ValueError(f"Duplicate decomposition row: {scene} / {row['method']}")
        by_method[row["method"]] = row

    reference = {}
    for row in factorial:
        if row["region"] != "all":
            continue
        if row["scene"] in reference:
            raise ValueError(f"Duplicate factorial row: {row['scene']}")
        reference[row["scene"]] = row
    if not scenes or set(scenes) != set(reference):
        raise ValueError("The region=all scene sets differ, or no scenes exist")

    results = []
    comparisons = []
    decomposition_residuals = []
    for scene in sorted(scenes):
        by_method = scenes[scene]
        absent = set(METHODS.values()) - set(by_method)
        if absent:
            raise ValueError(f"{scene}: missing methods {sorted(absent)}")
        cells = {label: by_method[method] for label, method in METHODS.items()}
        original = reference[scene]
        for key in ("model", "snr_db", "Q_percent", "R"):
            values = [row[key] for row in cells.values()] + [original[key]]
            if len(set(values)) != 1:
                raise ValueError(f"{scene}: inconsistent metadata {key}: {values}")

        result = {key: original[key] for key in ("scene", "model", "snr_db", "Q_percent", "R", "region")}
        for metric, suffix in METRICS.items():
            values = {label: number(row, metric) for label, row in cells.items()}
            for label, value in values.items():
                result[f"{label}_{metric}"] = value
            raw_delta = values["C"] - values["A"]
            coh_delta = values["D"] - values["B"]
            interaction = coh_delta - raw_delta
            for prefix, value in (("raw", raw_delta), ("coh", coh_delta), ("int", interaction)):
                column = f"{prefix}_{suffix}"
                archived = number(original, column)
                difference = value - archived
                matches = math.isclose(value, archived, rel_tol=RTOL, abs_tol=ATOL)
                result[column] = value
                result[f"archived_{column}"] = archived
                result[f"difference_{column}"] = difference
                comparisons.append({
                    "scene": scene,
                    "column": column,
                    "recomputed": value,
                    "archived": archived,
                    "difference": difference,
                    "matches_tolerance": matches,
                })

        result["rhoV_numerator"] = result["B_var"] - result["D_var"]
        result["rhoV_denominator"] = result["A_var"] - result["C_var"]
        result["rhoV"] = ratio(result["rhoV_numerator"], result["rhoV_denominator"])
        result["rhoB_numerator"] = result["D_bias2"] - result["B_bias2"]
        result["rhoB_denominator"] = result["C_bias2"] - result["A_bias2"]
        result["rhoB"] = ratio(result["rhoB_numerator"], result["rhoB_denominator"])
        result["factorial_matches"] = all(
            comparison["matches_tolerance"]
            for comparison in comparisons
            if comparison["scene"] == scene
        )
        for label in METHODS:
            residual = result[f"{label}_mse"] - result[f"{label}_bias2"] - result[f"{label}_var"]
            result[f"{label}_decomposition_residual"] = residual
            decomposition_residuals.append(abs(residual))
        results.append(result)

    mismatches = [row for row in comparisons if not row["matches_tolerance"]]
    report = {
        "audit_type": "arithmetic_recalculation_of_archived_aggregate_tables",
        "scope": "region=all only",
        "limitations": [
            "This is an archive arithmetic audit, not a new experiment or independent scientific replication.",
            "Aggregate CSVs do not permit re-estimation of sampling uncertainty, pairing, or significance.",
            "rhoB denotes a ratio of squared-bias contrasts (the archived bias2), not signed bias.",
            "Ratios are derived here; they are not present as columns in the archived factorial table.",
            "A zero ratio denominator is reported as null in JSON and an empty CSV cell.",
        ],
        "sources": [source_d, source_f],
        "methods": METHODS,
        "observed_methods": observed_methods,
        "formulas": {
            "rawDelta": "C - A, separately for mse, bias2, and var",
            "cohDelta": "D - B, separately for mse, bias2, and var",
            "I": "(D - B) - (C - A), separately for mse, bias2, and var",
            "rhoV": "(varB - varD) / (varA - varC)",
            "rhoB": "(bias2D - bias2B) / (bias2C - bias2A)",
        },
        "tolerance": {"absolute": ATOL, "relative": RTOL},
        "scene_count": len(results),
        "comparison_count": len(comparisons),
        "mismatch_count": len(mismatches),
        "all_factorial_contrasts_match": not mismatches,
        "max_absolute_factorial_difference": max(abs(row["difference"]) for row in comparisons),
        "max_absolute_mse_minus_bias2_minus_var": max(decomposition_residuals),
        "mismatches": mismatches,
        "scenes": results,
        "factorial_comparisons": comparisons,
    }
    output_json = ROOT / "analysis" / "archive_audit.json"
    output_csv = ROOT / "analysis" / "archive_audit.csv"
    writejson(output_json, report)
    with output_csv.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    print(f"Archive arithmetic audit: {len(results)} scenes; {len(mismatches)} contrast mismatches.")
    print(output_json)
    print(output_csv)
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
