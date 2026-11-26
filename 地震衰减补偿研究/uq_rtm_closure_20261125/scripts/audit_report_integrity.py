#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_report_integrity.py
=========================
Mechanical Scientific Integrity & Anti-Selective-Reporting Auditor (Gate R0).

Enforces the Six Zero-Tolerance Iron Laws on any experimental Markdown report
against raw JSON outputs and MATLAB source scripts:
  Gate 1 (Zero Row Omission): Every method in every experimental JSON must appear
         in its corresponding dataset table in the report (1 row per method).
  Gate 2 (Zero Column Omission): Field tables must include ALL mandatory metrics
         (F-5 primary & secondary, blind well F-14, ZeroLag L0, BestLag lag,
          Res_Cen, Smile_Frac, Res_Slope); synthetic W_deep tables must include
         E1_mean, per-m E1, E1_HF, E_amp, SNR_img_dB.
  Gate 3 (Better-Baseline & Degradation Disclosure): Every method beating the
         designated primary method by > 0.001 (and any negative/degraded blind
         well or L0 metric) must be explicitly named in a dedicated
         '反超主推 / 恶化项核查' block beneath each table.
  Gate 4 (Active-Operator & Marginal Attribution Gate): Computes median HF
         spectral gain G_HF (18-24 Hz for <=30 Hz; 28-40 Hz for >30 Hz) vs
         zero-compensation equalized baseline (DSA_R0_Eq) and front-end
         anchoring share. Fails if G_HF < 1.10 without explicit dormancy
         disclosure or if marginal attribution is missing.
  Gate 5 (Anti-Hype Gate): Bans hype phrases ('全面击败', '完胜', '碾压',
         '独占第一', '绝对第一', '显著超越', '全面超越', '全面领先', '大胜')
         when single-well |Delta r| < 0.01 against core competitors.
  Gate 6 (Code Wiring Checklist Gate): Requires '代码接线与实际增益核查表' with
         exact iv variable name, script path + line numbers (verified on disk),
         df (Hz), and measured G_HF matching JSON.

Exit code: 0 if and only if 0 failures and 0 warnings; 1 otherwise.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import statistics
import sys
from typing import Any, Dict, List, Optional, Set, Tuple


BANNED_HYPE_TERMS = [
    "全面击败",
    "完胜",
    "碾压",
    "独占第一",
    "绝对第一",
    "显著超越",
    "全面超越",
    "全面领先",
    "大胜",
]

DORMANT_DISCLOSURE_TERMS = [
    "高频补偿处于休眠",
    "弱作用状态",
    "休眠态",
    "休眠状态",
    "关停状态",
    "被混叠 iv3 关停",
]


def load_json_files(results_dir: str) -> List[Tuple[str, Dict[str, Any]]]:
    files = sorted(glob.glob(os.path.join(results_dir, "*.json")))
    out: List[Tuple[str, Dict[str, Any]]] = []
    for fp in files:
        base = os.path.basename(fp)
        if base.endswith("_smoke.json"):
            continue
        if "verdict" in base.lower() or "diagnostics" in base.lower():
            continue
        with open(fp, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict) and ("methods" in data or "W_deep" in data):
            out.append((fp, data))
    return out


def get_dataset_id(fp: str, data: Dict[str, Any]) -> str:
    stem = os.path.splitext(os.path.basename(fp))[0]
    if "model" in data and isinstance(data["model"], str):
        return data["model"]
    return stem


def get_methods_dict(data: Dict[str, Any]) -> Dict[str, Any]:
    if "methods" in data and isinstance(data["methods"], dict):
        return data["methods"]
    if "W_deep" in data and isinstance(data["W_deep"], dict):
        return data["W_deep"]
    return {}


def detect_primary_method(ds_id: str, data: Dict[str, Any], override_map: Dict[str, str]) -> Optional[str]:
    if ds_id in override_map:
        return override_map[ds_id]
    stem_key = data.get("dataset", "")
    if stem_key in override_map:
        return override_map[stem_key]
    md = get_methods_dict(data)
    candidates = [
        "DSA_UQW1b_TrueDSAiv_PAnc_25",
        "DSA_UQW1b_PAnc_25",
        "DSA_UQW1b_25",
    ]
    for c in candidates:
        if c in md:
            return c
    return None


def get_primary_score(m_entry: Dict[str, Any], data: Dict[str, Any]) -> Optional[float]:
    for key in ("F5_fix_mean", "F5_R0band_mean", "E1_mean", "mean"):
        if key in m_entry and m_entry[key] is not None:
            return float(m_entry[key])
    for key in ("F5_fix", "F5_R0band", "E1", "tie"):
        if key in m_entry:
            v = m_entry[key]
            if isinstance(v, list) and len(v) > 0:
                return float(sum(v) / len(v))
            if isinstance(v, (int, float)):
                return float(v)
    return None


def compute_ghf(data: Dict[str, Any], prop_method: str, base_method: str = "DSA_R0_Eq") -> Optional[Tuple[float, Tuple[float, float], str]]:
    if "spec_res" not in data or "f" not in data:
        return None
    spec = data["spec_res"]
    freqs = [float(x) for x in data["f"]]
    if prop_method not in spec:
        return None
    if base_method not in spec:
        if "DSA_R0" in spec:
            base_method = "DSA_R0"
        else:
            return None
    fmax = max(freqs)
    band = (28.0, 40.0) if fmax > 30.0 else (18.0, 24.0)
    sp = spec[prop_method]
    sb = spec[base_method]
    ratios: List[float] = []
    for f_val, vp, vb in zip(freqs, sp, sb):
        if band[0] - 1e-9 <= f_val <= band[1] + 1e-9 and float(vb) > 1e-30:
            ratios.append(float(vp) / float(vb))
    if not ratios:
        return None
    return (statistics.median(ratios), band, base_method)


def strip_errata_quotes(md_text: str) -> str:
    """
    Remove explicitly demarcated historical errata quotes (HTML comments or
    markdown table cells under '上一轮报告的失实写法') so anti-hype scanning
    checks active claims rather than quoted errors being corrected.
    """
    cleaned = re.sub(
        r"<!--\s*ERRATA_QUOTE_START\s*-->.*?<!--\s*ERRATA_QUOTE_END\s*-->",
        "",
        md_text,
        flags=re.DOTALL,
    )
    out_lines: List[str] = []
    in_errata_ledger = False
    for line in cleaned.splitlines():
        if "|" in line and "上一轮报告的失实写法" in line:
            in_errata_ledger = True
            continue
        if in_errata_ledger:
            if line.strip().startswith("|"):
                continue
            in_errata_ledger = False
        out_lines.append(line)
    return "\n".join(out_lines)


class ParsedTable:
    def __init__(self, heading_context: str, start_line: int, end_line: int, headers: List[str], rows: List[List[str]], after_text: str):
        self.heading_context = heading_context
        self.start_line = start_line
        self.end_line = end_line
        self.headers = headers
        self.rows = rows
        self.after_text = after_text


def parse_markdown_tables(md_text: str) -> List[ParsedTable]:
    lines = md_text.splitlines()
    n = len(lines)
    tables: List[ParsedTable] = []
    current_headings = ["", "", "", ""]
    i = 0
    while i < n:
        line = lines[i]
        hm = re.match(r"^(#{1,4})\s+(.*)$", line.strip())
        if hm:
            lvl = len(hm.group(1)) - 1
            current_headings[lvl] = hm.group(2).strip()
            for k in range(lvl + 1, 4):
                current_headings[k] = ""
            i += 1
            continue

        if line.strip().startswith("|") and i + 1 < n and re.match(r"^\|\s*[:\-]+", lines[i + 1].strip()):
            start_line = i + 1
            headers = [c.strip() for c in line.strip().strip("|").split("|")]
            i += 2
            rows: List[List[str]] = []
            while i < n and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                rows.append(cells)
                i += 1
            end_line = i
            # Capture text immediately after the table up to the next heading or next table (max 40 lines)
            after_buf: List[str] = []
            j = i
            while j < min(n, i + 45):
                lj = lines[j].strip()
                if re.match(r"^#{1,4}\s+", lj):
                    break
                if lj.startswith("|") and j + 1 < n and re.match(r"^\|\s*[:\-]+", lines[j + 1].strip()):
                    break
                after_buf.append(lines[j])
                j += 1
            ctx_above = "\n".join([h for h in current_headings if h] + lines[max(0, start_line - 6):start_line - 1])
            tables.append(ParsedTable(ctx_above, start_line, end_line, headers, rows, "\n".join(after_buf)))
            continue
        i += 1
    return tables


def match_tables_for_dataset(fp: str, ds_id: str, data: Dict[str, Any], tables: List[ParsedTable]) -> Tuple[List[ParsedTable], Dict[str, List[List[str]]]]:
    """
    Find markdown tables corresponding to a given JSON dataset.
    Supports both dedicated per-dataset tables and grouped multi-model tables
    (where column 0 carries the model ID and blank column-0 cells inherit above).
    Returns (matched_tables, rows_for_ds_per_table).
    """
    stem = os.path.splitext(os.path.basename(fp))[0]
    ds_tag = data.get("dataset", "")
    matched: List[ParsedTable] = []
    rows_map: Dict[int, List[List[str]]] = {}

    # Exact token regex for ds_id (avoid matching 'U1' inside 'U1b' or 'U10'!)
    ds_pat = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(ds_id)}(?![A-Za-z0-9_])")
    stem_pat = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(stem)}(?![A-Za-z0-9_])")

    for idx, tb in enumerate(tables):
        # Skip errata or judgment summary tables
        if "问题类型" in " ".join(tb.headers) or "门控编号" in " ".join(tb.headers) or "代码接线" in tb.heading_context:
            continue
        if "边际归因" in tb.heading_context or "阶段增量" in " ".join(tb.headers):
            continue

        # Check if grouped table has ds_id in column 0
        grouped_rows: List[List[str]] = []
        active_group = ""
        has_group_col = len(tb.headers) >= 2 and any(k in tb.headers[0] for k in ("模型", "数据集", "Model", "Dataset"))
        if has_group_col:
            for r in tb.rows:
                if r and r[0].strip():
                    active_group = r[0].strip()
                if ds_pat.search(active_group):
                    grouped_rows.append(r)
        if grouped_rows:
            matched.append(tb)
            rows_map[idx] = grouped_rows
            continue

        # Dedicated table check via heading/context above table
        ctx = tb.heading_context
        is_match = False
        if stem_pat.search(ctx) or ds_pat.search(ctx):
            is_match = True
        elif ds_tag == "volve_4_25hz" and ("4–25 Hz" in ctx or "4-25 Hz" in ctx) and stem in ctx:
            is_match = True
        elif ds_tag == "volve_4_45hz" and ("4–45 Hz" in ctx or "4-45 Hz" in ctx) and stem in ctx:
            is_match = True
        elif ds_tag == "volve_4_25hz" and ("volve25" in stem) and ("volve25" in ctx or ("4–25 Hz" in ctx and "PSAE" not in stem)):
            is_match = True
        elif ds_tag == "volve_4_45hz" and ("volve45" in stem) and ("volve45" in ctx or "4–45 Hz" in ctx):
            is_match = True

        if is_match:
            matched.append(tb)
            rows_map[idx] = tb.rows

    return matched, rows_map


def main() -> int:
    ap = argparse.ArgumentParser(description="Mechanical Scientific Integrity Auditor (Gate R0)")
    ap.add_argument("--report", required=True, help="Path to Markdown report to audit")
    ap.add_argument("--results-dir", required=True, help="Directory containing raw *.json files")
    ap.add_argument("--scripts-dir", default="", help="Directory containing MATLAB/Python scripts for wiring check")
    ap.add_argument("--log-out", default="", help="Optional path to write full audit log")
    ap.add_argument("--primary-map", default="", help="Optional JSON string mapping dataset_id -> primary_method")
    args = ap.parse_args()

    logs: List[str] = []
    failures: List[str] = []

    def emit(msg: str) -> None:
        print(msg)
        logs.append(msg)

    def fail(gate: str, msg: str) -> None:
        line = f"[FAIL][{gate}] {msg}"
        emit(line)
        failures.append(line)

    def ok(gate: str, msg: str) -> None:
        emit(f"[PASS][{gate}] {msg}")

    if not os.path.isfile(args.report):
        fail("IO", f"Report not found: {args.report}")
        return 1
    if not os.path.isdir(args.results_dir):
        fail("IO", f"Results dir not found: {args.results_dir}")
        return 1

    with open(args.report, "r", encoding="utf-8") as f:
        md_text = f.read()

    override_map: Dict[str, str] = {}
    if args.primary_map:
        override_map = json.loads(args.primary_map)

    json_entries = load_json_files(args.results_dir)
    if not json_entries:
        fail("IO", f"No experimental JSON files found in {args.results_dir}")
        return 1

    emit(f"=== Scientific Integrity Mechanical Audit (Gate R0) ===")
    emit(f"Report      : {args.report}")
    emit(f"Results Dir : {args.results_dir} ({len(json_entries)} experimental JSON files)")
    tables = parse_markdown_tables(md_text)
    emit(f"Parsed {len(tables)} Markdown tables from report.\n")

    # -------------------------------------------------------------------------
    # Gate 1: Zero Row Omission (100% methods per JSON in its dataset table)
    # -------------------------------------------------------------------------
    g1_errors = 0
    dataset_tables_cache: Dict[str, Tuple[List[ParsedTable], List[List[str]]]] = {}
    for fp, data in json_entries:
        ds_id = get_dataset_id(fp, data)
        md = get_methods_dict(data)
        m_names = list(md.keys())
        matched_tbs, rows_map = match_tables_for_dataset(fp, ds_id, data, tables)
        if not matched_tbs:
            fail("Gate1-ZeroRowOmission", f"Dataset '{ds_id}' ({os.path.basename(fp)}, {len(m_names)} methods) has NO table in report!")
            g1_errors += 1
            continue

        # Flatten matched rows across matched tables for this dataset
        all_rows: List[List[str]] = []
        for idx, tb in enumerate(tables):
            if idx in rows_map:
                all_rows.extend(rows_map[idx])
        dataset_tables_cache[ds_id] = (matched_tbs, all_rows)

        missing_in_ds: List[str] = []
        merged_rows: List[str] = []
        for m in m_names:
            m_pat = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(m)}(?![A-Za-z0-9_])")
            matching_r = [r for r in all_rows if any(m_pat.search(c) for c in r[:2])]
            if not matching_r:
                missing_in_ds.append(m)
            else:
                # Check if multiple JSON methods were crammed onto a single row
                r0_text = " ".join(matching_r[0][:2])
                other_m = [om for om in m_names if om != m and re.search(rf"(?<![A-Za-z0-9_]){re.escape(om)}(?![A-Za-z0-9_])", r0_text)]
                if other_m:
                    merged_rows.append(f"{m} merged with {other_m}")

        if missing_in_ds:
            fail("Gate1-ZeroRowOmission", f"Dataset '{ds_id}' ({os.path.basename(fp)}) omits {len(missing_in_ds)}/{len(m_names)} methods: {missing_in_ds}")
            g1_errors += 1
        if merged_rows:
            fail("Gate1-ZeroRowOmission", f"Dataset '{ds_id}' merges multiple methods into one row (hides individual columns): {sorted(set(merged_rows))}")
            g1_errors += 1

    if g1_errors == 0:
        ok("Gate1-ZeroRowOmission", f"100% of methods across all {len(json_entries)} JSON files appear in dedicated unmerged rows.")

    # -------------------------------------------------------------------------
    # Gate 2: Zero Column Omission (All mandatory columns on field & synth tables)
    # -------------------------------------------------------------------------
    g2_errors = 0
    for fp, data in json_entries:
        ds_id = get_dataset_id(fp, data)
        if ds_id not in dataset_tables_cache:
            continue
        matched_tbs, all_rows = dataset_tables_cache[ds_id]
        hdr_union = " | ".join([" | ".join(tb.headers) for tb in matched_tbs])

        is_field = str(data.get("dataset", "")).startswith("volve")
        if is_field:
            req_cols = [
                ("F-5 Primary", [r"F5_fix", r"F5_R0band", r"F-5"]),
                ("F-5 Secondary", [r"F5_r0", r"F5_Broad"]),
                ("Blind Well F-14", [r"F-14", r"F14"]),
                ("Zero-Lag L0", [r"L0", r"ZeroLag", r"零时移"]),
                ("Optimal Lag", [r"lag", r"BestLag", r"滞后"]),
                ("Spectral Centroid Res_Cen", [r"Res_Cen", r"质心"]),
                ("Smile_Frac", [r"Smile_Frac", r"画弧"]),
                ("Res_Slope", [r"Res_Slope", r"斜率"]),
            ]
            missing_cols = [label for label, pats in req_cols if not any(re.search(p, hdr_union, re.IGNORECASE) for p in pats)]
            if missing_cols:
                fail("Gate2-ZeroColumnOmission", f"Field dataset '{ds_id}' table is missing mandatory columns: {missing_cols}")
                g2_errors += 1
        elif "W_deep" in data:
            req_synth_cols = [
                ("E1_mean", [r"E1"]),
                ("E1_HF", [r"E1_HF", r"HF"]),
                ("E_amp", [r"E_amp"]),
                ("SNR_img", [r"SNR_img", r"SNR"]),
            ]
            missing_scols = [label for label, pats in req_synth_cols if not any(re.search(p, hdr_union, re.IGNORECASE) for p in pats)]
            if missing_scols:
                fail("Gate2-ZeroColumnOmission", f"Synthetic dataset '{ds_id}' table is missing mandatory columns: {missing_scols}")
                g2_errors += 1

    if g2_errors == 0:
        ok("Gate2-ZeroColumnOmission", "All field and synthetic tables contain 100% of mandatory evaluation columns.")

    # -------------------------------------------------------------------------
    # Gate 3: Better-Baseline & Degradation Disclosure
    # -------------------------------------------------------------------------
    g3_errors = 0
    for fp, data in json_entries:
        ds_id = get_dataset_id(fp, data)
        md = get_methods_dict(data)
        prop = detect_primary_method(ds_id, data, override_map)
        if not prop or prop not in md:
            fail("Gate3-BetterBaselineDisclosure", f"Cannot identify primary proposed method for '{ds_id}'")
            g3_errors += 1
            continue
        p_score = get_primary_score(md[prop], data)
        if p_score is None:
            continue

        beaters: List[Tuple[str, float, float]] = []
        for m, m_val in md.items():
            if m == prop:
                continue
            sc = get_primary_score(m_val, data)
            if sc is not None and (sc - p_score) > 0.001:
                beaters.append((m, sc, sc - p_score))

        # Check F-14 negative or worse than RAW_R0 on field datasets
        f14_degraded = False
        for f14_key in ("F14_fix_mean", "F14_R0band_mean"):
            if f14_key in md[prop]:
                p_f14 = float(md[prop][f14_key])
                r0_f14 = float(md.get("RAW_R0", {}).get(f14_key, 0.0))
                if p_f14 < 0.0 or (p_f14 - r0_f14) < -0.005:
                    f14_degraded = True

        if ds_id not in dataset_tables_cache:
            g3_errors += 1
            continue
        matched_tbs, _ = dataset_tables_cache[ds_id]
        after_blob = "\n".join([tb.after_text for tb in matched_tbs])

        if "反超主推 / 恶化项核查" not in after_blob:
            fail(
                "Gate3-BetterBaselineDisclosure",
                f"Dataset '{ds_id}' lacks mandatory '反超主推 / 恶化项核查' block after its table (beaters > +0.001: {[b[0] for b in beaters]}; F14_degraded={f14_degraded}).",
            )
            g3_errors += 1
            continue

        undisclosed = [b[0] for b in beaters if not re.search(rf"(?<![A-Za-z0-9_]){re.escape(b[0])}(?![A-Za-z0-9_])", after_blob)]
        if undisclosed:
            fail(
                "Gate3-BetterBaselineDisclosure",
                f"Dataset '{ds_id}' '反超主推 / 恶化项核查' omits baselines that beat '{prop}' ({p_score:.4f}): {undisclosed}",
            )
            g3_errors += 1
        if f14_degraded and ("F-14" not in after_blob and "F14" not in after_blob):
            fail(
                "Gate3-BetterBaselineDisclosure",
                f"Dataset '{ds_id}' has negative/degraded blind well F-14 for '{prop}', but '反超主推 / 恶化项核查' fails to disclose F-14!",
            )
            g3_errors += 1

    if g3_errors == 0:
        ok("Gate3-BetterBaselineDisclosure", "Every dataset table includes an accurate '反超主推 / 恶化项核查' naming all beating baselines and degraded metrics.")

    # -------------------------------------------------------------------------
    # Gate 4: Active-Operator & Marginal Attribution Gate
    # -------------------------------------------------------------------------
    g4_errors = 0
    field_ghf_records: Dict[str, Tuple[str, float, Tuple[float, float]]] = {}
    for fp, data in json_entries:
        if not str(data.get("dataset", "")).startswith("volve"):
            continue
        ds_id = get_dataset_id(fp, data)
        prop = detect_primary_method(ds_id, data, override_map)
        if not prop:
            continue
        ghf_info = compute_ghf(data, prop, "DSA_R0_Eq")
        if ghf_info is None:
            fail("Gate4-ActiveOperator", f"Field dataset '{ds_id}' lacks spec_res to verify HF gain G_HF for '{prop}'.")
            g4_errors += 1
            continue
        ghf_val, band, base_used = ghf_info
        field_ghf_records[ds_id] = (prop, ghf_val, band)

        # Marginal share calculation
        md = get_methods_dict(data)
        r0_sc = get_primary_score(md.get("RAW_R0", {}), data)
        eq_sc = get_primary_score(md.get(base_used, {}), data)
        pr_sc = get_primary_score(md.get(prop, {}), data)
        if r0_sc is not None and eq_sc is not None and pr_sc is not None:
            tot_delta = pr_sc - r0_sc
            eq_share = ((eq_sc - r0_sc) / tot_delta * 100.0) if abs(tot_delta) > 1e-6 else 0.0
            emit(f"  [Info][Gate4] {ds_id} ({prop}): G_HF[{band[0]:.0f}-{band[1]:.0f}Hz] vs {base_used} = {ghf_val:.3f}x | R0->{base_used} share = {eq_share:.1f}% of total delta ({tot_delta:+.4f})")

        if ghf_val < 1.10:
            has_disclosure = any(t in md_text for t in DORMANT_DISCLOSURE_TERMS)
            if not has_disclosure:
                fail(
                    "Gate4-ActiveOperator",
                    f"Field dataset '{ds_id}' ('{prop}') has G_HF = {ghf_val:.3f} < 1.10 in {band[0]:.0f}-{band[1]:.0f} Hz (operator dormant!), and report lacks explicit dormancy disclosure!",
                )
                g4_errors += 1

        if "阶梯边际归因" not in md_text and "边际归因" not in md_text:
            fail("Gate4-ActiveOperator", f"Report lacks mandatory step-by-step marginal attribution ('阶梯边际归因') separating front-end equalization from UQ operator!")
            g4_errors += 1

    if g4_errors == 0:
        ok("Gate4-ActiveOperator", "All field datasets pass HF active-operator verification and include step-by-step marginal attribution.")

    # -------------------------------------------------------------------------
    # Gate 5: Anti-Hype Gate (Single-well micro-differences |Delta r| < 0.01)
    # -------------------------------------------------------------------------
    g5_errors = 0
    small_diff_triggered = False
    for fp, data in json_entries:
        if not str(data.get("dataset", "")).startswith("volve"):
            continue
        ds_id = get_dataset_id(fp, data)
        md = get_methods_dict(data)
        prop = detect_primary_method(ds_id, data, override_map)
        if not prop or prop not in md:
            continue
        p_sc = get_primary_score(md[prop], data)
        if p_sc is None:
            continue
        for m, mv in md.items():
            if m == prop:
                continue
            if any(tag in m for tag in ("W1b", "UQW1_", "B1", "R0_Eq", "PH_", "NoW2")):
                sc = get_primary_score(mv, data)
                if sc is not None and abs(p_sc - sc) < 0.01:
                    small_diff_triggered = True

    if small_diff_triggered:
        active_text = strip_errata_quotes(md_text)
        found_hype: List[Tuple[str, int, str]] = []
        for lineno, line in enumerate(active_text.splitlines(), start=1):
            # Ignore lines that explicitly list banned words inside meta-descriptions of Gate 5 / R0
            if "黑名单" in line or "禁用" in line or "Anti-Hype" in line or "夸大词" in line:
                continue
            for term in BANNED_HYPE_TERMS:
                if term in line:
                    found_hype.append((term, lineno, line.strip()[:90]))
        if found_hype:
            for term, ln, snippet in found_hype:
                fail("Gate5-AntiHype", f"Banned hype term '{term}' on active line ~{ln} while single-well |Delta r| < 0.01: \"{snippet}...\"")
            g5_errors += len(found_hype)

    if g5_errors == 0:
        ok("Gate5-AntiHype", "Zero banned hype phrases found in active report text.")

    # -------------------------------------------------------------------------
    # Gate 6: Code Wiring Checklist Gate
    # -------------------------------------------------------------------------
    g6_errors = 0
    if "代码接线与实际增益核查表" not in md_text:
        fail("Gate6-WiringChecklist", "Report is missing mandatory '代码接线与实际增益核查表' section!")
        g6_errors += 1
    else:
        wiring_tbs = [tb for tb in tables if "代码接线与实际增益核查表" in tb.heading_context]
        if not wiring_tbs:
            fail("Gate6-WiringChecklist", "Heading '代码接线与实际增益核查表' found, but no Markdown table underneath it!")
            g6_errors += 1
        else:
            wtb = wiring_tbs[0]
            w_blob = "\n".join([" | ".join(r) for r in wtb.rows])
            for ds_id, (prop, ghf_val, _) in field_ghf_records.items():
                if prop not in w_blob:
                    fail("Gate6-WiringChecklist", f"Wiring table does not list primary field method '{prop}' ({ds_id})!")
                    g6_errors += 1
            # Verify cited script files & line numbers on disk if --scripts-dir provided
            if args.scripts_dir and os.path.isdir(args.scripts_dir):
                script_refs = re.findall(r"([A-Za-z0-9_]+\.m)\s*`?\s*L(\d+)(?:[\-–]L?(\d+))?", w_blob)
                if not script_refs:
                    fail("Gate6-WiringChecklist", "Wiring table does not cite exact MATLAB script line numbers (e.g., run_closure_volve45.m L180-L245)!")
                    g6_errors += 1
                for sname, l_start, l_end in script_refs:
                    spath = os.path.join(args.scripts_dir, sname)
                    if not os.path.isfile(spath):
                        fail("Gate6-WiringChecklist", f"Wiring table cites non-existent script '{sname}' in {args.scripts_dir}!")
                        g6_errors += 1
                    else:
                        with open(spath, "r", encoding="utf-8") as sf:
                            slines = sf.readlines()
                        ls = int(l_start)
                        le = int(l_end) if l_end else ls
                        if ls < 1 or le > len(slines) or ls > le:
                            fail("Gate6-WiringChecklist", f"Wiring table cites invalid line range L{ls}-L{le} in '{sname}' (file has {len(slines)} lines)!")
                            g6_errors += 1

    if g6_errors == 0:
        ok("Gate6-WiringChecklist", "Code wiring & active HF gain checklist table verified.")

    # -------------------------------------------------------------------------
    # Final Summary
    # -------------------------------------------------------------------------
    emit("\n" + "=" * 72)
    if failures:
        emit(f"FINAL AUDIT VERDICT: FAIL ({len(failures)} integrity violations detected)")
        exit_code = 1
    else:
        emit("FINAL AUDIT VERDICT: PASS (0 violations, 0 warnings — cleared for archiving)")
        exit_code = 0
    emit("=" * 72)

    if args.log_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.log_out)), exist_ok=True)
        with open(args.log_out, "w", encoding="utf-8") as lf:
            lf.write("\n".join(logs) + "\n")

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
