#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_report_integrity.py
=========================
Mechanical Scientific Integrity & Anti-Selective-Reporting Auditor (Gate X0).
Enforces all Six Zero-Tolerance Iron Laws against raw *.json outputs and MATLAB scripts.
Exit code: 0 if and only if 0 failures and 0 warnings; 1 otherwise.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import statistics
import sys
from typing import Any, Dict, List, Optional, Tuple

BANNED_HYPE_TERMS = [
    "全面击败", "完胜", "碾压", "独占第一", "绝对第一",
    "显著超越", "全面超越", "全面领先", "大胜",
]

DORMANT_DISCLOSURE_TERMS = [
    "高频补偿处于休眠", "弱作用状态", "休眠态", "休眠状态", "关停状态", "被混叠 iv3 关停",
]


def load_json_files(results_dir: str) -> List[Tuple[str, Dict[str, Any]]]:
    files = sorted(glob.glob(os.path.join(results_dir, "*.json")))
    out: List[Tuple[str, Dict[str, Any]]] = []
    for fp in files:
        base = os.path.basename(fp)
        if base.endswith("_smoke.json"):
            continue
        if any(k in base.lower() for k in ("verdict", "diagnostics", "physics_verify")):
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
        "DSA_Bent_UQW1b_PAnc_25",
        "DSA_Bent_UQW1b_25",
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


def compute_ghf(data: Dict[str, Any], prop_method: str) -> Optional[Tuple[float, Tuple[float, float], str]]:
    if "spec_res" not in data or "f" not in data:
        return None
    spec = data["spec_res"]
    freqs = [float(x) for x in data["f"]]
    if prop_method not in spec:
        return None
    base_method = ""
    for cand in ("DSA_Bent_R0_Eq", "DSA_R0_Eq", "DSA_Str_R0_Eq", "DSA_Bent_R0", "DSA_R0"):
        if cand in spec:
            base_method = cand
            break
    if not base_method:
        return None
    fmax = max(freqs)
    band = (28.0, 40.0) if fmax > 30.0 else (18.0, 24.0)
    sp = spec[prop_method]
    sb = spec[base_method]
    ratios = [
        float(vp) / float(vb)
        for f_val, vp, vb in zip(freqs, sp, sb)
        if band[0] - 1e-9 <= f_val <= band[1] + 1e-9 and float(vb) > 1e-30
    ]
    if not ratios:
        return None
    return (statistics.median(ratios), band, base_method)


def strip_errata_quotes(md_text: str) -> str:
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


def match_tables_for_dataset(fp: str, ds_id: str, data: Dict[str, Any], tables: List[ParsedTable]) -> Tuple[List[ParsedTable], Dict[int, List[List[str]]]]:
    stem = os.path.splitext(os.path.basename(fp))[0]
    matched: List[ParsedTable] = []
    rows_map: Dict[int, List[List[str]]] = {}
    ds_pat = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(ds_id)}(?![A-Za-z0-9_])")
    stem_pat = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(stem)}(?![A-Za-z0-9_])")
    for idx, tb in enumerate(tables):
        hdr_s = " ".join(tb.headers)
        if any(skip in hdr_s for skip in ("问题类型", "门控编号", "阶段边际增量")) or "代码接线" in tb.heading_context or "边际归因" in tb.heading_context:
            continue
        if "成本" in tb.heading_context or "LU" in hdr_s:
            continue
        ctx = tb.heading_context
        if stem_pat.search(ctx) or ds_pat.search(ctx):
            matched.append(tb)
            rows_map[idx] = tb.rows
    return matched, rows_map


def main() -> int:
    ap = argparse.ArgumentParser(description="Mechanical Scientific Integrity Auditor (Gate X0)")
    ap.add_argument("--report", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--scripts-dir", default="")
    ap.add_argument("--log-out", default="")
    ap.add_argument("--primary-map", default="")
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

    if not os.path.isfile(args.report) or not os.path.isdir(args.results_dir):
        fail("IO", f"Missing report or results dir: {args.report} / {args.results_dir}")
        return 1

    with open(args.report, "r", encoding="utf-8") as f:
        md_text = f.read()

    override_map: Dict[str, str] = json.loads(args.primary_map) if args.primary_map else {}
    json_entries = load_json_files(args.results_dir)
    if not json_entries:
        fail("IO", f"No experimental JSON files found in {args.results_dir}")
        return 1

    emit("=== Scientific Integrity Mechanical Audit (Gate X0) ===")
    emit(f"Report      : {args.report}")
    emit(f"Results Dir : {args.results_dir} ({len(json_entries)} experimental JSON files)")
    tables = parse_markdown_tables(md_text)
    emit(f"Parsed {len(tables)} Markdown tables from report.\n")

    # Gate 1: Zero Row Omission
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
                r0_text = " ".join(matching_r[0][:2])
                other_m = [om for om in m_names if om != m and re.search(rf"(?<![A-Za-z0-9_]){re.escape(om)}(?![A-Za-z0-9_])", r0_text)]
                if other_m:
                    merged_rows.append(f"{m} merged with {other_m}")
        if missing_in_ds:
            fail("Gate1-ZeroRowOmission", f"Dataset '{ds_id}' omits {len(missing_in_ds)}/{len(m_names)} methods: {missing_in_ds}")
            g1_errors += 1
        if merged_rows:
            fail("Gate1-ZeroRowOmission", f"Dataset '{ds_id}' merges methods: {sorted(set(merged_rows))}")
            g1_errors += 1
    if g1_errors == 0:
        ok("Gate1-ZeroRowOmission", f"100% of methods across all {len(json_entries)} JSON files appear in dedicated unmerged rows.")

    # Gate 2: Zero Column Omission
    g2_errors = 0
    for fp, data in json_entries:
        ds_id = get_dataset_id(fp, data)
        if ds_id not in dataset_tables_cache:
            continue
        matched_tbs, _ = dataset_tables_cache[ds_id]
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
                fail("Gate2-ZeroColumnOmission", f"Field dataset '{ds_id}' missing columns: {missing_cols}")
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
                fail("Gate2-ZeroColumnOmission", f"Synthetic dataset '{ds_id}' missing columns: {missing_scols}")
                g2_errors += 1
    if g2_errors == 0:
        ok("Gate2-ZeroColumnOmission", "All field and synthetic tables contain 100% of mandatory evaluation columns.")

    # Gate 3: Better-Baseline & Degradation Disclosure
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
        beaters = []
        for m, mv in md.items():
            if m == prop:
                continue
            sc = get_primary_score(mv, data)
            if sc is not None and (sc - p_score) > 0.001:
                beaters.append((m, sc, sc - p_score))
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
            fail("Gate3-BetterBaselineDisclosure", f"Dataset '{ds_id}' lacks '反超主推 / 恶化项核查' block.")
            g3_errors += 1
            continue
        undisclosed = [b[0] for b in beaters if not re.search(rf"(?<![A-Za-z0-9_]){re.escape(b[0])}(?![A-Za-z0-9_])", after_blob)]
        if undisclosed:
            fail("Gate3-BetterBaselineDisclosure", f"Dataset '{ds_id}' omits beating baselines: {undisclosed}")
            g3_errors += 1
        if f14_degraded and ("F-14" not in after_blob and "F14" not in after_blob):
            fail("Gate3-BetterBaselineDisclosure", f"Dataset '{ds_id}' omits degraded F-14 disclosure.")
            g3_errors += 1
    if g3_errors == 0:
        ok("Gate3-BetterBaselineDisclosure", "Every dataset table includes an accurate '反超主推 / 恶化项核查' naming all beating baselines and degraded metrics.")

    # Gate 4: Active-Operator & Marginal Attribution Gate
    g4_errors = 0
    field_ghf_records: Dict[str, Tuple[str, float, Tuple[float, float]]] = {}
    for fp, data in json_entries:
        if not str(data.get("dataset", "")).startswith("volve"):
            continue
        ds_id = get_dataset_id(fp, data)
        prop = detect_primary_method(ds_id, data, override_map)
        if not prop:
            continue
        ghf_info = compute_ghf(data, prop)
        if ghf_info is None:
            fail("Gate4-ActiveOperator", f"Field dataset '{ds_id}' lacks spec_res for G_HF.")
            g4_errors += 1
            continue
        ghf_val, band, base_used = ghf_info
        field_ghf_records[ds_id] = (prop, ghf_val, band)
        md = get_methods_dict(data)
        r0_sc = get_primary_score(md.get("RAW_R0", {}), data)
        eq_sc = get_primary_score(md.get(base_used, {}), data)
        pr_sc = get_primary_score(md.get(prop, {}), data)
        if r0_sc is not None and eq_sc is not None and pr_sc is not None:
            tot_delta = pr_sc - r0_sc
            eq_share = ((eq_sc - r0_sc) / tot_delta * 100.0) if abs(tot_delta) > 1e-6 else 0.0
            emit(f"  [Info][Gate4] {ds_id} ({prop}): G_HF[{band[0]:.0f}-{band[1]:.0f}Hz] vs {base_used} = {ghf_val:.3f}x | R0->{base_used} share = {eq_share:.1f}% of total delta ({tot_delta:+.4f})")
        if ghf_val < 1.10 and not any(t in md_text for t in DORMANT_DISCLOSURE_TERMS):
            fail("Gate4-ActiveOperator", f"Field dataset '{ds_id}' has G_HF={ghf_val:.3f} < 1.10 without dormancy disclosure.")
            g4_errors += 1
        if "阶梯边际归因" not in md_text and "边际归因" not in md_text:
            fail("Gate4-ActiveOperator", "Report lacks mandatory '阶梯边际归因'.")
            g4_errors += 1
    if g4_errors == 0:
        ok("Gate4-ActiveOperator", "All field datasets pass HF active-operator verification and include step-by-step marginal attribution.")

    # Gate 5: Anti-Hype Gate
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
            if m != prop:
                sc = get_primary_score(mv, data)
                if sc is not None and abs(p_sc - sc) < 0.01:
                    small_diff_triggered = True
    if small_diff_triggered:
        active_text = strip_errata_quotes(md_text)
        for lineno, line in enumerate(active_text.splitlines(), start=1):
            if any(k in line for k in ("黑名单", "禁用", "Anti-Hype", "夸大词")):
                continue
            for term in BANNED_HYPE_TERMS:
                if term in line:
                    fail("Gate5-AntiHype", f"Banned hype term '{term}' on line ~{lineno}: \"{line.strip()[:90]}...\"")
                    g5_errors += 1
    if g5_errors == 0:
        ok("Gate5-AntiHype", "Zero banned hype phrases found in active report text.")

    # Gate 6: Code Wiring Checklist Gate
    g6_errors = 0
    if "代码接线与实际增益核查表" not in md_text:
        fail("Gate6-WiringChecklist", "Report is missing mandatory '代码接线与实际增益核查表' section!")
        g6_errors += 1
    else:
        wiring_tbs = [tb for tb in tables if "代码接线与实际增益核查表" in tb.heading_context]
        if not wiring_tbs:
            fail("Gate6-WiringChecklist", "Heading '代码接线与实际增益核查表' has no Markdown table!")
            g6_errors += 1
        else:
            w_blob = "\n".join([" | ".join(r) for r in wiring_tbs[0].rows])
            for ds_id, (prop, _, _) in field_ghf_records.items():
                if prop not in w_blob:
                    fail("Gate6-WiringChecklist", f"Wiring table does not list primary field method '{prop}' ({ds_id})!")
                    g6_errors += 1
            if args.scripts_dir and os.path.isdir(args.scripts_dir):
                script_refs = re.findall(r"([A-Za-z0-9_]+\.m)\s*`?\s*L(\d+)(?:[\-–]L?(\d+))?", w_blob)
                if not script_refs:
                    fail("Gate6-WiringChecklist", "Wiring table does not cite exact MATLAB script line numbers!")
                    g6_errors += 1
                for sname, l_start, l_end in script_refs:
                    spath = os.path.join(args.scripts_dir, sname)
                    if not os.path.isfile(spath):
                        fail("Gate6-WiringChecklist", f"Wiring table cites non-existent script '{sname}'!")
                        g6_errors += 1
                    else:
                        slines = open(spath, "r", encoding="utf-8").readlines()
                        ls, le = int(l_start), int(l_end) if l_end else int(l_start)
                        if ls < 1 or le > len(slines) or ls > le:
                            fail("Gate6-WiringChecklist", f"Invalid line range L{ls}-L{le} in '{sname}' (len={len(slines)})!")
                            g6_errors += 1
    if g6_errors == 0:
        ok("Gate6-WiringChecklist", "Code wiring & active HF gain checklist table verified.")

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
