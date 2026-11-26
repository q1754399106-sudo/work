#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_full_tables.py
=====================
Mechanically renders 100% zero-omission Markdown tables, step-by-step marginal
attribution tables, and mandatory '反超主推 / 恶化项核查' disclosure blocks
directly from raw *.json files (eliminating manual copy-paste and row/column
cherry-picking).
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import statistics
from typing import Any, Dict, List, Optional, Tuple


def mean_or_val(v: Any) -> float:
    if isinstance(v, list):
        return float(sum(v) / len(v)) if v else float("nan")
    return float(v)


def idx_or_val(v: Any, idx: int) -> float:
    if isinstance(v, list):
        return float(v[min(idx, len(v) - 1)]) if v else float("nan")
    return float(v)


def fmt_m_vec(v: Any) -> str:
    if isinstance(v, list):
        return " / ".join(f"{float(x):.3f}" for x in v)
    return f"{float(v):.3f} (const)"


def compute_method_ghf(data: Dict[str, Any], method: str, base_method: str = "DSA_R0_Eq") -> Optional[float]:
    if "spec_res" not in data or "f" not in data:
        return None
    spec = data["spec_res"]
    freqs = [float(x) for x in data["f"]]
    if method not in spec:
        return None
    if base_method not in spec:
        base_method = "DSA_R0" if "DSA_R0" in spec else ""
    if not base_method:
        return None
    band = (28.0, 40.0) if max(freqs) > 30.0 else (18.0, 24.0)
    ratios = [
        float(vp) / float(vb)
        for f_val, vp, vb in zip(freqs, spec[method], spec[base_method])
        if band[0] - 1e-9 <= f_val <= band[1] + 1e-9 and float(vb) > 1e-30
    ]
    return statistics.median(ratios) if ratios else None


def render_field_table(fp: str, data: Dict[str, Any], prop: str) -> str:
    stem = os.path.splitext(os.path.basename(fp))[0]
    ds = data.get("dataset", stem)
    is_45 = "45" in ds or max(data.get("f", [25])) > 30
    p_key = "F5_fix" if is_45 else "F5_R0band"
    s_key = "F5_r0" if is_45 else "F5_Broad"
    f14_key = "F14_fix" if is_45 else "F14_R0band"
    hf_label = "G_HF[28-40Hz]" if is_45 else "G_HF[18-24Hz]"
    ms = data.get("m", [])
    m_hdr = "/".join(str(x) for x in ms)
    m1_idx = ms.index(1.0) if 1.0 in ms else (ms.index(1) if 1 in ms else len(ms) // 2)

    md = data["methods"]
    slopes = data.get("res_slope", {})

    lines: List[str] = []
    lines.append(f"### 数据集 `{stem}`（`{ds}`，全 {len(md)} 种方法零省略主表）")
    lines.append("")
    lines.append(
        f"| 方法 | `{p_key}` 均值 | `{p_key}` (`m={m_hdr}`) | `{s_key}` 均值 | `F-5 L0` (`m=1`) | `F-5 lag` (`m=1`, m) | 盲井 `F-14` 均值 (`m=1`) | `F-14 L0 / lag` (`m=1`) | `Res_Cen` (cyc/km) | `Smile_Frac` | `Res_Slope` (/Hz) | `{hf_label}` |"
    )
    lines.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")

    p_score = mean_or_val(md[prop][f"{p_key}_mean"] if f"{p_key}_mean" in md[prop] else md[prop][p_key])
    beaters: List[Tuple[str, float, float]] = []
    degradations: List[str] = []

    for m, mv in md.items():
        pm = float(mv.get(f"{p_key}_mean", mean_or_val(mv[p_key])))
        sm = float(mv.get(f"{s_key}_mean", mean_or_val(mv[s_key])))
        l0 = idx_or_val(mv["F5_ZeroLag"], m1_idx)
        lag = idx_or_val(mv["F5_BestLag_m"], m1_idx)
        f14m = float(mv.get(f"{f14_key}_mean", mean_or_val(mv[f14_key])))
        f14_m1 = idx_or_val(mv[f14_key], m1_idx)
        f14_l0 = idx_or_val(mv["F14_ZeroLag"], m1_idx) if "F14_ZeroLag" in mv else float("nan")
        f14_lag = idx_or_val(mv["F14_BestLag_m"], m1_idx) if "F14_BestLag_m" in mv else float("nan")
        f14_l0_str = f"`{f14_l0:+.3f}` / `{f14_lag:+.0f}m`" if not (f14_l0 != f14_l0) else "—"
        cen = idx_or_val(mv["Res_Cen"], m1_idx)
        smf = idx_or_val(mv["Smile_Frac"], m1_idx)
        sl = float(slopes.get(m, float("nan")))
        ghf = compute_method_ghf(data, m, "DSA_R0_Eq")
        ghf_str = f"`{ghf:.3f}x`" if ghf is not None else "—"

        m_disp = f"**`{m}`（主推）**" if m == prop else f"`{m}`"
        lines.append(
            f"| {m_disp} | `{pm:.4f}` | `{fmt_m_vec(mv[p_key])}` | `{sm:.4f}` | `{l0:+.4f}` | `{lag:+.1f}` | `{f14m:+.4f}` (`{f14_m1:+.4f}`) | {f14_l0_str} | `{cen:.2f}` | `{smf*100:.1f}%` | `{sl:+.4f}` | {ghf_str} |"
        )
        if m != prop and (pm - p_score) > 0.001:
            beaters.append((m, pm, pm - p_score))

    # Check degradations on prop
    r0_f14 = float(md["RAW_R0"].get(f"{f14_key}_mean", mean_or_val(md["RAW_R0"][f14_key]))) if "RAW_R0" in md else 0.0
    pr_f14 = float(md[prop].get(f"{f14_key}_mean", mean_or_val(md[prop][f14_key])))
    if pr_f14 < 0.0 or pr_f14 < r0_f14:
        degradations.append(f"盲井 `F-14` 主推 `{prop}` 为 `{pr_f14:+.4f}`（对比 `RAW_R0 = {r0_f14:+.4f}`，差值 `{pr_f14 - r0_f14:+.4f}`）")
    r0_l0 = idx_or_val(md["RAW_R0"]["F5_ZeroLag"], 0) if "RAW_R0" in md else 0.0
    pr_l0 = idx_or_val(md[prop]["F5_ZeroLag"], m1_idx)
    if pr_l0 < r0_l0:
        degradations.append(f"`F-5` 零时移 `L0` 主推 `{prop}` 为 `{pr_l0:+.4f}`（对比 `RAW_R0 = {r0_l0:+.4f}`）")

    lines.append("")
    if beaters:
        b_str = "；".join([f"`{b[0]} = {b[1]:.4f}`（较主推 `{b[2]:+.4f}`）" for b in sorted(beaters, key=lambda x: -x[1])])
    else:
        b_str = f"本表内无任何基线在 `{p_key}` 均值上超出主推 `{prop}` 达 `+0.001` 以上"
    d_str = "；".join(degradations) if degradations else "无盲井 `F-14` 变负或零时移 `L0` 低于 `RAW_R0` 现象"
    lines.append(f"- **反超主推 / 恶化项核查**：① 高于主推 `{prop}`（`{p_score:.4f}`）> `+0.001` 的方法：{b_str}；② 恶化/负值项核查：{d_str}。")
    lines.append("")

    # Marginal attribution sub-table if RAW_R0, DSA_R0, DSA_R0_Eq exist
    if all(k in md for k in ("RAW_R0", "DSA_R0", "DSA_R0_Eq", prop)):
        s0 = float(md["RAW_R0"].get(f"{p_key}_mean", mean_or_val(md["RAW_R0"][p_key])))
        s1 = float(md["DSA_R0"].get(f"{p_key}_mean", mean_or_val(md["DSA_R0"][p_key])))
        s2 = float(md["DSA_R0_Eq"].get(f"{p_key}_mean", mean_or_val(md["DSA_R0_Eq"][p_key])))
        s3 = p_score
        tot = s3 - s0
        pct = lambda d: f"{(d / tot * 100.0):+.1f}%" if abs(tot) > 1e-6 else "N/A"
        ghf_p = compute_method_ghf(data, prop, "DSA_R0_Eq")
        lines.append(f"#### `{stem}` 阶梯边际归因表（`RAW_R0 → DSA_R0 → DSA_R0_Eq → {prop}`）")
        lines.append("")
        lines.append(f"| 阶梯阶段 | 方法 | `{p_key}` 均值 | 阶段边际增量 `Δr` | 占总增量 (`{tot:+.4f}`) 比例 | `{hf_label}` |")
        lines.append("|---|---|---:|---:|---:|---:|")
        lines.append(f"| Step 0（无补偿原图） | `RAW_R0` | `{s0:.4f}` | `0.0000` | `0.0%` | — |")
        lines.append(f"| Step 1（+DSA 镜面防假频门控） | `DSA_R0` | `{s1:.4f}` | `{s1-s0:+.4f}` | `{pct(s1-s0)}` | — |")
        lines.append(f"| Step 2（+PSAE 前置谱锚定） | `DSA_R0_Eq` | `{s2:.4f}` | `{s2-s1:+.4f}` | `{pct(s2-s1)}` | `1.000x` |")
        lines.append(f"| Step 3（+同域 UQW1b 宽频 Q 补偿） | **`{prop}`** | **`{s3:.4f}`** | **`{s3-s2:+.4f}`** | **`{pct(s3-s2)}`** | **`{ghf_p:.3f}x`** |")
        lines.append("")

    return "\n".join(lines)


def render_synth_wdeep_table(fp: str, data: Dict[str, Any], prop: str) -> str:
    ds_id = data.get("model", os.path.splitext(os.path.basename(fp))[0])
    md = data["W_deep"]
    ms = data.get("m", [])
    m_hdr = "/".join(str(x) for x in ms)
    na = data.get("noise_acc", {})
    dip = data.get("dip_stats", {})

    lines: List[str] = []
    lines.append(
        f"### 合成模型 `{ds_id}`（全 {len(md)} 种方法零省略表，`RAW_R0 SNR_img = {float(md['RAW_R0']['SNR_img_dB_mean']):+.2f} dB`，盲倾角中位数 `{dip.get('med_abs_phi_deg', 0):.2f}°` / `p90={dip.get('p90_abs_phi_deg', 0):.2f}°`，估噪 `DSA B300 log2err={na.get('DSA_B300_log2err', 0):+.2f}` vs `B1={na.get('DSA_B1_log2err', 0):+.2f}`）"
    )
    lines.append("")
    lines.append(f"| 方法 | `E1` 均值 | `E1` (`m={m_hdr}`) | 高频 `E1_HF` 均值 | 振幅误差 `E_amp` 均值 (`m=1`) | `SNR_img` 均值 (dB) |")
    lines.append("|---|---:|---|---:|---:|---:|")

    p_score = float(md[prop]["E1_mean"])
    p_eamp = float(md[prop]["E_amp_mean"])
    beaters: List[Tuple[str, float, float]] = []
    eamp_beaters: List[Tuple[str, float]] = []
    m1_idx = ms.index(1.0) if 1.0 in ms else (ms.index(1) if 1 in ms else len(ms) // 2)

    for m, mv in md.items():
        e1m = float(mv["E1_mean"])
        e1h = float(mv["E1_HF_mean"])
        eam = float(mv["E_amp_mean"])
        ea1 = idx_or_val(mv["E_amp"], m1_idx)
        snr = float(mv["SNR_img_dB_mean"])
        m_disp = f"**`{m}`（主推）**" if m == prop else f"`{m}`"
        lines.append(f"| {m_disp} | `{e1m:.4f}` | `{fmt_m_vec(mv['E1'])}` | `{e1h:.4f}` | `{eam:.4f}` (`{ea1:.4f}`) | `{snr:+.2f}` |")
        if m != prop and (e1m - p_score) > 0.001:
            beaters.append((m, e1m, e1m - p_score))
        if m != prop and (p_eamp - eam) > 0.005:
            eamp_beaters.append((m, eam))

    lines.append("")
    if beaters:
        b_str = "；".join([f"`{b[0]} = {b[1]:.4f}`（较主推 `{b[2]:+.4f}`）" for b in sorted(beaters, key=lambda x: -x[1])])
    else:
        b_str = f"本模型无基线在 `E1` 均值上超出主推 `{prop}` 达 `+0.001` 以上"
    ea_str = (
        "；".join([f"`{x[0]} E_amp={x[1]:.4f}`" for x in sorted(eamp_beaters, key=lambda t: t[1])])
        if eamp_beaters
        else f"无基线 `E_amp` 低于主推 `{prop}` (`{p_eamp:.4f}`) 超过 `0.005`"
    )
    lines.append(f"- **反超主推 / 恶化项核查**：① `E1` 均值高于主推 `{prop}`（`{p_score:.4f}`）> `+0.001` 的方法：{b_str}；② `E_amp` 振幅误差优于主推（`{p_eamp:.4f}`）的方法：{ea_str}。")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description="Render 100% full Markdown tables from raw JSON files")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.results_dir, "*.json")))
    chunks: List[str] = []
    for fp in files:
        base = os.path.basename(fp)
        if base.endswith("_smoke.json") or "verdict" in base.lower() or "diagnostics" in base.lower():
            continue
        with open(fp, "r", encoding="utf-8") as f:
            data = json.load(f)
        if str(data.get("dataset", "")).startswith("volve"):
            md = data["methods"]
            prop = "DSA_UQW1b_TrueDSAiv_PAnc_25" if "DSA_UQW1b_TrueDSAiv_PAnc_25" in md else "DSA_UQW1b_PAnc_25"
            chunks.append(render_field_table(fp, data, prop))
        elif "W_deep" in data:
            chunks.append(render_synth_wdeep_table(fp, data, "DSA_UQW1b_25"))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as out_f:
        out_f.write("\n\n".join(chunks) + "\n")
    print(f"Rendered {len(chunks)} full tables to {args.out}")


if __name__ == "__main__":
    main()
