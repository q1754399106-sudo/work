#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_full_tables.py
=====================
Mechanically renders 100% zero-omission Markdown tables (Noise-Free m=1 Oracle +
Noisy All-m Synthetic Tables S1-S5, Volve 4-25 Hz 20-Method Full Table, Marginal
Attribution, and '反超主推 / 恶化项核查' disclosure blocks) directly from raw *.json.
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


def compute_method_ghf(data: Dict[str, Any], method: str) -> Optional[float]:
    if "spec_res" not in data or "f" not in data:
        return None
    spec = data["spec_res"]
    freqs = [float(x) for x in data["f"]]
    if method not in spec:
        return None
    base = next((c for c in ("DSA_Bent_R0_Eq", "DSA_R0_Eq", "DSA_Str_R0_Eq") if c in spec), "")
    if not base:
        return None
    band = (28.0, 40.0) if max(freqs) > 30.0 else (18.0, 24.0)
    ratios = [
        float(vp) / float(vb)
        for fv, vp, vb in zip(freqs, spec[method], spec[base])
        if band[0] - 1e-9 <= fv <= band[1] + 1e-9 and float(vb) > 1e-30
    ]
    return statistics.median(ratios) if ratios else None


def render_synth_table(fp: str, data: Dict[str, Any], prop: str = "DSA_Bent_UQW1b_25") -> str:
    ds_id = data.get("model", os.path.splitext(os.path.basename(fp))[0])
    md = data["W_deep"]
    nf = data.get("NoiseFree_m1", {})
    ms = data.get("m", [])
    m_hdr = "/".join(str(x) for x in ms)
    m1_idx = ms.index(1.0) if 1.0 in ms else (ms.index(1) if 1 in ms else len(ms) // 2)
    na = data.get("noise_acc", {})
    dip = data.get("dip_stats", {})
    xm = data.get("xmir_stats", {})

    lines: List[str] = []
    lines.append(
        f"### 合成模型 `{ds_id}`（全 {len(md)} 种方法零省略表，`RAW_R0 SNR_img = {float(md['RAW_R0']['SNR_img_dB_mean']):+.2f} dB`，"
        f"盲倾角中位 `{dip.get('med_abs_phi_deg', 0):.2f}°` / `p90={dip.get('p90_abs_phi_deg', 0):.2f}°`，"
        f"`|xmir_bent-xmir_str|` 中位 `{xm.get('W_deep_med_dxmir_m', 0):.1f} m` / `p90={xm.get('W_deep_p90_dxmir_m', 0):.1f} m`（门控翻转 `{xm.get('W_deep_gate_flip_pct', 0):.1f}%`），"
        f"估噪 `B300={na.get('DSA_B300_log2err', 0):+.2f}` vs `B1={na.get('DSA_B1_log2err', 0):+.2f}`）"
    )
    lines.append("")
    lines.append(
        f"| 方法 | 无噪 `m=1` (`E1` / `E_amp`) | 含噪 `E1` 均值 | 含噪 `E1` (`m={m_hdr}`) | 高频 `E1_HF` 均值 | 振幅误差 `E_amp` 均值 (`m=1`) | `SNR_img` 均值 (dB) | `Smile_Frac` |"
    )
    lines.append("|---|---:|---:|---|---:|---:|---:|---:|")

    p_score = float(md[prop]["E1_mean"])
    p_eamp = float(md[prop]["E_amp_mean"])
    beaters: List[Tuple[str, float, float]] = []
    eamp_beaters: List[Tuple[str, float]] = []
    nf_beaters: List[Tuple[str, float]] = []
    p_nf = float(nf[prop]["E1"]) if prop in nf else float("nan")

    for m, mv in md.items():
        e1m = float(mv["E1_mean"])
        e1h = float(mv["E1_HF_mean"])
        eam = float(mv["E_amp_mean"])
        ea1 = idx_or_val(mv["E_amp"], m1_idx)
        snr = float(mv["SNR_img_dB_mean"])
        smf = float(mv.get("Smile_Frac_mean", 0.0))
        if m in nf:
            nf_str = f"`{float(nf[m]['E1']):.4f}` / `{float(nf[m]['E_amp']):.4f}`"
            if m != prop and (float(nf[m]["E1"]) - p_nf) > 0.001:
                nf_beaters.append((m, float(nf[m]["E1"])))
        else:
            nf_str = "—"
        m_disp = f"**`{m}`（主推）**" if m == prop else f"`{m}`"
        lines.append(
            f"| {m_disp} | {nf_str} | `{e1m:.4f}` | `{fmt_m_vec(mv['E1'])}` | `{e1h:.4f}` | `{eam:.4f}` (`{ea1:.4f}`) | `{snr:+.2f}` | `{smf*100:.1f}%` |"
        )
        if m != prop and (e1m - p_score) > 0.001:
            beaters.append((m, e1m, e1m - p_score))
        if m != prop and (p_eamp - eam) > 0.005:
            eamp_beaters.append((m, eam))

    lines.append("")
    b_str = "；".join([f"`{b[0]} = {b[1]:.4f}`（较主推 `{b[2]:+.4f}`）" for b in sorted(beaters, key=lambda x: -x[1])]) if beaters else f"本模型无基线在含噪 `E1` 均值上超出主推 `{prop}` 达 `+0.001` 以上"
    ea_str = "；".join([f"`{x[0]} E_amp={x[1]:.4f}`" for x in sorted(eamp_beaters, key=lambda t: t[1])]) if eamp_beaters else f"无基线 `E_amp` 低于主推 `{prop}` (`{p_eamp:.4f}`) 超过 `0.005`"
    nf_str_b = "；".join([f"`{x[0]}={x[1]:.4f}`" for x in sorted(nf_beaters, key=lambda t: -t[1])]) if nf_beaters else "无"
    lines.append(
        f"- **反超主推 / 恶化项核查**：① 含噪 `E1` 均值高于主推 `{prop}`（`{p_score:.4f}`）> `+0.001` 的方法：{b_str}；"
        f"② `E_amp` 振幅误差优于主推（`{p_eamp:.4f}`）的方法：{ea_str}；"
        f"③ 无噪 `m=1` 理想真 $Q$ 基准高于主推（`{p_nf:.4f}`）的方法：{nf_str_b}。"
    )
    lines.append("")
    return "\n".join(lines)


def render_field_table(fp: str, data: Dict[str, Any], prop: str = "DSA_Bent_UQW1b_PAnc_25") -> str:
    stem = os.path.splitext(os.path.basename(fp))[0]
    ds = data.get("dataset", stem)
    ms = data.get("m", [])
    m_hdr = "/".join(str(x) for x in ms)
    m1_idx = ms.index(1.0) if 1.0 in ms else (ms.index(1) if 1 in ms else len(ms) // 2)
    md = data["methods"]
    slopes = data.get("res_slope", {})

    lines: List[str] = []
    lines.append(f"### 数据集 `{stem}`（`{ds}`，全 {len(md)} 种方法零省略主表）")
    lines.append("")
    lines.append(
        f"| 方法 | `F5_R0band` 均值 | `F5_R0band` (`m={m_hdr}`) | `F5_Broad` 均值 | `F-5 L0` (`m=1`) | `F-5 lag` (`m=1`, m) | 盲井 `F-14` 均值 (`m=1`) | `F-14 L0 / lag` (`m=1`) | `Res_Cen` (cyc/km) | `Smile_Frac` | `Res_Slope` (/Hz) | `G_HF[18-24Hz]` |"
    )
    lines.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")

    p_score = float(md[prop]["F5_R0band_mean"])
    beaters: List[Tuple[str, float, float]] = []
    degradations: List[str] = []

    for m, mv in md.items():
        pm = float(mv["F5_R0band_mean"])
        sm = float(mv["F5_Broad_mean"])
        l0 = idx_or_val(mv["F5_ZeroLag"], m1_idx)
        lag = idx_or_val(mv["F5_BestLag_m"], m1_idx)
        f14m = float(mv["F14_R0band_mean"])
        f14_m1 = idx_or_val(mv["F14_R0band"], m1_idx)
        f14_l0 = idx_or_val(mv["F14_ZeroLag"], m1_idx)
        f14_lag = idx_or_val(mv["F14_BestLag_m"], m1_idx)
        cen = idx_or_val(mv["Res_Cen"], m1_idx)
        smf = idx_or_val(mv["Smile_Frac"], m1_idx)
        sl = float(slopes.get(m, float("nan")))
        ghf = compute_method_ghf(data, m)
        ghf_str = f"`{ghf:.3f}x`" if ghf is not None else "—"
        m_disp = f"**`{m}`（主推）**" if m == prop else f"`{m}`"
        lines.append(
            f"| {m_disp} | `{pm:.4f}` | `{fmt_m_vec(mv['F5_R0band'])}` | `{sm:.4f}` | `{l0:+.4f}` | `{lag:+.1f}` | `{f14m:+.4f}` (`{f14_m1:+.4f}`) | `{f14_l0:+.3f}` / `{f14_lag:+.0f}m` | `{cen:.2f}` | `{smf*100:.1f}%` | `{sl:+.4f}` | {ghf_str} |"
        )
        if m != prop and (pm - p_score) > 0.001:
            beaters.append((m, pm, pm - p_score))

    r0_f14 = float(md["RAW_R0"]["F14_R0band_mean"])
    pr_f14 = float(md[prop]["F14_R0band_mean"])
    if pr_f14 < 0.0 or pr_f14 < r0_f14:
        degradations.append(f"盲井 `F-14` 主推 `{prop}` 为 `{pr_f14:+.4f}`（对比 `RAW_R0 = {r0_f14:+.4f}`，差值 `{pr_f14 - r0_f14:+.4f}`）")
    r0_l0 = idx_or_val(md["RAW_R0"]["F5_ZeroLag"], 0)
    pr_l0 = idx_or_val(md[prop]["F5_ZeroLag"], m1_idx)
    if pr_l0 < r0_l0:
        degradations.append(f"`F-5` 零时移 `L0` 主推 `{prop}` 为 `{pr_l0:+.4f}`（对比 `RAW_R0 = {r0_l0:+.4f}`）")

    lines.append("")
    b_str = "；".join([f"`{b[0]} = {b[1]:.4f}`（较主推 `{b[2]:+.4f}`）" for b in sorted(beaters, key=lambda x: -x[1])]) if beaters else f"本表内无任何基线在 `F5_R0band` 均值上超出主推 `{prop}` 达 `+0.001` 以上"
    d_str = "；".join(degradations) if degradations else "无盲井 `F-14` 变负或零时移 `L0` 低于 `RAW_R0` 现象"
    lines.append(f"- **反超主推 / 恶化项核查**：① 高于主推 `{prop}`（`{p_score:.4f}`）> `+0.001` 的方法：{b_str}；② 恶化/负值项核查：{d_str}。")
    lines.append("")

    # Marginal attribution table
    if all(k in md for k in ("RAW_R0", "DSA_Bent_R0", "DSA_Bent_R0_Eq", prop)):
        s0 = float(md["RAW_R0"]["F5_R0band_mean"])
        s1 = float(md["DSA_Bent_R0"]["F5_R0band_mean"])
        s2 = float(md["DSA_Bent_R0_Eq"]["F5_R0band_mean"])
        s3 = p_score
        tot = s3 - s0
        pct = lambda d: f"{(d / tot * 100.0):+.1f}%" if abs(tot) > 1e-6 else "N/A"
        ghf_p = compute_method_ghf(data, prop)
        lines.append(f"#### `{stem}` 阶梯边际归因表（`RAW_R0 → DSA_Bent_R0 → DSA_Bent_R0_Eq → {prop}`）")
        lines.append("")
        lines.append(f"| 阶梯阶段 | 方法 | `F5_R0band` 均值 | 阶段边际增量 `Δr` | 占总增量 (`{tot:+.4f}`) 比例 | `G_HF[18-24Hz]` |")
        lines.append("|---|---|---:|---:|---:|---:|")
        lines.append(f"| Step 0（无补偿原图） | `RAW_R0` | `{s0:.4f}` | `0.0000` | `0.0%` | — |")
        lines.append(f"| Step 1（+升级弯曲射线 `xmir_bent` DSA） | `DSA_Bent_R0` | `{s1:.4f}` | `{s1-s0:+.4f}` | `{pct(s1-s0)}` | — |")
        lines.append(f"| Step 2（+PSAE 前置谱锚定） | `DSA_Bent_R0_Eq` | `{s2:.4f}` | `{s2-s1:+.4f}` | `{pct(s2-s1)}` | `1.000x` |")
        lines.append(f"| Step 3（+同域 `UQW1b` Q 补偿） | **`{prop}`** | **`{s3:.4f}`** | **`{s3-s2:+.4f}`** | **`{pct(s3-s2)}`** | **`{ghf_p:.3f}x`** |")
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
        if base.endswith("_smoke.json") or any(k in base.lower() for k in ("verdict", "diagnostics", "physics_verify")):
            continue
        with open(fp, "r", encoding="utf-8") as f:
            data = json.load(f)
        if str(data.get("dataset", "")).startswith("volve"):
            chunks.append(render_field_table(fp, data, "DSA_Bent_UQW1b_PAnc_25"))
        elif "W_deep" in data:
            chunks.append(render_synth_table(fp, data, "DSA_Bent_UQW1b_25"))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as out_f:
        out_f.write("\n\n".join(chunks) + "\n")
    print(f"Rendered {len(chunks)} full tables to {args.out}")


if __name__ == "__main__":
    main()
