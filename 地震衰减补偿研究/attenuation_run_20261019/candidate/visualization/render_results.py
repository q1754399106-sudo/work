"""Truthful, deterministic P1 scientific figures from saved metrics/examples.

Usage: python render_results.py --metrics RUN/metrics.csv --examples RUN/examples
                               --output NEW_DIRECTORY [--worst-cases RUN/worst_cases]
No solver or MATLAB graphics are called. Results are never selected by appearance.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import shutil
import sys

RUN_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RUN_ROOT / "dependencies"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch

METHODS = ["Baseline", "Frozen", "NoiseUncertaintyOnly", "QModelErrorOnly",
           "CombinedUncertainty", "ScalarBudgetDamping", "NoiseCovarianceOnly",
           "NoStructureCombined"]
LABELS = ["常规 Baseline", "主版 Frozen", "噪声不确定 Noise U", "Q误差 Q error",
          "联合 Combined", "普通强阻尼 Damping", "噪声协方差 Noise C", "无结构 No H"]
COLORS = ["#555555", "#B35400", "#0072B2", "#9B4E93", "#007B59", "#827300", "#466C96", "#A62653"]
MARKERS = ["s", "o", "^", "D", "P", "v", "X", "h"]
STYLE = {m: {"label": LABELS[i], "color": COLORS[i], "marker": MARKERS[i]} for i, m in enumerate(METHODS)}
PAIR_KEYS = ["Model", "Suite", "Scenario", "Seed", "CalibrationSeed", "Route"]
LOCAL = [("EventAmplitudeError", "事件振幅 / Event amplitude"),
         ("WeakEventError", "弱反射 / Weak event"),
         ("PinchoutTailError", "尖灭尾部 / Pinchout tail"),
         ("AVOInterceptError", "AVO截距 / AVO intercept"),
         ("AVOGradientError", "AVO梯度 / AVO gradient"),
         ("FaultNeighborhoodError", "断层邻域 / Fault neighborhood")]
EXAMPLE_METHODS = ["Frozen_Split", "CombinedUncertainty_Direct", "CombinedUncertainty_Split"]
EXAMPLE_LABELS = {"truth": "真值 / Truth", "observation": "观测 / Observed",
                  "Frozen_Split": "主版分步 / Frozen Split", "CombinedUncertainty_Direct": "联合一步 / Combined Direct",
                  "CombinedUncertainty_Split": "联合分步 / Combined Split"}
EXAMPLE_STYLE = {"truth": ("#000000", "-"), "observation": ("#999999", ":"),
                 "Frozen_Split": ("#B35400", "-"), "CombinedUncertainty_Direct": ("#0072B2", "--"),
                 "CombinedUncertainty_Split": ("#007B59", "-.")}
CANDIDATE_BRANCHES = ["NoiseUncertaintyOnly", "QModelErrorOnly", "CombinedUncertainty"]


def series_label(name):
    if name in EXAMPLE_LABELS:
        return EXAMPLE_LABELS[name]
    method, _, route = name.rpartition("_")
    return f'{STYLE.get(method, {"label": method})["label"]} / {route}'


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def plain(value):
    if value is pd.NA:
        return None
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def setup_font() -> Path:
    for p in [Path("C:/Windows/Fonts/msyh.ttc"), Path("C:/Windows/Fonts/simhei.ttf")]:
        if p.exists():
            font_manager.fontManager.addfont(str(p))
            name = font_manager.FontProperties(fname=str(p)).get_name()
            plt.rcParams.update({"font.family": name, "font.size": 9, "axes.titlesize": 10,
                                 "figure.titlesize": 12, "legend.fontsize": 8,
                                 "axes.unicode_minus": False, "pdf.fonttype": 42,
                                 "ps.fonttype": 42, "savefig.facecolor": "white"})
            return p
    raise RuntimeError("A Chinese font is required: Microsoft YaHei or SimHei was not found")


class Renderer:
    def __init__(self, metrics: Path, examples: Path, output: Path, worst_cases: Path | None = None):
        if output.exists():
            raise FileExistsError(f"Output must be fresh: {output}")
        output.mkdir(parents=True)
        self.output = output
        self.data = output / "source_data"
        self.data.mkdir()
        self.font = setup_font()
        self.metrics_path = metrics.resolve()
        self.df = pd.read_csv(metrics, keep_default_na=True)
        self.df.insert(0, "SourceRow", np.arange(1, len(self.df) + 1))
        required = set(PAIR_KEYS + ["Method", "RelativeError"])
        if required - set(self.df):
            raise ValueError(f"Missing metrics columns: {sorted(required - set(self.df))}")
        self.pair_keys = PAIR_KEYS + (["SNRdB"] if "SNRdB" in self.df else [])
        self.methods = [m for m in METHODS if m in set(self.df.Method)]
        self.methods += sorted(set(self.df.Method) - set(self.methods))
        for m in self.methods:
            if m not in STYLE:
                STYLE[m] = {"label": m, "color": "#333333", "marker": "o"}
        self.models = list(dict.fromkeys(self.df.Model.astype(str)))
        for c in ["RelativeError", "NormalResidual", "FinalFlag"] + [p[0] for p in LOCAL]:
            if c in self.df:
                self.df[c] = pd.to_numeric(self.df[c], errors="coerce")
        self.df.to_csv(self.data / "metrics_all_rows.csv", index=False)
        self.inputs = [{"path": str(metrics.resolve()), "sha256": digest(metrics)}]
        self.examples = []
        for path in sorted(examples.glob("*.json")):
            ex = json.loads(path.read_text(encoding="utf-8-sig"))
            if "truth" in ex and "observation" in ex and "methods" in ex:
                ex["_path"] = path
                self.examples.append(ex)
        self.worst_cases = []
        self.worst_cases_path = worst_cases if worst_cases is not None else examples.parent / "worst_cases"
        for path in sorted(self.worst_cases_path.glob("*.json")):
            ex = json.loads(path.read_text(encoding="utf-8-sig"))
            if "truth" in ex and "observation" in ex and "methods" in ex:
                ex["_path"] = path
                self.worst_cases.append(ex)
        self.worst_selection = None
        self.figures = []
        self.exported_examples = set()
        self.summary = self.make_summary()
        self.summary.to_csv(self.data / "distribution_summary_all_rows.csv", index=False)
        self.paired = self.make_pairs()
        self.paired.to_csv(self.data / "paired_all_rows.csv", index=False)

    def save(self, fig, stem, title, description, transformations, sources):
        # Constrained-layout rect is (left,bottom,width,height), not top/right.
        # Reserve a real top margin; several panels also reserve a footer.
        engine = fig.get_layout_engine()
        if engine is not None and hasattr(engine, 'get'):
            left, bottom, width, height = engine.get().get('rect', (0, 0, 1, 1))
            engine.set(rect=(left, bottom, width, min(height, .92-bottom)))
        fig.suptitle(title, y=.985)
        files = []
        for suffix in ["png", "pdf"]:
            target = self.output / f"{stem}.{suffix}"
            fig.savefig(target, dpi=180, facecolor="white", transparent=False,
                        metadata={"Title": title, "Creator": "P1 render_results.py"})
            files.append({"path": target.name, "sha256": digest(target)})
        self.figures.append({"id": stem, "title": title, "alt_text": description,
                             "transformations": transformations, "source_data": sources,
                             "files": files, "size_inches": list(fig.get_size_inches())})
        plt.close(fig)

    def make_summary(self):
        rows = []
        for keys, g in self.df.groupby(["Suite", "Model", "Method", "Route"], dropna=False, sort=False):
            x = g.RelativeError.to_numpy(float)
            f = x[np.isfinite(x)]
            row = dict(zip(["Suite", "Model", "Method", "Route"], keys))
            row.update({"rows": len(g), "finite": len(f), "nonfinite": len(g) - len(f),
                        "mean": np.mean(f) if len(f) else np.nan, "median": np.median(f) if len(f) else np.nan,
                        "p90": np.quantile(f, .9) if len(f) else np.nan,
                        "p95": np.quantile(f, .95) if len(f) else np.nan,
                        "worst": np.max(f) if len(f) else np.nan,
                        "error_gt_1": np.sum(f > 1),
                        "failed_status": int(np.sum(g.Status.fillna("MISSING") != "SUCCESS")) if "Status" in g else 0})
            rows.append(row)
        return pd.DataFrame(rows)

    def make_pairs(self):
        fields = ["RelativeError"] + [c for c, _ in LOCAL if c in self.df]
        ref = self.df[self.df.Method == "Frozen"][self.pair_keys + fields].copy()
        if ref.duplicated(self.pair_keys).any():
            raise ValueError("Frozen pairing keys are duplicated; figures refuse ambiguous pairing")
        ref = ref.rename(columns={c: f"Frozen_{c}" for c in fields})
        merged = self.df.merge(ref, on=self.pair_keys, how="left", indicator="FrozenPair", validate="many_to_one")
        for c in fields:
            merged[f"Delta_{c}"] = merged[c] - merged[f"Frozen_{c}"]
        den = merged.Frozen_RelativeError
        merged["RelativeChange"] = np.where(np.isfinite(den) & (den != 0), merged.Delta_RelativeError / den, np.nan)
        return merged

    def legend(self, fig, methods=None, bottom=.03, columns=4):
        from matplotlib.lines import Line2D
        methods = self.methods if methods is None else methods
        handles = [Line2D([], [], color=STYLE[m]["color"], marker=STYLE[m]["marker"], linestyle="none",
                          label=STYLE[m]["label"]) for m in methods]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, bottom), ncol=columns, frameon=False)

    def model_errors(self):
        fig, axs = plt.subplots(2, 2, figsize=(15, 9), layout="constrained")
        fig.get_layout_engine().set(rect=(0, .12, 1, .94))
        offsets = np.linspace(-.32, .32, len(self.methods))
        for row, suite in enumerate(["Main", "LowerSNR"]):
            for col, route in enumerate(["Direct", "Split"]):
                ax = axs[row, col]
                sub = self.df[(self.df.Suite == suite) & (self.df.Route == route)]
                missing = int((~np.isfinite(sub.RelativeError)).sum())
                for i, m in enumerate(self.methods):
                    for j, model in enumerate(self.models):
                        raw = sub[(sub.Method == m) & (sub.Model == model)].RelativeError.to_numpy(float)
                        x = raw[np.isfinite(raw)]
                        if len(x):
                            mean = np.mean(x); low, high = np.quantile(x, [.1, .9])
                            ax.vlines(j + offsets[i], low, high, color=STYLE[m]["color"], lw=1)
                            ax.scatter(j + offsets[i], mean, marker=STYLE[m]["marker"], color=STYLE[m]["color"], s=28)
                ax.set_title(f"{suite} | {route} | 全部行 / all rows={len(sub)}; 非有限 / missing={missing}")
                ax.set_xticks(range(len(self.models)), self.models, rotation=17, ha="right")
                ax.set_ylabel("相对波形误差 / Relative waveform error")
                ax.set_ylim(bottom=0); ax.grid(axis="y", alpha=.2)
        self.legend(fig)
        self.save(fig, "01_model_errors", "逐模型波形误差 / Waveform error by model",
                  "Four panels compare every method for each model in Main and LowerSNR, separately for Direct and Split. Points are finite-row means, vertical lines span the empirical 10th–90th percentiles. Numerical failures with finite scores remain included; nonfinite counts are explicit.",
                  ["Main/LowerSNR are displayed separately", "Mean and empirical p10–p90 within each model; descriptive, not a confidence interval", "No success-status filter; nonfinite values counted, never converted to zero"],
                  ["source_data/metrics_all_rows.csv", "source_data/distribution_summary_all_rows.csv"])

    def paired_changes(self):
        methods = [m for m in self.methods if m != "Frozen"]
        fig, axs = plt.subplots(2, 2, figsize=(15, 9), layout="constrained")
        fig.get_layout_engine().set(rect=(0, .13, 1, .94))
        offsets = np.linspace(-.3, .3, len(methods))
        grouped = []
        for row, suite in enumerate(["Main", "LowerSNR"]):
            for col, route in enumerate(["Direct", "Split"]):
                ax = axs[row, col]; count = 0
                sub = self.paired[(self.paired.Suite == suite) & (self.paired.Route == route)]
                for i, m in enumerate(methods):
                    for j, model in enumerate(self.models):
                        z = sub[(sub.Method == m) & (sub.Model == model)].RelativeChange.to_numpy(float)
                        finite = z[np.isfinite(z)]; count += len(z) - len(finite)
                        grouped.append({"Suite": suite, "Route": route, "Method": m, "Model": model,
                                        "rows": len(z), "missing_or_zero_denominator": len(z)-len(finite),
                                        "mean_relative_change": np.mean(finite) if len(finite) else np.nan})
                        if len(finite):
                            x = j + offsets[i]
                            ax.scatter(np.full(len(finite), x), 100*finite, s=6, color=STYLE[m]["color"], alpha=.18, rasterized=True)
                            ax.scatter(x, 100*np.mean(finite), s=33, color=STYLE[m]["color"], marker=STYLE[m]["marker"])
                ax.axhline(0, color="black", lw=.8, label="与Frozen相同 / equal to Frozen")
                ax.set_title(f"{suite} | {route} | 缺失配对或分母零 / unresolved={count}")
                ax.set_xticks(range(len(self.models)), self.models, rotation=17, ha="right")
                ax.set_ylabel("配对误差变化 / Paired error change (%)\n负值较好 / negative is better")
                ax.grid(axis="y", alpha=.2)
        pd.DataFrame(grouped).to_csv(self.data / "paired_model_changes.csv", index=False)
        self.legend(fig, methods)
        self.save(fig, "02_paired_changes", "相对Frozen的全部配对差 / All paired changes versus Frozen",
                  "Individual paired percentage changes are shown as faint points and each model mean as a larger symbol. Zero means equality with the same-route Frozen. Every paired finite score is included, with unresolved counts disclosed.",
                  ["Pairing uses " + "/".join(self.pair_keys), "100*(method error-Frozen error)/Frozen error; a zero reference denominator is unresolved", "All finite paired rows retained, including failed-status rows"],
                  ["source_data/paired_all_rows.csv", "source_data/paired_model_changes.csv"])

    def local_metrics(self, route, stem):
        methods = [m for m in self.methods if m != "Frozen"]
        fig, axs = plt.subplots(2, 3, figsize=(16, 9), layout="constrained")
        fig.get_layout_engine().set(rect=(0, .1, 1, .94))
        table = []
        for ax, (field, label) in zip(axs.flat, LOCAL):
            missing = 0
            for i, m in enumerate(methods):
                values = []
                for model in self.models:
                    z = self.paired[(self.paired.Route == route) & self.paired.Suite.isin(["Main", "LowerSNR"]) &
                                    (self.paired.Method == m) & (self.paired.Model == model)][f"Delta_{field}"].to_numpy(float)
                    finite = z[np.isfinite(z)]; missing += len(z)-len(finite)
                    value = np.mean(finite) if len(finite) else np.nan
                    table.append({"Route": route, "Metric": field, "Method": m, "Model": model,
                                  "finite_applicable_pairs": len(finite), "nonfinite_or_not_applicable": len(z)-len(finite), "mean_delta": value})
                    if len(finite):
                        values.append(value)
                        ax.scatter(i, value, color=STYLE[m]["color"], marker="o", s=18, alpha=.5)
                if values:
                    ax.scatter(i, np.mean(values), color=STYLE[m]["color"], marker=STYLE[m]["marker"], s=65)
            ax.axhline(0, color="black", lw=.8)
            ax.set_title(f"{label}\n不可数或不适用 / missing or N/A={missing}")
            ax.set_xticks(range(len(methods)), [STYLE[m]["label"] for m in methods], rotation=48, ha="right", fontsize=7)
            ax.set_ylabel("误差差值 / Error delta\n负值较好 / negative is better")
            ax.grid(axis="y", alpha=.2)
        name = f"local_metrics_{route}.csv"
        pd.DataFrame(table).to_csv(self.data / name, index=False)
        from matplotlib.lines import Line2D
        fig.legend(handles=[Line2D([], [], color="#555555", marker="o", ls="none", markersize=4, label="各模型均值 / Model mean"),
                            Line2D([], [], color="#555555", marker="D", ls="none", markersize=7, label="模型等权均值 / Equal-model mean")],
                   loc="lower center", ncol=2, frameon=False)
        self.save(fig, stem, f"{route}局部保幅代价 / Local amplitude costs ({route})",
                  "Six panels display paired absolute changes in event amplitude, weak-event, pinchout-tail, AVO intercept, AVO gradient and fault-neighborhood errors. Small circles show model means, larger method symbols show the equal-model average. Missing and inapplicable entries are counted together and retained in source data.",
                  ["Main and LowerSNR only; all finite applicable pairs", "Absolute method-minus-Frozen metric difference", "Small circles=model means; larger symbols=equal-model mean; no inferential intervals", "NaN can mean metric inapplicable; no NaN is turned into zero"],
                  ["source_data/paired_all_rows.csv", f"source_data/{name}"])

    def failures(self):
        fig, axs = plt.subplots(1, 2, figsize=(15, 6), layout="constrained")
        fig.get_layout_engine().set(rect=(0, .12, 1, .93))
        table = []
        colors = ["#A62653", "#827300", "#007B59"]
        labels = ["求解失败或非有限 / numerical or invalid", "有限误差>1 / finite error>1", "有限误差≤1 / finite error≤1"]
        for ax, route in zip(axs, ["Direct", "Split"]):
            for i, m in enumerate(self.methods):
                g = self.df[(self.df.Route == route) & (self.df.Method == m)]
                x = g.RelativeError.to_numpy(float)
                bad = (~np.isfinite(x)) | (x < 0)
                if "Status" in g: bad |= g.Status.fillna("MISSING").to_numpy() != "SUCCESS"
                if "FinalFlag" in g: bad |= (~np.isfinite(g.FinalFlag.to_numpy(float))) | (g.FinalFlag.to_numpy(float) != 0)
                if "NormalResidual" in g: bad |= (~np.isfinite(g.NormalResidual.to_numpy(float))) | (g.NormalResidual.to_numpy(float) > 2e-7)
                counts = np.array([np.sum(bad), np.sum(~bad & (x > 1)), np.sum(~bad & (x <= 1))])
                bottom = 0
                for count, color in zip(counts, colors):
                    ax.bar(i, count/max(len(g), 1), bottom=bottom, color=color, width=.7, edgecolor="white", lw=.4)
                    bottom += count/max(len(g), 1)
                ax.text(i, 1.02, str(len(g)), ha="center", fontsize=8)
                table.append({"Route": route, "Method": m, "all_rows": len(g), "numerical_or_invalid": counts[0],
                              "valid_error_gt_1": counts[1], "valid_error_le_1": counts[2], "nonfinite_error": np.sum(~np.isfinite(x))})
            ax.set_title(f"{route} | 全部suite / all suites")
            ax.set_xticks(range(len(self.methods)), [STYLE[m]["label"] for m in self.methods], rotation=40, ha="right")
            ax.set_ylim(0, 1.12); ax.set_ylabel("全部行比例 / Fraction of all rows")
            ax.grid(axis="y", alpha=.2)
        pd.DataFrame(table).to_csv(self.data / "failure_counts_all_rows.csv", index=False)
        fig.legend(handles=[Patch(color=c, label=l) for c, l in zip(colors, labels)], loc="lower center", ncol=3, frameon=False)
        self.save(fig, "05_failures", "不隐藏失败 / Failures retained in the full denominator",
                  "Two stacked bar charts classify every row from every suite as numerical/invalid, a valid waveform error above one, or a valid error at most one. Numbers above bars give all-row denominators; nonfinite errors remain in the invalid category.",
                  ["All suites included, with no removed rows", "Numerical/invalid includes non-success Status, nonzero FinalFlag, nonfinite scores/residuals or normal residual>2e-7", "Stacked fractions have a zero baseline and full all-row denominator"],
                  ["source_data/metrics_all_rows.csv", "source_data/failure_counts_all_rows.csv"])

    def distributions(self):
        fig, axs = plt.subplots(1, 2, figsize=(15, 6), layout="constrained")
        fig.get_layout_engine().set(rect=(0, .14, 1, .93))
        for ax, route in zip(axs, ["Direct", "Split"]):
            nonfinite = 0
            for m in self.methods:
                raw = self.df[(self.df.Route == route) & (self.df.Method == m)].RelativeError.to_numpy(float)
                x = np.sort(raw[np.isfinite(raw)]); nonfinite += len(raw)-len(x)
                if len(x):
                    ax.step(x, np.arange(1, len(x)+1)/len(raw), where="post", color=STYLE[m]["color"],
                            marker=STYLE[m]["marker"], markevery=max(1, len(x)//10), ms=4, lw=1.2)
            ax.axvline(1, color="black", ls="--", lw=.8)
            ax.set_title(f"{route} | 非有限 / nonfinite={nonfinite}")
            ax.set_xlabel("相对误差（完整范围） / Relative error (full range)")
            ax.set_ylabel("有限误差累计数 / 全部行数\nFinite cumulative count / all rows")
            all_route = self.df[self.df.Route == route].RelativeError.to_numpy(float)
            finite_route = all_route[np.isfinite(all_route)]
            ax.set_ylim(0, 1.02); ax.set_xlim(left=min(0, np.min(finite_route)) if len(finite_route) else 0); ax.grid(alpha=.2)
        self.legend(fig)
        self.save(fig, "06_distributions", "全部分布及完整最差尾部 / Full distributions including the worst tail",
                  "Empirical cumulative curves include every finite RelativeError, regardless of solver status, across all suites. Their denominator includes nonfinite rows, so an incomplete curve can end below one. The full observed error range is shown and the error=1 threshold is marked.",
                  ["All suites and all finite scores; no status filter or tail clipping", "ECDF denominator is all rows, including nonfinite", "Pooled descriptive distribution; repeated noise rows are not independent geology", "Mean/median/p90/p95/worst exported per model/suite/method/route"],
                  ["source_data/metrics_all_rows.csv", "source_data/distribution_summary_all_rows.csv"])

    def stress_and_worst(self):
        fig, axs = plt.subplots(1, 2, figsize=(16, 7), layout="constrained")
        fig.get_layout_engine().set(rect=(0, .04, 1, .94))
        sub = self.paired[(self.paired.Route == "Split") & (self.paired.Method == "CombinedUncertainty") & ~self.paired.Suite.isin(["Main", "LowerSNR"])]
        suites = list(dict.fromkeys(sub.Suite.astype(str)))
        matrix = np.full((len(self.models), max(len(suites), 1)), np.nan); records=[]
        for i, model in enumerate(self.models):
            for j, suite in enumerate(suites):
                z = sub[(sub.Model == model) & (sub.Suite == suite)].RelativeChange.to_numpy(float); finite=z[np.isfinite(z)]
                value=np.mean(finite) if len(finite) else np.nan; matrix[i,j]=100*value
                records.append({"Model":model,"Suite":suite,"mean_paired_relative_change":value,"rows":len(z),"nonfinite":len(z)-len(finite)})
        finite=matrix[np.isfinite(matrix)]; lim=max(float(np.max(np.abs(finite))) if len(finite) else 1, 1e-9)
        cmap=matplotlib.colormaps["BrBG_r"].with_extremes(bad="#BBBBBB")
        im=axs[0].imshow(matrix,cmap=cmap,norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),aspect="auto",interpolation="nearest")
        axs[0].set_xticks(range(len(suites)),suites,rotation=55,ha="right",fontsize=8)
        axs[0].set_yticks(range(len(self.models)),self.models,fontsize=8)
        axs[0].set_title("联合分步压力配对变化 / Combined Split stress changes\n灰色=无有效配对 / gray=no finite pair")
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                axs[0].text(j,i,f"{matrix[i,j]:.1f}" if np.isfinite(matrix[i,j]) else "N/A",ha="center",va="center",fontsize=7)
        fig.colorbar(im,ax=axs[0],label="配对变化 / Paired change (%)")
        worst=[]
        for m in self.methods:
            g=self.df[self.df.Method==m]; finite=g[np.isfinite(g.RelativeError)]
            if len(finite):
                row=finite.loc[finite.RelativeError.idxmax()]
                worst.append([STYLE[m]["label"],f"{row.RelativeError:.4g}",str(row.Model),f"{row.Suite}/{row.Scenario}",str(row.Route)])
        if not worst:
            worst = [["无有限结果 / No finite result", "N/A", "", "", ""]]
        axs[1].axis("off")
        table=axs[1].table(cellText=worst,colLabels=["方法 / Method","最差 / Worst","模型 / Model","Suite / Scene","路线 / Route"],loc="center",cellLoc="left",colWidths=[.25,.11,.23,.3,.11])
        table.auto_set_font_size(False); table.set_fontsize(7); table.scale(1,2)
        axs[1].set_title("每方法全部行最差案例 / Worst finite case for every method\n非有限失败另见失败图 / nonfinite failures remain in failure chart")
        pd.DataFrame(records).to_csv(self.data/"stress_paired_changes.csv",index=False)
        self.df[np.isfinite(self.df.RelativeError)].sort_values("RelativeError",ascending=False).head(30).to_csv(self.data/"worst_30_finite_rows.csv",index=False)
        self.save(fig,"07_stress_worst","压力情景与最差案例 / Stress suites and worst cases",
                  "The left heatmap contains mean paired percentage changes for Combined Split versus Frozen Split in every stress suite and model, with symmetric full-range colors and gray missing cells. The right table names the worst finite error for every method across all rows and routes. Nonfinite failures are reported separately rather than omitted from the study.",
                  ["Stress heatmap uses every paired finite stress row with no success-status filter", "Heatmap limits are ±maximum absolute displayed change; missing=gray", "Worst cases selected explicitly by maximum finite RelativeError, not appearance"],
                  ["source_data/stress_paired_changes.csv","source_data/worst_30_finite_rows.csv"])

    def choose_example(self,suite,scenario):
        candidates=[e for e in self.examples if e.get("Suite")==suite and e.get("Scenario")==scenario]
        # Fixed selection order uses IDs and seeds only, never errors/amplitudes.
        return sorted(candidates,key=lambda e:(str(e.get("Model")),int(e.get("CalibrationSeed",0)),int(e.get("Seed",0)),str(e["_path"])))[:1]

    def arrays(self,ex):
        arrays={"truth":np.asarray(ex["truth"],float),"observation":np.asarray(ex["observation"],float)}
        for key in EXAMPLE_METHODS:
            if key in ex["methods"]: arrays[key]=np.asarray(ex["methods"][key],float)
        return arrays

    def export_example(self,ex,arrays=None,prefix=""):
        path=ex["_path"]; stem=prefix+path.stem
        if stem in self.exported_examples: return f"source_data/{stem}.csv"
        self.exported_examples.add(stem); self.inputs.append({"path":str(path.resolve()),"sha256":digest(path)})
        t=np.asarray(ex["time"],float).ravel(); arrays=self.arrays(ex) if arrays is None else arrays; records=[]
        for series,S in arrays.items():
            for j in range(S.shape[1]):
                records.append(pd.DataFrame({"Series":series,"Trace":j+1,"Time_s":t,"Amplitude":S[:,j]}))
        pd.concat(records,ignore_index=True).to_csv(self.data/f"{stem}.csv",index=False)
        return f"source_data/{stem}.csv"

    def example_caption(self,ex):
        return f"{ex['Model']} | {ex['Suite']} | {ex['Scenario']} | cal={ex['CalibrationSeed']} | seed={ex['Seed']}"

    def missing_example(self,stem,title,suite,scenario):
        fig,ax=plt.subplots(figsize=(10,4),layout="constrained"); ax.axis("off")
        ax.text(.5,.5,f"未保存指定案例 / Designated example unavailable\n{suite} / {scenario}\n不以其他结果替代 / No result-based replacement",ha="center",va="center",fontsize=13)
        self.save(fig,stem,title,"The predesignated example was not saved. This placeholder records the missing artifact without selecting an alternative result.",["No replacement by result quality"],[])

    def gather(self,suite,scenario,stem):
        selected=self.choose_example(suite,scenario); title=f"指定道集 / Predesignated gather: {suite} / {scenario}"
        if not selected: return self.missing_example(stem,title,suite,scenario)
        ex=selected[0]; arrays=self.arrays(ex); t=np.asarray(ex["time"],float).ravel(); source=self.export_example(ex)
        finite=np.concatenate([s[np.isfinite(s)] for s in arrays.values()]); lim=max(np.max(np.abs(finite)) if len(finite) else 1,1e-15)
        keys=["truth", "observation"]+EXAMPLE_METHODS
        fig,axs=plt.subplots(1,len(keys),figsize=(17,5.6),layout="constrained",sharey=True)
        fig.get_layout_engine().set(rect=(0,.07,1,.87))
        for ax,name in zip(np.atleast_1d(axs),keys):
            if name not in arrays:
                ax.axis("off"); ax.set_title(EXAMPLE_LABELS[name]); ax.text(.5,.5,"方法失败或未保存\nMethod failed or not saved",ha="center",va="center",transform=ax.transAxes)
                continue
            S=arrays[name]
            im=ax.imshow(np.ma.masked_invalid(S),origin="upper",extent=(.5,S.shape[1]+.5,t[-1],t[0]),aspect="auto",interpolation="nearest",cmap=matplotlib.colormaps["RdBu_r"].with_extremes(bad="#888888"),vmin=-lim,vmax=lim)
            ax.set_title(EXAMPLE_LABELS[name]); ax.set_xlabel("道 / Trace"); ax.set_xticks([1,6,11]); ax.set_ylabel("时间 / Time (s)")
        fig.colorbar(im,ax=np.atleast_1d(axs).tolist(),label="原始振幅（共享完整范围） / Amplitude (shared full range)",shrink=.8)
        fig.text(.5,.025,self.example_caption(ex),ha="center",fontsize=9)
        self.save(fig,stem,title,"Five same-gather panels show truth, observation, Frozen Split and Combined Direct/Split for the predesignated first model/calibration/observation. Every panel uses the same symmetric amplitude range, set by the largest absolute value across all shown arrays, so no sample is clipped.",
                  ["Predesignated example selected by suite/scenario/model ID/seeds only", "No per-panel normalization, filtering, smoothing or amplitude clipping",f"Shared amplitude limits ±{lim:.12g}; gray=nonfinite"],[source])

    def waveforms(self):
        selected=self.choose_example("Main","SmoothAVO"); stem="11_fixed_traces_waveforms_spectra"; title="固定1/6/11道 / Fixed traces 1/6/11: waveform and spectrum"
        if not selected: return self.missing_example(stem,title,"Main","SmoothAVO")
        ex=selected[0]; arrays=self.arrays(ex); t=np.asarray(ex["time"],float).ravel(); dt=float(np.mean(np.diff(t))); source=self.export_example(ex)
        fig,axs=plt.subplots(3,2,figsize=(15,10),layout="constrained"); fig.get_layout_engine().set(rect=(0,.09,1,.93))
        spectra=[]
        for row,trace in enumerate([1,6,11]):
            for k,(name,S) in enumerate(arrays.items()):
                y=S[:,trace-1]; f=np.fft.rfftfreq(len(y),dt); a=np.abs(np.fft.rfft(y))/len(y)
                color,style=EXAMPLE_STYLE[name]
                axs[row,0].plot(t,y,color=color,ls=style,lw=1,label=EXAMPLE_LABELS[name])
                axs[row,1].plot(f,a,color=color,ls=style,lw=1,label=EXAMPLE_LABELS[name])
                spectra.append(pd.DataFrame({"Trace":trace,"Series":name,"Frequency_Hz":f,"Magnitude_FFT_div_N":a}))
            axs[row,0].set_title(f"第{trace}道波形 / Trace {trace} waveform"); axs[row,0].set_xlabel("时间 / Time (s)"); axs[row,0].set_ylabel("原始振幅 / Amplitude")
            axs[row,1].set_title(f"第{trace}道频谱 / Trace {trace} spectrum"); axs[row,1].set_xlabel("频率 / Frequency (Hz)"); axs[row,1].set_ylabel("|FFT|/N"); axs[row,1].set_xlim(0,.5/dt); axs[row,1].set_ylim(bottom=0)
            for ax in axs[row]: ax.grid(alpha=.2)
        handles,labels=axs[0,0].get_legend_handles_labels(); fig.legend(handles,labels,loc="lower center",ncol=3,frameon=False)
        missing=[name for name in EXAMPLE_METHODS if name not in arrays]
        if missing:
            fig.text(.5,.065,"失败或未保存 / Failed or not saved: "+", ".join(missing),ha="center",fontsize=8)
        pd.concat(spectra,ignore_index=True).to_csv(self.data/"fixed_trace_spectra.csv",index=False)
        self.save(fig,stem,title,"Three fixed traces, 1, 6 and 11, are shown with their full time waveforms and unsmoothed Fourier magnitudes over the full Nyquist range. Truth, observed data, Frozen Split and both Combined routes retain raw amplitude scaling and distinct line styles.",
                  ["Fixed traces [1,6,11], selected before results", "Spectrum=abs(rfft(signal))/N; no taper, smoothing, interpolation, normalization to individual maxima or frequency truncation", "Full time extent and 0–Nyquist frequency extent"],[source,"source_data/fixed_trace_spectra.csv"])

    def local_zoom(self):
        fig,axs=plt.subplots(2,3,figsize=(15,8),layout="constrained"); fig.get_layout_engine().set(rect=(0,.08,1,.92))
        sources=[]
        for row,scenario in enumerate(["WeakSection","PinchoutSection"]):
            selected=self.choose_example("Main",scenario)
            if not selected:
                for ax in axs[row]: ax.axis("off"); ax.text(.5,.5,f"指定案例缺失 / Missing\nMain/{scenario}",ha="center",va="center")
                continue
            ex=selected[0]; sources.append(self.export_example(ex)); arrays=self.arrays(ex); t=np.asarray(ex["time"],float).ravel(); use=(t>=.45)&(t<=.85)
            traces=np.arange(1,12) if row==0 else np.arange(9,12); ids=traces-1; keys=["truth","Frozen_Split","CombinedUncertainty_Split"]
            pieces=[arrays[k][use][:,ids] for k in keys if k in arrays]; finite=np.concatenate([z[np.isfinite(z)] for z in pieces]); lim=max(np.max(np.abs(finite)) if len(finite) else 1,1e-15)
            for col,key in enumerate(keys):
                ax=axs[row,col]
                if key not in arrays: ax.axis("off"); ax.text(.5,.5,"方法失败或缺失 / Method failed or missing",ha="center",va="center"); continue
                im=ax.imshow(np.ma.masked_invalid(arrays[key][use][:,ids]),origin="upper",extent=(traces[0]-.5,traces[-1]+.5,t[use][-1],t[use][0]),aspect="auto",interpolation="nearest",cmap=matplotlib.colormaps["RdBu_r"].with_extremes(bad="#888888"),vmin=-lim,vmax=lim)
                ax.set_title(f"{scenario}\n{EXAMPLE_LABELS[key]}"); ax.set_xlabel("道 / Trace"); ax.set_ylabel("时间 / Time (s)"); ax.set_xticks(traces if row else [1,6,11])
            fig.colorbar(im,ax=axs[row].tolist(),label=f"本行共享振幅 / Row amplitude ±{lim:.3g}",shrink=.8)
        fig.text(.5,.02,"固定窗0.45–0.85 s；尖灭固定9–11道 / Fixed window 0.45–0.85 s; pinchout traces 9–11",ha="center",fontsize=9)
        self.save(fig,"12_weak_pinchout_zoom","弱反射与尖灭局部 / Predesignated weak-event and pinchout zooms",
                  "Two rows compare truth, Frozen Split and Combined Split for the predesignated weak-event and pinchout examples. The time window is fixed at 0.45–0.85 seconds; pinchout uses fixed tail traces 9–11. Each row shares its full absolute amplitude maximum across shown methods, with no clipping.",
                  ["Predesignated Main/WeakSection and Main/PinchoutSection examples", "Fixed time window [0.45,.85] s; no result-based peak picking", "Weak event all traces; pinchout tail traces [9,10,11]", "Each row uses a shared unclipped symmetric range across all shown arrays"],sources)

    def choose_worst_candidate(self):
        eligible=[]
        for ex in self.worst_cases:
            if ex.get("Method") not in CANDIDATE_BRANCHES or ex.get("Route") != "Split":
                continue
            value=ex.get("RelativeError")
            if value is None or not np.isfinite(float(value)):
                continue
            target=f"{ex['Method']}_Split"
            if target in ex["methods"]:
                eligible.append(ex)
        if not eligible:
            return None
        # This intentionally chooses the largest error, with fixed ID/seed tie breaks.
        ex=sorted(eligible,key=lambda e:(-float(e["RelativeError"]),CANDIDATE_BRANCHES.index(e["Method"]),
                                       str(e.get("Model")),str(e.get("Suite")),str(e.get("Scenario")),
                                       int(e.get("CalibrationSeed",0)),int(e.get("Seed",0))))[0]
        frozen=self.df[(self.df.Method=="Frozen") & (self.df.Route=="Split")].copy()
        candidate=self.df[(self.df.Method==ex["Method"]) & (self.df.Route=="Split")].copy()
        for key in ["Model","Suite","Scenario","Seed","CalibrationSeed"]:
            frozen=frozen[frozen[key]==ex[key]]
            candidate=candidate[candidate[key]==ex[key]]
        self.worst_selection={k:ex.get(k) for k in ["Method","Route","Model","Suite","Scenario","Seed","CalibrationSeed","Status","RelativeError","FailureMessage"]}
        self.worst_selection.update({"input_path":str(ex["_path"].resolve()),
                                    "rule":"Largest finite RelativeError among NoiseUncertaintyOnly/QModelErrorOnly/CombinedUncertainty Split; includes finite solver-failed results; ties use declared branch order then IDs/seeds",
                                    "frozen_comparison_rows":frozen.to_dict(orient="records"),
                                    "candidate_metric_rows":candidate.to_dict(orient="records"),
                                    "nonfinite_candidate_rows_not_visually_reconstructable":int((~np.isfinite(self.df[self.df.Method.isin(CANDIDATE_BRANCHES) & (self.df.Route=="Split")].RelativeError)).sum())})
        return ex

    def worst_caption(self,ex):
        frozen_rows=self.worst_selection["frozen_comparison_rows"] if self.worst_selection else []
        fvalue=frozen_rows[0].get("RelativeError") if len(frozen_rows)==1 else None
        frozen_text=f"Frozen error={fvalue:.5g}" if fvalue is not None and np.isfinite(fvalue) else "Frozen error unavailable/ambiguous"
        return (f"{ex['Method']}_Split | {self.example_caption(ex)}\n"
                f"Status={ex.get('Status','MISSING')} | candidate error={float(ex['RelativeError']):.5g} | {frozen_text}")

    def worst_failure_figures(self):
        ex=self.choose_worst_candidate()
        title="实际最差候选案例 / Actual worst candidate case"
        if ex is None:
            self.missing_example("13_worst_candidate_gather",title+": gather","worst_cases","three candidate Split branches")
            self.missing_example("14_worst_candidate_waveforms_spectra",title+": waveform/spectrum","worst_cases","three candidate Split branches")
            return
        target=f"{ex['Method']}_Split"
        arrays={"truth":np.asarray(ex["truth"],float),"observation":np.asarray(ex["observation"],float)}
        for name in ["Frozen_Split",target]:
            if name in ex["methods"]: arrays[name]=np.asarray(ex["methods"][name],float)
        source=self.export_example(ex,arrays=arrays,prefix="worst_case_")
        t=np.asarray(ex["time"],float).ravel(); dt=float(np.mean(np.diff(t)))
        finite=np.concatenate([S[np.isfinite(S)] for S in arrays.values()]); lim=max(float(np.max(np.abs(finite))) if len(finite) else 1,1e-15)
        fig,axs=plt.subplots(1,4,figsize=(16,6.2),layout="constrained",sharey=True)
        fig.get_layout_engine().set(rect=(0,.12,1,.9))
        for ax,name in zip(axs,["truth","observation","Frozen_Split",target]):
            ax.set_title(series_label(name))
            if name not in arrays:
                ax.axis("off"); ax.text(.5,.5,"方法失败或未保存\nMethod failed or not saved",ha="center",va="center",transform=ax.transAxes)
                continue
            S=arrays[name]
            im=ax.imshow(np.ma.masked_invalid(S),origin="upper",extent=(.5,S.shape[1]+.5,t[-1],t[0]),aspect="auto",interpolation="nearest",
                         cmap=matplotlib.colormaps["RdBu_r"].with_extremes(bad="#888888"),vmin=-lim,vmax=lim)
            ax.set_xlabel("道 / Trace"); ax.set_xticks([1,6,11]); ax.set_ylabel("时间 / Time (s)")
        fig.colorbar(im,ax=axs.tolist(),label="共享完整原始振幅 / Shared full raw amplitude",shrink=.8)
        fig.text(.5,.025,self.worst_caption(ex),ha="center",fontsize=9)
        self.save(fig,"13_worst_candidate_gather",title+": 道集 / gather",
                  "The actual largest finite RelativeError among the three declared candidate Split branches is shown beside its same-observation Frozen Split, truth and observation. Selection uses maximum error, including finite solver-failed results, with deterministic tie breaks; it does not choose an attractive example. The key, status and errors are printed. All panels share an unclipped absolute amplitude maximum.",
                  [self.worst_selection["rule"],"Same observation, calibration seed and observation seed for candidate and Frozen",f"Shared limits ±{lim:.12g}; no amplitude clipping or per-panel normalization", "Nonfinite failures remain in CSV/failure counts and cannot supply a finite recovered gather"],[source])
        fig,axs=plt.subplots(3,2,figsize=(15,10.8),layout="constrained")
        fig.get_layout_engine().set(rect=(0,.16,1,.93))
        spectra=[]
        styles={"truth":("#000000","-"),"observation":("#999999",":"),"Frozen_Split":("#B35400","--"),target:(STYLE[ex["Method"]]["color"],"-.")}
        for row,trace in enumerate([1,6,11]):
            for name,S in arrays.items():
                color,style=styles[name]; y=S[:,trace-1]; f=np.fft.rfftfreq(len(y),dt); a=np.abs(np.fft.rfft(y))/len(y)
                axs[row,0].plot(t,y,color=color,ls=style,lw=1,label=series_label(name))
                axs[row,1].plot(f,a,color=color,ls=style,lw=1,label=series_label(name))
                spectra.append(pd.DataFrame({"Trace":trace,"Series":name,"Frequency_Hz":f,"Magnitude_FFT_div_N":a}))
            axs[row,0].set_title(f"最差案例固定第{trace}道 / Worst-case fixed trace {trace}")
            axs[row,0].set_xlabel("时间 / Time (s)"); axs[row,0].set_ylabel("原始振幅 / Amplitude")
            axs[row,1].set_title(f"第{trace}道完整频谱 / Trace {trace} full spectrum")
            axs[row,1].set_xlabel("频率 / Frequency (Hz)"); axs[row,1].set_ylabel("|FFT|/N"); axs[row,1].set_xlim(0,.5/dt); axs[row,1].set_ylim(bottom=0)
            for ax in axs[row]: ax.grid(alpha=.2)
        handles,labels=axs[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.085),ncol=2,frameon=False)
        fig.text(.5,.018,self.worst_caption(ex),ha="center",fontsize=8)
        pd.concat(spectra,ignore_index=True).to_csv(self.data/"worst_candidate_fixed_trace_spectra.csv",index=False)
        self.save(fig,"14_worst_candidate_waveforms_spectra",title+": 波形与频谱 / waveform and spectrum",
                  "For the same actual worst finite candidate Split case, fixed traces 1, 6 and 11 show the full-time raw waveforms and unsmoothed full-Nyquist spectra of truth, observation, the candidate and its matched Frozen Split. The case identifiers, solver status and candidate/Frozen errors are printed. Trace choice is fixed rather than selected by the appearance of failure.",
                  [self.worst_selection["rule"],"Fixed traces [1,6,11]; full time and frequency extent", "Spectrum=abs(rfft)/N; no smoothing, taper, frequency clipping or individual normalization"],
                  [source,"source_data/worst_candidate_fixed_trace_spectra.csv"])

    def finish(self):
        manifest={"status":"RENDERED_REQUIRES_VISUAL_INSPECTION","audience":"local engineering research report; publisher requirements not specified",
                  "script":{"path":str(Path(__file__).resolve()),"sha256":digest(Path(__file__))},
                  "inputs":self.inputs,"font":{"path":str(self.font),"sha256":digest(self.font)},
                  "dependencies":{"python":sys.version,"numpy":np.__version__,"pandas":pd.__version__,"matplotlib":matplotlib.__version__},
                  "row_count":len(self.df),"all_rows_preserved_in_source_csv":True,"figure_count":len(self.figures),
                  "pairing_keys":self.pair_keys,"worst_case_selection":self.worst_selection,
                  "figure_selection":"First saved predesignated IDs/seeds for named suite/scenario; no score/appearance selection",
                  "missing_policy":"Numeric summaries use finite values, including solver-failed finite scores; nonfinite counts explicit. Original NaN/N/A retained. No value replaced by zero.",
                  "replicate_policy":"Descriptive observations nested within geology/calibration. No rowwise confidence intervals or significance claims.",
                  "source_data_files":[{"path":str(p.relative_to(self.output)),"sha256":digest(p)} for p in sorted(self.data.glob("*"))],
                  "figures":self.figures}
        (self.output/"FIGURE_MANIFEST.json").write_text(json.dumps(plain(manifest),ensure_ascii=False,indent=2),encoding="utf-8")
        (self.output/"ALT_TEXT.json").write_text(json.dumps([{k:f[k] for k in ["id","title","alt_text"]} for f in self.figures],ensure_ascii=False,indent=2),encoding="utf-8")
        entries=[]
        for f in self.figures:
            entries.append(f'<section><h2>{html.escape(f["title"])}</h2><a href="{f["id"]}.pdf">PDF</a><img src="{f["id"]}.png" alt="{html.escape(f["alt_text"],quote=True)}"><p>{html.escape(f["alt_text"])}</p></section>')
        page='<!doctype html><html lang="zh"><meta charset="utf-8"><title>P1 figures</title><style>body{font-family:Microsoft YaHei,Arial,sans-serif;max-width:1500px;margin:auto;padding:24px}img{width:100%;height:auto}section{margin:32px 0}p{line-height:1.5}</style><h1>P1结果图 / P1 result figures</h1><p>全部原始行保留于 <a href="source_data/metrics_all_rows.csv">metrics_all_rows.csv</a>。图中区间为描述分位数，不代表统计置信度。All rows remain available; quantiles are descriptive. <a href="FIGURE_MANIFEST.json">Manifest</a></p>'+''.join(entries)+'</html>'
        (self.output/"index.html").write_text(page,encoding="utf-8")

    def run(self):
        self.model_errors(); self.paired_changes()
        self.local_metrics("Direct","03_local_Direct"); self.local_metrics("Split","04_local_Split")
        self.failures(); self.distributions(); self.stress_and_worst()
        self.gather("Main","SmoothAVO","08_gather_Main")
        self.gather("LowerSNR","SmoothAVO","09_gather_LowerSNR")
        self.gather("StressSourceShape","SmoothAVO","10_gather_SourceShape")
        self.waveforms(); self.local_zoom(); self.worst_failure_figures(); self.finish()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics",type=Path,required=True)
    parser.add_argument("--examples",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--worst-cases",type=Path,default=None,help="Defaults to examples.parent/worst_cases")
    args=parser.parse_args()
    renderer=Renderer(args.metrics,args.examples,args.output,args.worst_cases)
    renderer.run()
    print(json.dumps({"status":"RENDERED_REQUIRES_VISUAL_INSPECTION","output":str(args.output.resolve()),"figures":len(renderer.figures)},ensure_ascii=False))


if __name__=="__main__":
    main()
