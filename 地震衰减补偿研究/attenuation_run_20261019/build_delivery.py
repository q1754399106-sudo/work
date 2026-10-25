"""Assemble a portable, evidence-bound report; never run an experiment.

python build_delivery.py --test-run EXPERIMENT/matlab --evaluation AUDIT \
    --figures FIGURES --output NEW_DIRECTORY [--numerical-audit POSTHOC_AUDIT]

All numerical conclusions come from a completed actual held-out MATLAB run and
the strict evaluator whose input hashes match that run. An incomplete/invalid
audit produces an explicitly partial report with effectiveness claims withheld.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone, timedelta
import hashlib
import html
import json
import math
from pathlib import Path
import shutil
from urllib.parse import quote

ROOT = Path(__file__).resolve().parent
BRANCHES = ["NoiseUncertaintyOnly", "QModelErrorOnly", "CombinedUncertainty"]
CONTROL_NAMES = ["ScalarBudgetDamping", "NoiseCovarianceOnly", "NoStructureCombined"]
KEYS = ["Model", "Suite", "SNRdB", "CalibrationSeed", "Seed", "Scenario"]
LOCAL_METRICS = ["EventAmplitudeError", "WeakEventError", "FaultNeighborhoodError",
                 "PinchoutTailError", "AVOInterceptError", "AVOGradientError"]
SUCCESS = {"COMPLETED", "COMPLETE", "SUCCESS", "SUCCEEDED", "OK", "PASS"}
LABELS = {"NoiseUncertaintyOnly": "仅噪声不确定性", "QModelErrorOnly": "仅Q模型误差",
          "CombinedUncertainty": "联合不确定性", "ScalarBudgetDamping": "普通强阻尼",
          "NoiseCovarianceOnly": "仅噪声协方差", "NoStructureCombined": "联合但无结构项"}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path, default=None):
    if not path.is_file():
        return {} if default is None else default
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def num(value):
    try:
        answer = float(value)
        return answer if math.isfinite(answer) else None
    except (ValueError, TypeError):
        return None


def fmt(value, digits=5):
    value = num(value)
    return "不可计算" if value is None else f"{value:.{digits}g}"


def pct(value, signed=False):
    value = num(value)
    return "不可计算" if value is None else (f"{value:+.2%}" if signed else f"{value:.2%}")


def esc(value):
    return html.escape(str(value), quote=True)


def link(path, label=None):
    name = str(path).replace("\\", "/")
    return f'<a href="{quote(name, safe="/#")}">{esc(label or name)}</a>'


def table(headers, rows, raw_columns=()):
    head = "".join(f"<th>{esc(x)}</th>" for x in headers)
    body = "".join("<tr>" + "".join(f"<td>{x if i in raw_columns else esc(x)}</td>"
        for i, x in enumerate(row)) + "</tr>" for row in rows)
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def copy_tree(source: Path, target: Path):
    if source.is_dir():
        shutil.copytree(source, target, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))


def observations(rows):
    return {tuple(row.get(k, "") for k in KEYS) for row in rows}


def actual_scope(rows):
    suites = []
    for suite in sorted({r.get("Suite", "") for r in rows}):
        sub = [r for r in rows if r.get("Suite") == suite]
        per_model = []
        for model in sorted({r.get("Model", "") for r in sub}):
            z = [r for r in sub if r.get("Model") == model]
            cals = sorted({r.get("CalibrationSeed", "") for r in z})
            per_model.append(dict(model=model, observations=len(observations(z)), rows=len(z),
                calibration_seeds=cals, observation_seeds_by_calibration={
                    c: sorted({r.get("Seed", "") for r in z if r.get("CalibrationSeed") == c}) for c in cals}))
        suites.append(dict(suite=suite, snr_db=sorted({r.get("SNRdB", "") for r in sub}),
            scenarios=sorted({r.get("Scenario", "") for r in sub}), models=per_model,
            observations=len(observations(sub)), rows=len(sub)))
    return dict(rows=len(rows), observations=len(observations(rows)),
        models=sorted({r.get("Model", "") for r in rows}),
        methods=sorted({r.get("Method", "") for r in rows}),
        routes=sorted({r.get("Route", "") for r in rows}), suites=suites,
        statuses=dict(Counter(r.get("Status", "MISSING") for r in rows)),
        failed_status_rows=sum(r.get("Status", "MISSING").upper() not in SUCCESS for r in rows),
        finite_waveform_errors=sum(num(r.get("RelativeError")) is not None for r in rows),
        waveform_error_gt_1_rows=sum((num(r.get("RelativeError")) or 0) > 1 for r in rows))


def resolve_spec(metadata):
    wanted = metadata.get("frozen_experiment_spec_sha256")
    if not wanted:
        return None
    candidates = list((ROOT / "experiments").rglob("frozen_experiment_spec.json"))
    candidates += [ROOT / "candidate" / "experiment" / "experiment_spec.json"]
    candidates += list((ROOT / "experiments").rglob("config_snapshot.json"))
    matches = [p for p in candidates if p.is_file() and sha(p) == wanted]
    return sorted(matches, key=lambda p: (p.name != "frozen_experiment_spec.json", str(p)))[0] if matches else None


def verify_inputs(test_run, evaluation, figures, metadata, audit, rows, specfile):
    reasons = []
    if metadata.get("stage") != "test":
        reasons.append("输入不是独立测试阶段；开发结果不能充当测试结果。")
    if metadata.get("status") != "COMPLETED":
        reasons.append("测试运行尚未标记 COMPLETED。")
    if metadata.get("actual_matlab_execution") is not True:
        reasons.append("缺少实际 MATLAB 执行证据。")
    if metadata.get("artificial_fixture") or audit.get("artificial_fixture"):
        reasons.append("人工测试夹具不能作为科学效果证据。")
    if audit.get("status") not in {"FAIL", "MIXED_RESULTS", "PROVISIONAL_PASS"}:
        reasons.append("严格审计未形成可用数值结论：" + str(audit.get("status", "MISSING")))
    if audit.get("source_and_spec_hashes_verified") is not True:
        reasons.append("严格审计未确认源码和配置哈希。")
    for field, path in [("metrics_sha256", test_run / "metrics.csv"),
                        ("metadata_sha256", test_run / "RUN_METADATA.json")]:
        if not path.is_file() or audit.get(field) != sha(path):
            reasons.append("严格审计输入绑定不一致：" + field)
    if not specfile or audit.get("frozen_experiment_spec_sha256") != sha(specfile):
        reasons.append("找不到与实际运行和审计匹配的冻结配置。")
    if not rows:
        reasons.append("没有实际 metrics.csv 观测行。")
    if rows and (metadata.get("rows") != len(rows) or audit.get("rows") != len(rows)):
        reasons.append("运行、审计与 CSV 行数不一致。")
    count = len(observations(rows))
    if rows and (metadata.get("observations") != count or audit.get("observations") != count):
        reasons.append("运行、审计与 CSV 配对观测数不一致。")
    if set(audit.get("branches", {})) != set(BRANCHES):
        reasons.append("严格审计缺少预声明的三条候选分支。")
    # Recheck recorded scientific sources before copying them into the delivery.
    for record in metadata.get("source_manifest", []):
        source = Path(record["path"])
        if not source.is_absolute():
            source = test_run / source
        if not source.is_file() or sha(source) != record["sha256"]:
            reasons.append("实际执行后源码改变或缺失：" + str(source))
    manifest = read_json(figures / "FIGURE_MANIFEST.json")
    if not reasons:
        if not manifest:
            reasons.append("缺少图件输入清单；先完成真实结果绘图。")
        elif manifest.get("row_count") != len(rows):
            reasons.append("图件行数与实际测试 CSV 不一致。")
        elif not any(r.get("sha256") == sha(test_run / "metrics.csv") for r in manifest.get("inputs", [])):
            reasons.append("图件输入未绑定到本次测试 metrics.csv。")
    for record in manifest.get("figures", []):
        for item in record.get("files", []):
            path = figures / item["path"]
            if not path.is_file() or sha(path) != item["sha256"]:
                reasons.append("图件文件哈希不一致：" + str(path))
    for item in manifest.get("source_data_files", []):
        path = figures / item["path"]
        if not path.is_file() or sha(path) != item["sha256"]:
            reasons.append("图件来源数据哈希不一致：" + str(path))
    return reasons, manifest


def local_costs(comparison):
    records = []
    for check in comparison.get("checks", []):
        metric = check["name"].rsplit("/", 1)[-1]
        if metric not in LOCAL_METRICS:
            continue
        detail = check.get("detail", {})
        b, c = num(detail.get("baseline")), num(detail.get("candidate"))
        if b is None or c is None:
            continue
        records.append(dict(gate=check["name"], metric=metric, baseline=b, candidate=c,
            delta=c-b, relative_change=(c / b - 1) if b > 0 else None,
            allowance=detail.get("allowance"), n=detail.get("n"), passed=check["passed"]))
    return sorted(records, key=lambda x: (-x["delta"], x["gate"]))


def route_summary(method, route, comp):
    local = local_costs(comp)
    worst_group = max(comp.get("groups", []), key=lambda g: g["candidate"]["mean"] /
        max(g["baseline"]["mean"], 1e-300), default=None)
    worst = local[0] if local else None
    groups = comp.get("groups", [])
    worst_value = max((g["candidate"]["worst"] for g in groups), default=None)
    return dict(method=method, route=route, status=comp.get("status", "UNAVAILABLE"),
        main_gain=comp.get("equal_model_main_improvement"),
        failed_gates=sum(not c["passed"] for c in comp.get("checks", [])),
        gates=len(comp.get("checks", [])), worst_local=worst,
        worst_waveform_error=worst_value, worst_mean_group=worst_group,
        diagnostics=comp.get("diagnostics", {}), local_costs=local)


def worst_observations(rows):
    output = []
    frozen = {tuple(r[k] for k in KEYS) + (r["Route"],): r for r in rows if r.get("Method") == "Frozen"}
    for method in BRANCHES:
        for route in ["Split", "Direct"]:
            sub = [r for r in rows if r.get("Method") == method and r.get("Route") == route
                   and num(r.get("RelativeError")) is not None]
            if not sub:
                continue
            worst = sorted(sub, key=lambda r: (-float(r["RelativeError"]), tuple(r[k] for k in KEYS)))[0]
            base = frozen.get(tuple(worst[k] for k in KEYS) + (route,), {})
            output.append(dict(method=method, route=route, key={k: worst[k] for k in KEYS},
                candidate_error=float(worst["RelativeError"]), frozen_error=num(base.get("RelativeError")),
                status=worst.get("Status"), selection="maximum finite RelativeError; all finite failed rows eligible"))
    return output


def source_snapshot(metadata, test_run, output):
    records = []
    for index, item in enumerate(metadata.get("source_manifest", []), 1):
        source = Path(item["path"])
        if not source.is_absolute():
            source = test_run / source
        source = source.resolve()
        try:
            relative = source.relative_to(ROOT)
        except ValueError:
            relative = Path("external") / f"{index:04d}_{source.name}"
        target = output / "sources" / relative
        valid = source.is_file() and sha(source) == item["sha256"]
        record = dict(original_path=str(source), delivery_path=str(Path("sources") / relative).replace("\\", "/"),
            recorded_sha256=item["sha256"], matches_execution=valid)
        if source.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            record["copied_sha256"] = sha(target)
        records.append(record)
    write_json(output / "SOURCE_MANIFEST.json", records)
    return records


def load_numerical_audit(test_run, directory=None):
    """Bind post-hoc diagnostics without touching frozen CSV or metadata bytes."""
    directory = (directory if directory is not None else
        test_run.parent.parent / (test_run.parent.name + "_numerical_audit")).resolve()
    checks = read_json(directory / "audit_checks.json")
    reasons = []
    if not checks:
        reasons.append("未找到后验逐搜索求解/校准日志审计。")
    else:
        if checks.get("status") not in {"COMPLETED_POSTHOC_AUDIT", "COMPLETED_WITH_FINDINGS"}:
            reasons.append("后验日志审计尚未完成。")
        inputs = checks.get("inputs", {})
        for name, path in [("metrics_sha256", test_run / "metrics.csv"),
                           ("run_metadata_sha256", test_run / "RUN_METADATA.json")]:
            if not path.is_file() or inputs.get(name) != sha(path):
                reasons.append("后验审计输入绑定不一致：" + name)
    return dict(status="BOUND_POSTHOC_AUDIT" if not reasons else "MISSING_OR_UNBOUND_POSTHOC_AUDIT",
        path=str(directory), report_path="evidence/numerical_audit/NUMERICAL_AUDIT_ZH.md",
        checks_path="evidence/numerical_audit/audit_checks.json", binding_verified=not reasons,
        binding_issues=reasons, checks=checks)


def numerical_audit_html(numerical):
    checks = numerical.get("checks", {})
    body = ['<section id="numerical"><h2>逐搜索求解与校准日志审计</h2>',
        '<p>这项后验审计只读取完整运行后的日志，不改变冻结CSV、配置、参数或运行元数据。'
        '原始CSV的TargetMismatch在两条路线中均记录Direct目标偏差；本节从选中求解日志单独恢复Split目标偏差，'
        '并核对搜索过程的求解标志、残差、重启和校准边界。</p>']
    if not numerical.get("binding_verified"):
        body.append('<p class="notice">' + esc('；'.join(numerical.get("binding_issues", []))) +
            ' 后验诊断不参与本报告的数值解释。</p></section>')
        return ''.join(body)
    body.append('<p>审计状态：<strong>' + esc(checks.get("status")) + '</strong>。' +
        link(numerical["report_path"], "完整审计说明") + ' · ' + link(numerical["checks_path"], "审计JSON") + '</p>')
    s = checks.get("solver_summary", {})
    labels = [("records_total", "全部搜索求解记录"), ("nonzero_final_flags", "最终非零标志"),
        ("nonzero_initial_flags", "初次非零标志"), ("normal_residual_exceedances", "正规方程残差超限"),
        ("nonfinite_or_missing_records", "非有限或缺失记录"), ("restart_or_fallback_records", "重启或回退记录"),
        ("max_normal_residual", "最大正规方程相对残差")]
    body.append(table(["求解检查", "实际记录"], [(label, s.get(key, "未记录")) for key, label in labels]))
    counts = checks.get("counts", {})
    body.append('<p>日志覆盖：' + '；'.join(f'{esc(key)}={esc(counts.get(key, "未记录"))}' for key in
        ["missing_solver_logs", "corrupt_solver_logs", "missing_calibration_logs", "corrupt_calibration_logs"]) + '。</p>')
    selected = checks.get("selected_route_summary", [])
    body.append('<h3>选中路线目标偏差</h3><p>偏差阈值为0.005；端点或目标偏差属于诊断，不能与求解失败混为一谈。'
        '单位以审计记录为准。</p>' + table(["方法", "路线", "案例数", "最终非零标志", "目标偏差>0.005", "最大目标偏差", "目标单位"],
        [(r.get("Method", r.get("method", "未记录")), r.get("Route", r.get("route", "未记录")),
          r.get("cases", "未记录"), r.get("nonzero_final_flags", "未记录"), r.get("target_mismatches_gt_0_005", "未记录"),
          fmt(r.get("max_target_mismatch")), r.get("target_units", "未记录")) for r in selected]))
    c = checks.get("calibration_summary", {})
    body.append('<h3>校准点估计及bootstrap边界</h3>' + table(
        ["估计类型", "记录数", "Q边界", "主频边界", "优化器非正退出标志"],
        [(role, c.get(prefix + "_estimates", "未记录"), c.get(prefix + "_q_boundaries", "未记录"),
          c.get(prefix + "_frequency_boundaries", "未记录"), c.get(prefix + "_optimizer_nonpositive_exitflags", "未记录"))
         for prefix, role in [("reference", "参考点估计"), ("bootstrap", "bootstrap估计")]]) +
        '<p>独立的模型/校准种子块：' + esc(json.dumps(c.get("unique_model_calibration_seed_blocks", "未记录"), ensure_ascii=False)) +
        '。同一块跨套件重复记录；参考与bootstrap行数不能当作独立地质或校准样本数。</p>')
    issues = checks.get("issues", [])
    if issues:
        body.append('<details open><summary>后验审计发现</summary><pre>' + esc(json.dumps(issues, ensure_ascii=False, indent=2)) + '</pre></details>')
    notes = checks.get("interpretation_notes", [])
    if notes:
        body.append('<details><summary>诊断解释说明</summary><pre>' + esc(json.dumps(notes, ensure_ascii=False, indent=2)) + '</pre></details>')
    body.append('<p>' + ' · '.join(link("evidence/numerical_audit/" + name, label) for name, label in [
        ("solver_records.csv", "全部搜索求解记录"), ("selected_route_diagnostics.csv", "选中路线诊断"),
        ("calibration_estimates.csv", "校准估计记录"), ("case_diagnostics.csv", "逐案例诊断"),
        ("input_manifest.csv", "原始日志输入清单")]) + '</p></section>')
    return ''.join(body)


def context_evidence(output, specfile):
    context = output / "context"
    context.mkdir()
    for source in [ROOT / "LITERATURE_DIFFERENCES_ZH.md", ROOT / "RUN_LEDGER.json", ROOT / "audit_solver_logs.py",
                   ROOT / "README_EXECUTION_ZH.md", Path(__file__).resolve(),
                   ROOT / "candidate" / "matlab" / "P1_THEORY_AND_ASSUMPTIONS.md",
                   ROOT / "candidate" / "experiment" / "METHODS_ZH_EN.md",
                   ROOT / "candidate" / "experiment" / "README_ZH.md"]:
        if source.is_file():
            shutil.copy2(source, context / source.name)
    historical = ROOT / "baseline_snapshot" / "evidence" / "fourway" / "CONCLUSIONS_ZH.md"
    if historical.is_file():
        shutil.copy2(historical, context / "HISTORICAL_CONCLUSIONS_ZH.md")
    for source in list((ROOT / "candidate" / "audit").glob("*.py")) + \
                  list((ROOT / "candidate" / "visualization").glob("*.py")) + \
                  list((ROOT / "candidate").glob("*.py")):
        target = context / "report_tools" / source.relative_to(ROOT / "candidate")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    # Prior failure attempts are retained as evidence, never treated as successful runs.
    for source_root, name in [(ROOT / "preflight", "preflight"), (ROOT / "experiments", "experiment_receipts")]:
        if source_root.is_dir():
            for source in source_root.rglob("*"):
                if not source.is_file() or "__pycache__" in source.parts:
                    continue
                if source.is_relative_to(output):
                    continue
                rel = source.relative_to(source_root)
                if name == "experiment_receipts":
                    # Retain all run receipts/CSV/logs and freeze records; examples
                    # from the actual held-out run are already copied in full.
                    allowed = source.suffix.lower() in {".csv", ".log"} or source.name in {
                        "RUN_METADATA.json", "connection.json", "run_status.json", "FREEZE_RECEIPT.json",
                        "frozen_experiment_spec.json", "config_snapshot.json", "freeze_identity.txt"}
                    if not allowed:
                        continue
                target = context / name / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
    if specfile:
        copy_tree(specfile.parent, context / "frozen_configuration")


def summary_markdown(usable, reasons, metadata, scope, audit, summaries, worst, figure_manifest, numerical):
    status = audit.get("status", "MISSING") if usable else "PARTIAL_NO_EFFECTIVENESS_CONCLUSION"
    lines = ["# 衰减补偿研究：本次执行结果", "", f"**结论：{status}；保留 Frozen，不自动替换。**", "",
        f"实际记录：{scope['observations']} 个配对观测、{scope['rows']} 条方法/路线记录；"
        f"MATLAB 执行={metadata.get('actual_matlab_execution', False)}，运行状态={metadata.get('status', 'MISSING')}。"]
    if usable:
        lines += ["", "下表为 Main 中各合成参数模型等权的误差改善；正值较好，均与同路线 Frozen 配对。", "",
            "| 候选分支 | Split 改善 | Direct 改善 | Split / Direct 结论 | 最坏局部均值代价（Split） |",
            "|---|---:|---:|---|---|"]
        for method in BRANCHES:
            s, d = summaries[(method, "Split")], summaries[(method, "Direct")]
            local = s["worst_local"]
            cost = f"{local['metric']} Δ={fmt(local['delta'])}（{pct(local['relative_change'], True)}）" if local else "不可计算"
            lines.append(f"| {method} | {pct(s['main_gain'], True)} | {pct(d['main_gain'], True)} | {s['status']} / {d['status']} | {cost} |")
        lines += ["", "每条分支均同时审查逐模型均值、局部保幅、error>1失败比例及求解残差；完整未通过项、"
            "三组机制对照和实际最坏观测见 [交互式报告](effects.html)。"]
        splitworst = max((x for x in worst if x["route"] == "Split"), key=lambda x: x["candidate_error"], default=None)
        if splitworst:
            k = splitworst["key"]
            lines += ["", f"实际最差候选 Split：{splitworst['method']}，{k['Model']}/{k['Suite']}/{k['Scenario']}，"
                f"误差 {fmt(splitworst['candidate_error'])}，同观测 Frozen {fmt(splitworst['frozen_error'])}；"
                f"校准/观测种子 {k['CalibrationSeed']}/{k['Seed']}。"]
    else:
        lines += ["", "本次不形成方法效果结论：" + "；".join(reasons), "", "未完成记录保留于附件，未计算或补造测试改善。"]
    lines += ["", "候选 Split 使用完整非稳态映射 **M·Gref** 的精确求解定义；这是相对历史分步实现的定义变化，"
        "不假定非稳态算子可精确分解。历史 Frozen 的 [**6.48%**](context/HISTORICAL_CONCLUSIONS_ZH.md) 为原归档口径，不能当作本次新结果。",
        "", "相对 Frozen 的差异还包含目标函数加权与路线定义的变化；不能把改善全部归因于 Q 误差协方差。"
        "NoiseCovarianceOnly 使用点估计 MAD/AR(1)，不是严格的白噪声点尺度单因素对照。",
        "", "Main18/LowerSNR12 各含四个白噪声场景；CoherentNoiseAVO 独立放入 StressCoherent，避免无相干扰动时的重复场景。"
        "四模型属于预选合成参数留出，并非四处实测地质。压力测试仅每模型1个校准、6个观测，重复层级弱于主试验。",
        "", "分层 bootstrap 区间仅描述此模型/校准/观测设计下的条件波动，不是总体置信区间或显著性证据。"
        "数值门槛通过也只是暂定工程结果；真实地质独立性与真值隔离仍需外部审查。",
        "", "证据入口：[完整报告与证据索引](effects.html) · [运行元数据](RUN_METADATA.json) · [源码哈希](SOURCE_MANIFEST.json)"]
    if usable:
        lines += ["", "[严格审计](evidence/evaluation/acceptance.json) · [全部测试CSV](evidence/test_run/metrics.csv)"]
    if numerical.get("binding_verified"):
        s = numerical["checks"].get("solver_summary", {})
        lines += ["", f"[后验日志审计]({numerical['report_path']})：{numerical['checks'].get('status')}；"
            f"搜索求解 {s.get('records_total', '未记录')} 条，最终非零标志 {s.get('nonzero_final_flags', '未记录')}，"
            f"残差超限 {s.get('normal_residual_exceedances', '未记录')}。Split目标偏差、重启及校准边界见完整报告。"]
    else:
        lines += ["", "后验日志诊断未绑定；全搜索目标偏差与校准边界不作完成声明。"]
    return "\n".join(lines) + "\n"


def render_html(output, usable, reasons, metadata, scope, audit, summaries, worst, manifest, passport, numerical):
    status = audit.get("status", "MISSING") if usable else "PARTIAL_NO_EFFECTIVENESS_CONCLUSION"
    body = [f'<header><p class="eyebrow">ATTENUATION RESEARCH · 实际执行报告</p><h1>误差预算与有效 Q 补偿</h1>'
        f'<p class="lead">{esc(status)} · 保留 Frozen，不自动替换。</p>'
        f'<p>{scope["observations"]} 个配对观测 / {scope["rows"]} 条方法路线记录 · '
        f'MATLAB {esc(metadata.get("matlab_version", "版本未记录"))}</p></header>',
        '<nav><a href="#branches">分支与对照</a><a href="#models">逐模型分布</a><a href="#local">局部代价</a>'
        '<a href="#gates">门槛与失败</a><a href="#numerical">逐搜索诊断</a><a href="#figures">全部图件</a><a href="#evidence">证据索引</a></nav>']
    if not usable:
        body.append('<section class="notice"><h2>未形成效果结论</h2><ul>' + ''.join(f'<li>{esc(r)}</li>' for r in reasons) + '</ul></section>')
    body += ['<section><h2>如何解释本次结果</h2><p>主目标是带限未衰减波形。改善定义为 1−候选误差/Frozen误差；'
        '主试验先按模型求均值，再对模型等权。局部代价定义为候选误差−同路线Frozen误差，正值表示退化。</p>'
        '<p>候选 Split 使用完整非稳态 <strong>M·Gref</strong> 映射的精确求解定义，不能用两个非稳态算子的乘积近似来宣称等价；'
        '这与历史分步实现存在定义变化。历史归档 Frozen 的 <a href="context/HISTORICAL_CONCLUSIONS_ZH.md"><strong>6.48%</strong></a> 属于历史数据与口径，不能与下表当作同一试验。</p>'
        '<p>归因范围：候选 Direct 在目标函数中按 Cn⁻¹ 对道加权，正则项采用 λj/σref,j²；继承的 Frozen Direct 使用原始数据范数，'
        'μ=γ·median(λ)。候选 Split 同时改为精确 M·Gref 定义。因此候选相对Frozen的差异包含这些已声明的加权和路线变化，'
        '不能单凭该差值作出Q误差协方差的纯因果收益结论。NoiseCovarianceOnly 使用点估计 MAD/AR(1)，并非严格的白噪声点尺度单因素对照。'
        '三条分支、共享先验的Direct/Split及三组机制对照均须合并判断。</p>'
        '<p>Main18 和 LowerSNR12 各保留四个白噪声场景，CoherentNoiseAVO 独立放入 StressCoherent；执行前去除主套件中的重复白噪声场景。'
        '四个测试单元是预选的合成层参数留出，不能解释为四处独立实测地质。压力套件每模型仅1个校准×6个观测，不能估计跨校准稳定性。</p>'
        '<p>分层 bootstrap 区间是此设计下的条件性描述区间，不是总体置信区间、显著性检验或真实地质泛化证明。'
        '数值审计只给暂定工程判断；真值隔离、物理校准独立性和冻结时间线仍需外部核查。</p></section>']
    suite_rows = [(s["suite"], ", ".join(s["snr_db"]), ", ".join(s["scenarios"]),
                   s["observations"], s["rows"]) for s in scope["suites"]]
    body.append('<section><h2>实际执行范围</h2>' + table(["套件", "SNR dB", "场景", "配对观测", "方法路线行"], suite_rows) +
        f'<p>运行状态：{esc(metadata.get("status", "MISSING"))}；全部失败保留：{esc(metadata.get("all_failures_retained", False))}；'
        f'状态失败行：{scope["failed_status_rows"]}。原始 CSV 与种子层级见元数据。</p></section>')
    if usable:
        overview = []
        for method in BRANCHES + CONTROL_NAMES:
            for route in ["Split", "Direct"]:
                comp = audit["branches"][method]["routes"][route] if method in BRANCHES else audit.get("controls", {}).get(method, {}).get(route, {})
                s = summaries.get((method, route), route_summary(method, route, comp))
                overview.append((LABELS.get(method, method) + " / " + method, route, s["status"], pct(s["main_gain"], True),
                    f'{s["failed_gates"]}/{s["gates"]}', fmt(s["worst_waveform_error"]),
                    json.dumps(s["diagnostics"], ensure_ascii=False)))
        body.append('<section id="branches"><h2>三分支及三组机制对照</h2><p>Split 是预声明主路线。Direct 共享该方法的先验信息；'
            '对照采用相同评分规则，但不会改变候选资格。通过状态不授权替换 Frozen。表中选择诊断来自原始CSV；'
            'CSV的TargetMismatch两路线均为Direct目标偏差，Split日志诊断单列于后验审计。</p>' +
            table(["方法", "路线", "结论", "Main改善", "未过/总门槛", "最坏单行误差", "选择诊断"], overview) + '</section>')
        groups = []
        local_rows = []
        gates = []
        bootstrap_rows = []
        for method in BRANCHES + CONTROL_NAMES:
            for route in ["Split", "Direct"]:
                comp = audit["branches"][method]["routes"][route] if method in BRANCHES else audit.get("controls", {}).get(method, {}).get(route, {})
                for g in comp.get("groups", []):
                    for side, label in [("baseline", "Frozen"), ("candidate", method)]:
                        d = g[side]
                        groups.append((g["Suite"], g["Model"], g["SNRdB"], method, route, label,
                            d["n"], fmt(d["mean"]), fmt(d["median"]), fmt(d["p90"]), fmt(d["p95"]),
                            fmt(d["worst"]), pct(d["failure_fraction"])))
                costs = local_costs(comp)
                for metric in LOCAL_METRICS:
                    sub = [c for c in costs if c["metric"] == metric]
                    if sub:
                        c = sub[0]
                        local_rows.append((method, route, metric, c["gate"].rsplit("/", 1)[0], fmt(c["baseline"]),
                            fmt(c["candidate"]), fmt(c["delta"]), pct(c["relative_change"], True), fmt(c["allowance"]),
                            c["n"], "通过" if c["passed"] else "未过"))
                failed_checks = [c for c in comp.get("checks", []) if not c["passed"]]
                for c in failed_checks:
                    gates.append((method, route, c["name"], json.dumps(c["detail"], ensure_ascii=False)))
                u = comp.get("hierarchical_uncertainty", {})
                interval = u.get("relative_improvement_interval")
                bootstrap_rows.append((method, route, u.get("replicates", "N/A"), u.get("geological_models", "N/A"),
                    " ~ ".join(pct(x, True) for x in interval) if interval else "不可计算"))
        body.append('<section id="models"><h2>全部逐模型分布</h2><p>每个套件/模型/SNR均提供均值、中位数、p90、p95、最大值和误差&gt;1比例。'
            '相同Frozen分布在各比较中重复显示，以便核对每个比较。无失败行被删除。</p><details><summary>展开完整模型表</summary>' +
            table(["套件", "模型", "SNR", "比较分支", "路线", "被评分方法", "n", "均值", "中位数", "p90", "p95", "最大值", "误差>1比例"], groups) + '</details></section>')
        body.append('<section id="local"><h2>每项局部指标的最坏均值代价</h2><p>对每方法/路线/指标，选择候选减Frozen的绝对差最大的套件/模型/SNR/场景组。'
            '相对百分比仅辅助解释，接近零的参照使用冻结绝对容差。所有组的门槛值保存在严格审计JSON。</p>' +
            table(["方法", "路线", "指标", "最坏组", "Frozen", "候选", "绝对差", "相对变化", "允许绝对差", "n", "门槛"], local_rows) + '</section>')
        body.append('<section><h2>实际最坏观测</h2><p>每个候选/路线选最大有限RelativeError，固定ID/种子打破并列；'
            '有限求解失败同样有资格入选。非有限失败仅在CSV和失败记录中报告。</p>' + table(
            ["方法", "路线", "模型/套件/场景", "校准/观测种子", "候选误差", "Frozen误差", "状态"],
            [(w["method"], w["route"], "/".join(w["key"][k] for k in ["Model", "Suite", "Scenario"]),
              w["key"]["CalibrationSeed"] + "/" + w["key"]["Seed"], fmt(w["candidate_error"]), fmt(w["frozen_error"]), w["status"]) for w in worst]) + '</section>')
        body.append('<section id="gates"><h2>全部未通过门槛</h2>' + (table(["方法", "路线", "门槛", "原始详情"], gates) if gates else '<p>无未通过数值门槛；工程状态仍是暂定。</p>') +
            f'<p>严格审计记录的求解失败行数：{len(audit.get("failed_rows", []))}；Baseline数值检查：{esc(audit.get("baseline_control_numerically_valid"))}。</p>' +
            '<details><summary>冻结工程阈值与附加检查</summary><p>Main模型等权改善至少5%；各套件/模型均值退化不超过2%；'
            '事件、弱反射、AVO局部退化不超过5%，尖灭尾部10%，近零绝对容差0.001；误差&gt;1比例不得增加；'
            '最终求解标志为0且正规方程相对残差≤2×10⁻⁷。精确冻结值以封存配置为准。</p></details></section>')
        body.append('<section><h2>条件性描述 bootstrap</h2><p>先重采样合成模型，再校准层，再观测种子块；'
            '配对关系保留。下列2.5–97.5百分位范围仅描述这个设计，不构成总体置信区间。</p>' +
            table(["方法", "路线", "重采样次数", "模型数", "改善描述范围"], bootstrap_rows) + '</section>')
    body.append(numerical_audit_html(numerical))
    figure_links = []
    for target, label in [("evidence/figures/index.html", "独立图册"), ("evidence/figures/FIGURE_MANIFEST.json", "图件清单")]:
        if (output / target).is_file():
            figure_links.append(link(target, label))
    body.append('<section id="figures"><h2>全部真实结果图</h2><p>图件状态：' + esc(manifest.get("status", "MISSING")) + '。'
        '图件清单包含绘图变换、缺失值政策、输入哈希和来源数据。图件保留原始振幅与实际最坏案例。</p>' +
        ' · '.join(figure_links) + '</section>')
    for figure in manifest.get("figures", []):
        png = next((f["path"] for f in figure.get("files", []) if f["path"].endswith(".png")), None)
        pdf = next((f["path"] for f in figure.get("files", []) if f["path"].endswith(".pdf")), None)
        if png:
            body.append(f'<figure><h3>{esc(figure["title"])}</h3><img loading="lazy" '
                f'src="{quote("evidence/figures/" + png, safe="/")}" alt="{esc(figure.get("alt_text", ""))}">'
                f'<figcaption>{esc(figure.get("alt_text", ""))} ' + (link("evidence/figures/" + pdf, "PDF") if pdf else "") + '</figcaption></figure>')
    body.append('<section id="evidence"><h2>完整证据索引</h2><p>' +
        ' · '.join([link("EFFECTS_ZH.md", "一页结果"), link("RUN_METADATA.json", "交付运行元数据"),
        link("SOURCE_MANIFEST.json", "执行源码清单"), link("FILE_MANIFEST.json", "全部文件哈希"),
        link("context/LITERATURE_DIFFERENCES_ZH.md", "已核对文献与差异")]) + '</p>' +
        '<p>下面列出交付包中全部 CSV、JSON、日志、文本与源代码，包含预检失败尝试和人工夹具；'
        '人工夹具仅验证评价程序，开发集与历史回放均不计入本次测试效果。</p><details><summary>展开全部证据文件</summary><ul>')
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.suffix.lower() in {".csv", ".json", ".log", ".txt", ".m", ".py", ".md"}:
            relative = path.relative_to(output).as_posix()
            body.append('<li>' + link(relative) + '</li>')
    body.append('</ul></details><details><summary>Material Passport / 材料护照</summary><pre>' +
        esc(json.dumps(passport, ensure_ascii=False, indent=2)) + '</pre><p>UNVERIFIED指未完成外部学术完整性核验；'
        '本地字节哈希与数值矩阵核查的范围分别记录在元数据中，不能升级为科学真实性认证。</p></details></section>')
    css = '''*{box-sizing:border-box}body{margin:0;background:#f4f6fa;color:#182334;font:15px/1.65 "Microsoft YaHei",system-ui,sans-serif}main{max-width:1440px;margin:auto;padding:28px}header{background:#142d49;color:#fff;padding:36px;border-radius:16px}h1{font-size:34px;line-height:1.25;margin:8px 0 16px}h2{font-size:23px;color:#173c63}h3{font-size:18px}.eyebrow{font-size:12px;letter-spacing:1.6px;color:#a8d4f4}.lead{font-size:21px;font-weight:700}nav{display:flex;flex-wrap:wrap;gap:18px;padding:22px 4px}a{color:#0866a7;text-decoration-thickness:1px;text-underline-offset:3px}section,figure{margin:20px 0;padding:26px;background:#fff;border:1px solid #dce3ec;border-radius:12px}figure{margin-left:0;margin-right:0}img{display:block;width:100%;height:auto}figcaption{color:#4c5d70;font-size:13px;margin-top:12px}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:12px;font-variant-numeric:tabular-nums}th,td{border-bottom:1px solid #dce3ec;padding:9px 11px;text-align:left;vertical-align:top}th{background:#edf3f9;position:sticky;top:0}tr:nth-child(even){background:#f8fafc}td{overflow-wrap:anywhere;min-width:72px}summary{cursor:pointer;color:#0866a7;font-weight:700;padding:12px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}.notice{border-left:6px solid #b96800;background:#fff7e9}li{overflow-wrap:anywhere;margin:3px 0}@media(max-width:700px){main{padding:12px}section,figure,header{padding:18px}h1{font-size:27px}}@media print{nav,details{display:none}section,figure{break-inside:avoid;border:0}body{background:white}main{padding:0}}'''
    return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>衰减补偿 · 实际研究结果</title><style>' + css + '</style></head><body><main>' + ''.join(body) + '</main></body></html>'


def build(test_run: Path, evaluation: Path, figures: Path, output: Path, numerical_audit: Path | None = None):
    test_run, evaluation, figures, output = (p.resolve() for p in (test_run, evaluation, figures, output))
    if output.exists():
        raise ValueError("Output must be a fresh directory: " + str(output))
    # Avoid recursively copying a delivery into any of its own input trees.
    numerical = load_numerical_audit(test_run, numerical_audit)
    if any(output == p or output.is_relative_to(p) for p in [test_run, evaluation, figures, Path(numerical["path"])]):
        raise ValueError("Output must be outside all three input directories")
    metadata = read_json(test_run / "RUN_METADATA.json")
    audit = read_json(evaluation / "acceptance.json")
    metricfile = test_run / "metrics.csv"
    with metricfile.open(encoding="utf-8-sig", newline="") if metricfile.is_file() else _empty_stream() as stream:
        rows = list(csv.DictReader(stream))
    specfile = resolve_spec(metadata)
    reasons, manifest = verify_inputs(test_run, evaluation, figures, metadata, audit, rows, specfile)
    usable = not reasons
    scope = actual_scope(rows)
    summaries = {}
    if usable:
        for method in BRANCHES + CONTROL_NAMES:
            for route in ["Split", "Direct"]:
                comp = audit["branches"][method]["routes"][route] if method in BRANCHES else audit.get("controls", {}).get(method, {}).get(route, {})
                summaries[(method, route)] = route_summary(method, route, comp)
    worst = worst_observations(rows) if usable else []
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    copy_tree(test_run, evidence / "test_run")
    copy_tree(evaluation, evidence / "evaluation")
    copy_tree(figures, evidence / "figures")
    copy_tree(Path(numerical["path"]), evidence / "numerical_audit")
    # MATLAB connector logs live beside its matlab/ output, not inside it.
    for name in ["connection.json", "matlab_stdout.log", "matlab_stderr.log"]:
        source = test_run.parent / name
        if source.is_file():
            shutil.copy2(source, evidence / name)
    sources = source_snapshot(metadata, test_run, output)
    context_evidence(output, specfile)
    timestamp = datetime.now(timezone(timedelta(hours=8))).isoformat()
    passport = dict(origin_skill="academic-research-suite / experiment-agent", origin_mode="execution-report",
        origin_date=timestamp, verification_status="UNVERIFIED", version_label="attenuation_delivery_v1",
        repro_lock=None, upstream_dependencies=["actual_test_RUN_METADATA", "strict_acceptance", "FIGURE_MANIFEST"],
        note="Local numerical/hash audits only; no independent external scientific audit or human-read attestation.")
    result = dict(schema_version="attenuation-delivery-1", assembled_at=timestamp,
        status="COMPLETED_REPORT" if usable else "PARTIAL_NO_EFFECTIVENESS_CONCLUSION",
        numerical_acceptance=audit.get("status", "MISSING"), effectiveness_summary_permitted=usable,
        reasons_effectiveness_withheld=reasons, actual_matlab_execution=metadata.get("actual_matlab_execution", False),
        actual_test_run_status=metadata.get("status", "MISSING"), actual_test_stage=metadata.get("stage", "MISSING"),
        actual_scope=scope, source_run_metadata="evidence/test_run/RUN_METADATA.json",
        actual_run_metadata=metadata, strict_audit="evidence/evaluation/acceptance.json",
        audit_input_hashes={k: audit.get(k) for k in ["metrics_sha256", "metadata_sha256", "frozen_experiment_spec_sha256"]},
        input_directories=dict(test_run=str(test_run), evaluation=str(evaluation), figures=str(figures)),
        source_manifest="SOURCE_MANIFEST.json", source_files=len(sources),
        source_snapshot_matches_execution=bool(sources) and all(s["matches_execution"] for s in sources),
        figures_status=manifest.get("status", "MISSING"), figure_count=manifest.get("figure_count", 0),
        posthoc_numerical_audit=numerical,
        no_automatic_frozen_replacement=True, historical_frozen_gain=dict(value=.0648,
            scope="Historical archived data/definition; never a new held-out result"),
        split_definition="Exact full nonstationary M*Gref mapping; differs from historical split implementation",
        attribution_limit="Candidate Direct uses Cn^-1 trace weighting and lambda_j/sigma_ref_j^2 regularization; inherited Frozen Direct uses raw-data norm and mu=gamma*median(lambda). Candidate Split also changes to exact M*Gref. Differences versus Frozen cannot isolate a causal Q-error covariance benefit. NoiseCovarianceOnly uses point MAD/AR(1), not a strict white point-scale-only control.",
        independence="Preselected synthetic-parameter holdout; field-geology independence unestablished",
        resampling="Conditional descriptive hierarchical bootstrap; not a population confidence interval",
        main_scenario_deviation="Four white-noise scenes; coherent scene in separate StressCoherent to remove duplicates",
        stress_replication="One calibration and six observations per model, weaker than Main/LowerSNR",
        candidate_route_summaries=[summaries[k] for k in summaries], worst_observations=worst,
        material_passport=passport, packaging_exclusions=["dependency runtimes", "cache files"],
        manifest_policy="FILE_MANIFEST hashes all payload files except itself and SHA256SUMS.txt to avoid recursive hashes")
    write_json(output / "RUN_METADATA.json", result)
    (output / "EFFECTS_ZH.md").write_text(summary_markdown(usable, reasons, metadata, scope, audit, summaries, worst, manifest, numerical), encoding="utf-8")
    # Create manifest placeholders so all evidence links are included in HTML.
    write_json(output / "FILE_MANIFEST.json", [])
    (output / "SHA256SUMS.txt").write_text("", encoding="utf-8")
    (output / "effects.html").write_text(render_html(output, usable, reasons, metadata, scope, audit, summaries, worst, manifest, passport, numerical), encoding="utf-8")
    files = [dict(path=p.relative_to(output).as_posix(), bytes=p.stat().st_size, sha256=sha(p))
             for p in sorted(output.rglob("*")) if p.is_file() and p.name not in {"FILE_MANIFEST.json", "SHA256SUMS.txt"}]
    write_json(output / "FILE_MANIFEST.json", files)
    (output / "SHA256SUMS.txt").write_text(''.join(f'{f["sha256"]}  {f["path"]}\n' for f in files), encoding="utf-8")
    return dict(status=result["status"], numerical_acceptance=result["numerical_acceptance"],
        rows=scope["rows"], observations=scope["observations"], files=len(files), output=str(output),
        reasons_effectiveness_withheld=reasons)


def _empty_stream():
    import io
    return io.StringIO("")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["test-run", "evaluation", "figures", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--numerical-audit", type=Path, default=None,
        help="Post-hoc audit directory; defaults to test_run.parent.parent/<test-tag>_numerical_audit")
    args = parser.parse_args()
    try:
        result = build(args.test_run, args.evaluation, args.figures, args.output, args.numerical_audit)
    except (ValueError, OSError, KeyError, TypeError) as error:
        raise SystemExit("DELIVERY_REFUSED: " + str(error))
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
