# 方案：锚定理论归一化相图与子波提取深度标度律（T3，`anchor_norm_phasediag_20261123`，开跑前提交）

## 0. 物理动机与统一闭式推导

此前 `anchor_active_hessian_20261120` 与 `anchor_closedform_confirm_20261121` 在单一稀疏正则化规则（中位列范数加权 `R_med`）与固定浅层子波窗（$t_w = 1.4\text{ s}$）下验证了固定活跃集闭式敏感度 $\alpha_F = -J_{v,v_a}/J_{vv}$。本实验将该理论推广到 **L1 列归一化规则 × 子波提取窗深度** 二维相图，并检验柴新涛等（Chai et al. 2019）标准单位列范数 L1 (`R_unit`) 在 QScan 中的系统偏差来源。

令未归一化传播字典 $\tilde{B}(v, v_a) = G(a v) W_{\text{in}}(v_a)$，列范数 $s_k(v, v_a) = \|\tilde{B}_{:,k}\|_2$，单位列归一化字典 $B = \tilde{B}\,\mathrm{diag}(1/s)$。在未归一化反射系数坐标 $y = \mathrm{diag}(1/s) X$ 下，三种 L1 归一化规则统一对应单道目标函数：
$$F(v, v_a) = \frac{1}{2}\|d - \tilde{B} y\|_2^2 + \sigma \sum_{k=1}^n p_k(v, v_a)\, |y_k|$$
其中惩罚权重向量 $p_k(v, v_a)$ 在基准点 $(v_0, v_a=1)$ 的中位列 $k_{\text{med}}$ 处完全相等（均为 $\lambda\sqrt{2\ln n}\,s_{\text{ref}}$，$\lambda=0.5$），仅在随深度 $k$ 与 $(v, v_a)$ 的演化方式上正交分离：
1. **`R_med`（现行管线）**：$p_k(v, v_a) = \lambda\sqrt{2\ln n}\,\mathrm{median}_j s_j(v, v_a)$（对 $k$ 均匀，随 $(v,v_a)$ 通过中位列范数 $\mu(v,v_a)$ 漂移）；
2. **`R_fix`（固定参考范数）**：$p_k(v, v_a) \equiv \lambda\sqrt{2\ln n}\,s_{\text{ref}}$（对 $k$ 与 $(v,v_a)$ 均恒定，精确关断全部 $\partial\mu/\partial v$ 与 $\partial\mu/\partial v_a$ 耦合项）；
3. **`R_unit`（Chai et al. 2019 标准单位列范数 L1）**：$p_k(v, v_a) = \lambda\sqrt{2\ln n}\,s_k(v, v_a)$（等价于在单位列字典 $B$ 上施加均匀 L1 惩罚 $\lambda\sqrt{2\ln n}\,\sigma\|X\|_1$，深层原子因衰减导致 $s_k$ 随 $v$ 指数收缩）。

在固定活跃集 $A$ 与符号 $\mathbf{s}_A = \mathrm{sign}(y_A)$ 上，三种规则共享**完全精确的统一闭式能量公式**：
$$\mathbf{c}_A = \tilde{B}_A^\top d - \sigma (\mathbf{p}_A \odot \mathbf{s}_A), \qquad F(v, v_a) = \frac{1}{2}\|d\|_2^2 - \frac{1}{2}\mathbf{c}_A^\top (\tilde{B}_A^\top \tilde{B}_A)^{-1}\mathbf{c}_A$$
由此单次求解即可按二阶有限差分直接给出三规则通用的闭式斜率预测 $\alpha_F = -J_{v,v_a}/J_{vv}$。

## 1. 预注册实验矩阵（全新独立随机种子 `9931, 9932`）

- **3 个正交 $Q$ 剖面**（均在 $c_{\text{true}} = 1.0$ 生成）：
  - `S`（浅层强衰减）：`[30, 200, 200, 200]`
  - `D`（深层强衰减）：`[250, 250, 40, 40]`
  - `N3`（中深夹层衰减）：`[100, 300, 40, 300]`
- **3 种 L1 归一化规则**：`R_med`, `R_fix`, `R_unit`
- **3 个子波提取窗中心深度**：$t_w \in \{1.4, 2.0, 2.5\}\text{ s}$（窗 $[t_w - 0.2, t_w + 0.2]\text{ s}$，以 $a(t_w)$ 为锚定累积吸收量）
- **总单元数（36 组全新合成，不含任何已见旧种子）**：
  - 核心三维相图（`seed = 9931`）：`3 profiles × 3 rules × 3 tw = 27` 组；
  - 独立种子复现层（`seed = 9932`, `tw = 1.4 s`）：`3 profiles × 3 rules = 9` 组。
- 每组执行 3 点锚定全网格扫描（$c_a \in \{1.0, 0.7, 1.4\}$，得实测真值 $\hat{\kappa} = (v_{\text{sel}}(0.7) - v_{\text{sel}}(1.4))/(1/0.7 - 1/1.4)$ 与零偏差基准 $v_{\text{sel}}(1.0)$）+ 1 次在 $v_0 = v_{\text{sel}}(1.0)$ 的统一闭式 $\alpha_F$ 预测与能量闭合自检。

## 2. 预注册机械判定准则（`scripts/judge_phasediag.py`）

- **`N0_Closure`（数学自检）**：所有收敛单元中，固定活跃集闭式能量与 FISTA 目标值相对差 $|J_{\text{fix}} - J_{\text{solve}}|/|J_{\text{solve}}| \le 10^{-3}$。
- **`N1_ClosedForm`（三规则通用闭式预测精度）**：在所有非边界收敛单元上，$\mathrm{median}|\alpha_F - \hat{\kappa}| \le 0.05$；且在 `R_med` 与 `R_fix` 子集上 $\mathrm{median}|\alpha_F - \hat{\kappa}| \le 0.03$。
- **`N2_RuleShift`（归一化规则相位分离与 Chai 2019 偏差解释）**：
  - 量化 `R_fix` 相对 `R_med` 的 $\hat{\kappa}$ 偏移（验证关断 $\mu$-耦合后 $\alpha$ 的系统变化方向与幅值）；
  - 量化 `R_unit` 在真实锚点 $c_a=1.0$ 处相对真值 $v_{\text{true}}=1.0$ 及相对 `R_med` 的选值系统偏移 $\Delta v_{\text{sel}} = v_{\text{sel}}(\text{R\_unit}) - v_{\text{sel}}(\text{R\_med})$（检验标准单位列范数 L1 是否导致显著 $v_{\text{sel}}$ 系统偏移，从而解释小 $Q$/大 $Q$ 偏差）。
- **`N3_DepthScaling`（子波提取深度单调标度律）**：在 `R_med` 与 `R_fix` 下，随 $t_w$ 从 $1.4 \to 2.0 \to 2.5\text{ s}$ 增大，三个剖面的 $\hat{\kappa}$ 均单调上升（$\hat{\kappa}(2.5) > \hat{\kappa}(2.0) > \hat{\kappa}(1.4)$），且与 $a(t_w)$ 标度一致。

## 3. 计算与资源约束

- 纯 Python (`D:\ChatGPT\研究\.venv\Scripts\python.exe`)，每 worker `OPENBLAS_NUM_THREADS=1`，内存 $< 250\text{ MB/worker}$，可断点续跑（`results/<job>.json` + `.claim`）。
