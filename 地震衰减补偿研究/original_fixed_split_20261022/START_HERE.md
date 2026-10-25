# 最初固定参数分步补偿：v1 原始归档

这里保留最早保存的 attenuation_rebuild_v1.zip 中全部13个文件，源码逐字节不改。原 ZIP 另存于 original_archive/，来源和文件 SHA256 见 PROVENANCE.json、PACKAGE_MANIFEST.json。

这是一份依据毕设第3章思路建立的研究重构程序。毕设原代码未提供，因此不能称为本科原程序或原图的逐点复现。后续 Frozen、CalShrink025Light、联合正则解耦与物理Q惩罚都没有加入这里的源码。

## 在 MATLAB 中运行

解压后进入 attenuation_rebuild 文件夹，执行：

```matlab
run_rebuild
```

无需额外 MATLAB 工具箱。参数统一在 config_model.m 中修改，结果生成在同目录 results/ 中；重复运行会覆盖同名结果。本次打包没有重新运行 MATLAB。

## 你要找的最初方法：FixedSplit

原入口同时输出固定分步及若干早期对照。你的目标分支名称是 **FixedSplit**，两个正则化参数均固定：

```matlab
c.lambda_x = 1e-3;       % 第一步：炮检距附加吸收补偿
c.lambda_fixed = 1e-3;   % 第二步：深度吸收补偿
```

算法对应原 run_rebuild.m 中：

```matlab
Bx = ridge_operator(Gx,c.lambda_x);
y = Bx*d;
B0 = ridge_operator(G0,c.lambda_fixed);
fixed = B0*y;
```

其中 Gx 表示炮检距附加吸收算子，G0 表示垂直路径的深度吸收算子。ridge_operator 实现 Tikhonov 正则化逆。

运行完成后，只提取固定分步结果：

```matlab
load(fullfile('results','experiment.mat'),'c','D','names','T');
idx = find(strcmp(names,'FixedSplit'));
fixed_result = D(:,:,idx);  % 时间采样点 × 炮检距道
fixed_metrics = T(T.MethodIndex == idx,:);
```

其他名称中，Truth 是合成无吸收参考，Observed 是观测，Step1 是第一步输出；AdaptiveUnweighted、CovarianceSplit 和 Direct 是原归档附带的早期对照，不能混称为 FixedSplit。第一步的评分目标是 G0*s，完整补偿的评分目标才是 s。

## 默认合成参数

- 采样间隔2 ms，601点，11道，炮检距0:150:1500 m。
- 速度 [2100,2250,2450,2800] m/s；Q=[40,50,70,100]。
- 事件时间0.38、0.65、0.90 s；主频60 Hz；合成信噪比18 dB。
- 随机种子20261017；速度、Q等设置及若干自选参数的边界见原 README.md。

本版本使用配置中的 Q 构造补偿算子，不含未知Q估计。附带的早期自适应对照使用合成噪声强度；这与后来从资料估计噪声的版本有区别。

## 文件与验证状态

- run_rebuild.m：原始主入口，FixedSplit 与早期对照。
- config_model.m：原始模型和固定参数。
- ray_absorption.m、attenuation_matrix.m、ridge_operator.m：射线、吸收矩阵和正则化逆。
- make_reference.m、event_metrics.m：合成参考与事后评分。
- verify_numpy.py、verification/numpy_checks.json：当时保存的独立 NumPy 核验及结果，不等于执行了同一份 MATLAB 源码。
- 原 README.md、VALIDATION.md：原样保留，记录该归档形成时的实现范围与执行状态。

原归档的 VALIDATION.md 记载当时 MATLAB 启动未成功，因此这份原 ZIP 没有 MATLAB 运行产物。之后历史工程虽有成功 MATLAB 验证，本包不把后来的修改版结果混入最初快照，也不把历史状态改写成本次验证通过。本次已检查 ZIP 完整性与解压源码哈希。
