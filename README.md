# 研究工作进度归档

本仓库归档地震衰减补偿研究截至 **2026-10-23（Asia/Taipei）** 的本地工作进度，包括方法代码、实验配置、运行记录、指标、图件、研究报告、交付压缩包与参考资料。仓库原有的「研习」项目压缩包予以保留。

## 阅读入口

- [实验版本筛选与当前结论](EXPERIMENT_VERSION_SELECTION_20261023_ZH.md)：当前主版本选择、各轮结果、失败记录及验证边界。
- [联合正则研究方向](JOINT_REGULARIZATION_DIRECTIONS_20261022.md)。
- [最初研究交接](attenuation_research_handoff_20261019/START_HERE.md)。
- [自适应联合方法](adaptive_joint_cv_20261022/START_HERE.md)。
- [复杂地层实验](complex_fieldlike_20261023/START_HERE.md)。
- [去干扰与真实资料](denoise_and_field_20261023/START_HERE.md)。
- [实测资料试验](field_pilot_20261023/README.md)。
- [Q 误差稳定性](q_error_stability_20261023/README.md)。
- [归档文件清单](UPLOAD_MANIFEST_20261023.csv)：研究快照中每个文件的路径、字节数与 SHA-256。

各目录保留对应阶段的原始命名与结构；目录内的 README、START_HERE、RESULTS、PLAN 和日志记录具体进度。运行失败与未通过门槛的实验同样保留，上传归档不代表研究结果已获验收。

## 获取完整资料

较大的文件和交付压缩包使用 Git LFS 保存。完整下载请安装 Git LFS 后执行：

```bash
git lfs install
git clone https://github.com/q1754399106-sudo/work.git
cd work
git lfs pull
```

GitHub 的网页源码压缩包可能只包含 LFS 指针；需要实际数据时使用上面的命令。

## 归档范围

保存工作目录中的研究文件、结果和已有压缩包；排除可重建的 Python 虚拟环境、已安装的第三方依赖、字节码缓存、系统临时文件和本地凭证。排除规则见 [.gitignore](.gitignore)。原有压缩包按原样保存，没有重新打包。原始实验文件没有为上传而修改。

`q_dense_gathers_20261023` 归档时仍在运行，保存其暂存时刻的进度和日志快照。上传之后继续生成的内容需要后续再同步；本次没有中断实验。

这是文件快照，不包含 Codex 聊天记录。运行脚本可能仍需配置 MATLAB、本地路径及各阶段依赖；请先阅读对应目录说明。
