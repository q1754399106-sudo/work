# GitHub 归档说明

本目录包含连续软收缩代码、完整实验输入输出、38 张高清图、报告和数值检查记录。查看 index.html 或 RESULTS_ZH.md。

按 D 盘主目录 AGENTS.md，cache/、Python 缓存、依赖和凭证保留在本地，不上传。PACKAGE_MANIFEST.json 与 DELIVERY_CHECK.json 是此前本地完整交付包的历史记录，其中列出的缓存不包含在 GitHub 工作树中。GitHub 专用 ZIP 使用独立的包内清单；实验数据和方法规则不变。

复算前先在 MATLAB 中创建空 cache 文件夹（若不存在），并保留同级 q_dense_gathers_20261023/results 作为已归档的比较基线。公开调用入口为 soft_refit_compensate.m。
