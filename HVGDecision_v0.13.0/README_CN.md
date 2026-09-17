# HVGDecision 0.13.0

本版在固定 HVG 面板中评估整合风险，结合 donor-level 和 cell-level 两类证据及生物学保护，输出保留基因的 AnnData 和逐基因审计表。两类证据都计算，按重复充分度 A 连续加权，不再按 0.75 阈值选择分支。

## 安装和运行

```bash
git clone https://github.com/babySherlock/HVGdecision.git
cd HVGdecision/HVGDecision_v0.13.0
python -m pip install ./dist/hvgdecision-0.13.0-py3-none-any.whl
```

安装到 notebook 所用环境，随后重启内核。导入仍为 `import hvgdecision as hd`。

- [PBMC notebook](examples/01_PBMC_reference_query.ipynb)：服务器路径、6 个 Reference donor、2 个 Query donor，Reference span=0.5，Query span=0.3。
- [Lung notebook](examples/02_HumanLung_pooled.ipynb)：pooled discovery，所有保留细胞参与发现，不属于独立 Query 验证。
- [完整用法](README.md)：指定 counts、外部面板、保护基因及输出结构。

先用 Scanpy 读入数据和筛选队列，再用 `hd.find_raw_counts` 验证 counts，最后调用 `hd.refine`。指定 Reference/Query 时，仅删除 Reference 风险名单与 Query HVG 的交集；Query cell-type 标签不参与风险发现。两组 donor 必须覆盖全部输入细胞。省略两组时，全体细胞用于 pooled discovery。

`n_hvg=2000` 是 Scanpy Seurat v3 基础面板大小，不是推荐数量搜索，也不保证删除固定数量的基因。可传入 R 导出的有序基因名单，跳过 HVG 拟合；缺失于 counts 的基因不会静默丢弃。

## 查看结果

`result.removed_genes` 是实际删除名单，`result.reference_risk_genes` 是发现队列的风险名单。完整证据见 `gene_decisions.csv`、`reference_gene_decisions.csv`、`donor_evidence.csv`、`cell_evidence.csv`。A、连续权重与可辨识性见 `design_audit.csv`。参数、源码哈希、依赖版本和面板哈希见 `run_manifest.json`。

`result.adata` 保留全部输入细胞和最终基因，X 与 counts layer 均为 raw counts。下游整合仍需各软件相应预处理，不保留旧 UMAP 或图结构。已有非空目录不会覆盖，允许删除 0 个基因。

## 复现与边界

复现旧实验应冻结 counts、细胞筛选、基因顺序、Reference/Query 面板、seed 和参数。不能强制删掉几个已知基因来充当算法复现。自动 span 只在数值拟合失败后重试，并保存尝试记录。

默认对不可辨识的支持设计发出警告并继续；可用 `design_policy='error'` 停止。A 高不等于设计可辨识，融合分数不是显著性或 FDR。细胞置换也不能代替独立 donor 重复。

本次发布测试不是服务器五种整合软件完整 benchmark 的重跑。详见 [测试记录](docs/TEST_REPORT.md) 和 [迁移说明](docs/MIGRATION.md)。

