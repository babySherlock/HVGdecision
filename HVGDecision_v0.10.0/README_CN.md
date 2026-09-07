# HVGDecision 0.10.0

当前版本用于相近生物学背景下的多 donor 数据：检查 raw counts、按重复充分度与
可辨识性自动选择风险分支、筛除风险候选基因，最后返回保留特征的 AnnData。
默认固定 Seurat v3 batch-aware 2000 HVGs，不运行下游整合，也不搜索最佳 HVG 数量。

## 安装

在准备使用的 Python/conda 环境中，进入本文件夹后运行：

```bash
python -m pip install .
```

也可直接安装配套 wheel：

```bash
python -m pip install /path/to/hvgdecision-0.10.0-py3-none-any.whl
```

安装后重启 notebook kernel，检查实际导入版本：

```python
import sys
import hvgdecision as hd
print(sys.executable)
print(hd.__version__)
print(hd.__file__)
```

## 完整调用

```python
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad('/path/to/data.h5ad')
# 使用质控后的细胞，保留完整 raw-count 基因轴，不要先仅留下 2000 个基因。
# donor、cell_type 请改成当前 adata.obs 的真实列名。

raw_counts = hd.find_raw_counts(adata)
display(raw_counts.audit)
if not raw_counts.valid:
    raise ValueError(raw_counts.error)

result = hd.refine(
    adata,
    batch_key='donor',
    label_key='cell_type',
    counts=raw_counts,
    n_hvg=2000,
    output_dir='./HVGDecision_run01',  # 新文件夹，不覆盖旧结果
    return_details=True,
)

print(result)
print('实际分支：', result.selected_mode)
print('基础 HVG 数量：', result.base_n_hvg)
print('删除数量：', len(result.removed_genes))
print('删除基因：', result.removed_genes)
print('最终 HVG 数量：', result.final_n_hvg)
display(result.routing.audit)
display(result.decision_table)

adata_hvg = result.adata
```

只要 AnnData 时可简化成一句：

```python
adata_hvg = hd.refine(adata, batch_key='donor', label_key='cell_type')
```

## 先看会进入哪个分支

```python
design = hd.audit_design(adata, batch_key='donor', label_key='cell_type')
display(design.audit)
display(design.coverage)
```

可辨识且 `w >= 0.75` 时选择 `donor_aware`，否则选择 `cell_level`。
若设计不可辨识或支持不足，返回 `insufficient_confounded`，保留原面板，不强行删除。
例如 100 个 donor 各自只有不同的一种细胞，donor 与 cell type 无法区分；
100 个 donor 都只有同一种细胞，也不足以估计跨 cell type 的保护证据。
这两种情况不会仅因 donor 数多就开始去基因。

`w` 是固定规则计算的路由分数，不是 donor-aware 更优的概率。
路由默认 `tau=15`，`A` 转换区间为 2 到 4，判定界限为 0.75。
这些值没有为了复现 PBMC 删除 3 个或 Lung 删除 10 个而硬编码结果。

## counts 可以由用户指定

```python
raw_counts = hd.find_raw_counts(adata, source='raw')       # adata.raw.X
raw_counts = hd.find_raw_counts(adata, source='X')         # adata.X
raw_counts = hd.find_raw_counts(adata, source='raw.counts')  # 自定义 layer 名
raw_counts = hd.find_raw_counts(adata, source='/path/to/counts.csv')
```

外部 CSV/TSV 或 DataFrame 必须带细胞 ID 与基因 ID。工具按细胞 ID 对齐，
支持转置，不会只因行数相同就接受。裸矩阵可用 `(matrix, gene_names)`，
此时细胞顺序必须与 `adata.obs_names` 一致。更安全的形式：

```python
raw_counts = hd.find_raw_counts(adata, source={
    'matrix': my_counts,
    'gene_names': my_genes,
    'obs_names': my_cell_ids,
})
```

检查覆盖全部存储值，要求非负、有限、整数型数值，并在运行前拒绝零文库细胞。
整数型数值只能通过数值检查，不能证明其一定来自原始测序；不要把归一化值四舍五入后输入。
手动指定的来源不合格时会报错，不会换用其他来源。

## 使用 R 导出的 HVG 名单

```python
import pandas as pd
panel = pd.read_csv('/path/to/hvg2000.csv')['gene'].tolist()
result = hd.refine(
    adata, batch_key='donor', label_key='cell_type', counts=raw_counts,
    hvg_genes=panel, return_details=True,
)
```

名单按传入顺序记录 rank，直接绕过 Seurat v3。缺失基因、重复 ID 会报错，
不会自行补齐、截断、转换 gene symbol 或 Ensembl ID。
`protected_genes=[...]` 可手动保护基因，默认没有显式保护名单。

## Reference 的含义

默认使用全部输入细胞进行路由、HVG 选择和风险分析。
若传入 `reference=['D1', 'D2']`，这三步只使用指定 donor，输出仍保留全部细胞。
这不是自动完成了独立 Query 验证。如果某些 cell-type 标签已经用于特征筛选，
不能再把同批细胞上的分类成绩称为独立的 label-transfer 测试。

## 参数与输出

```python
cell_config = hd.CellLevelConfig()  # 100 permutations、20 bootstraps 等原始默认值
donor_config = hd.DonorAwareConfig()
route_config = hd.RoutingConfig(donor_aware_threshold=0.75)
```

三个对象可分别通过 `cell_level_config`、`donor_aware_config`、`routing_config` 传入。
显式指定 `mode='cell_level'` 或 `'donor_aware'` 仅用于敏感性分析，会记录 override，
且无法绕过不可辨识设计保护。不要为了得到预期的删除数调阈值。

最重要的文件是 `adata_hvg.h5ad`、`gene_decisions.csv`、`removed_genes.csv` 和
`routing_summary.csv`。完整参数、版本和哈希位于 `run_manifest.json`。
所有基因的原始证据也在 `adata_hvg.uns['hvgdecision']['gene_decisions']` 中。
已删除基因不在最终 `.var` 里，不能通过最终 `.var_names` 查它们。

最终 `.X` 与 `layers['counts']` 都是保留基因的 raw counts；原始 UMAP、PCA 和邻接图
不会复制。用于 Harmony、BBKNN、Scanorama、scVI、Seurat CCA 时，仍需运行各方法的预处理。
筛选与五方法效果评估是两步，筛选快不表示漏跑了整合。

原版本迁移见 [MIGRATION.md](docs/MIGRATION.md)，算法与全部阈值见
[METHODS.md](docs/METHODS.md)。本次测试范围见 [TESTING.md](docs/TESTING.md)。
