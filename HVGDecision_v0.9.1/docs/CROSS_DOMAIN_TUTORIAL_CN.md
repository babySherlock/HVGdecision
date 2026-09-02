# Cross-domain 完整教程：Reference h5ad → Query h5ad

本教程适用于：

```text
不同 dataset
不同测序 technology
同一组织但明显不同实验 domain
```

示例：

```text
Reference: Skin 10X
Query:     Skin Smart-seq2
```

## 1. 读取两个独立 AnnData

```python
import scanpy as sc
import hvgdecision as hd

ref = sc.read_h5ad("reference_10x.h5ad")
qry = sc.read_h5ad("query_smartseq2.h5ad")
```

## 2. 确认 metadata

Reference 至少需要：

```text
batch/donor column
biological label column
```

例如：

```python
print(ref.obs[["donor", "cell_type"]].head())
```

Query 至少需要 batch/donor column：

```python
print(qry.obs[["donor"]].head())
```

如果两边 batch 列名字不同，先统一：

```python
qry.obs["donor"] = qry.obs["sample_id"].astype(str)
```

Query 没有 `cell_type` 也没关系；合并后该列在 Query cells 中可以是 missing。Cross-domain feature selection 不读取 Query true cell-type labels。

## 3. 自动寻找两边 raw counts

```python
ref_check = hd.find_raw_counts(ref)
qry_check = hd.find_raw_counts(qry)

print(ref_check.location)
print(qry_check.location)
```

可以分别位于不同位置，例如：

```text
Reference: adata.layers['counts']
Query:     adata.raw.X
```

## 4. 自动构建共同 gene space

```python
combined, audit = hd.prepare_cross_domain_inputs(
    ref,
    qry,
    reference_counts="auto",
    query_counts="auto",
    domain_key="hvgdecision_domain",
    reference_name="reference",
    query_name="query",
    return_audit=True,
)

print(audit)
print(combined)
```

此步骤会：

1. 分别验证 Reference/Query raw counts；
2. 获取两边 raw-count gene IDs；
3. 取共同基因；
4. 按共同基因对齐矩阵；
5. 合并 cells；
6. 将 raw counts 放入 `combined.X` 和 `combined.layers['counts']`；
7. 建立 `combined.obs['hvgdecision_domain']`。

## 5. 建立 study

```python
study = hd.setup_reference_query(
    combined,
    mode="cross_domain",
    batch_key="donor",
    label_key="cell_type",
    split_key="hvgdecision_domain",
    reference=["reference"],
    query=["query"],
    counts_layer="auto",
    output_dir="HVGDecision_results/Skin_cross_domain",
)
```

## 6A. 让 HVGDecision 在 Query 内部选择 HVG2000

```python
result = study.run(
    n_hvg=2000,
    hvg_method="seurat_v3",
    cross_domain_delete_budget=5,
    return_details=True,
)
```

## 6B. 如果你已经有 Query HVG 表

例如：

```text
query_hvg.csv
```

内容：

```text
gene
GENE1
GENE2
GENE3
...
```

直接：

```python
result = study.refine_hvg(
    "query_hvg.csv",
    method_name="Seurat_v3",
    initial_n_hvg=2000,
    return_details=True,
)
```

Cross-domain 下，这个外部 HVG 文件被解释为 **Query base panel**，不会重新跑 HVG selection。

## 7. Rule V3

$$
S_g=0.70E_g^*+0.30D_g^*
$$

$$
B_g=0.60B_g^{celltype}+0.25M_g+0.15P_g
$$

$$
R_g^{cross}
=
\min(R_g,S_g)
\left[
0.75+0.25\left(1-\left|R_g-S_g\right|\right)
\right]
\left(1-0.75B_g\right)
$$

## 8. 最终输出

```python
print(result.harmful_genes)
print(result.final_hvg_genes)
print(result.final_n_hvg)
```

主要审计文件：

```text
REFERENCE_RULE_V1_full_audit_for_crossdomain.csv
CROSSDOMAIN_reference_query_shift_audit.csv
CROSSDOMAIN_RULE_V3_ranking.csv
CROSSDOMAIN_RULE_V3_budget_gene_membership.csv
CROSSDOMAIN_RULE_V3_selection_manifest.csv
CROSSDOMAIN_RULE_V3_removal_audit.csv
```
