# HVGDecision 0.9.1 中文使用说明

HVGDecision 是单细胞 RNA-seq integration 之前的 **integration-aware HVG refinement** 工具。

> **0.9.x 最重要的变化：用户必须显式选择 `mode`。**
>
> - `within_domain`：同一组织/相近实验体系，多 donor 或 batch。
> - `cross_domain`：Reference 与 Query 来自不同 dataset、technology 或明显不同的 domain。

HVGDecision 不是 Harmony、Seurat CCA、BBKNN、Scanorama 或 scVI 的替代品。它先决定用于这些下游整合方法的 feature panel。

---

# 1. 用户安装：直接安装 WHL

解压下载的 ZIP 后，进入 `HVGDecision_v0.9.1` 目录。普通用户直接安装根目录里的 wheel：

```bash
python -m pip install ./hvgdecision-0.9.1-py3-none-any.whl
```

如果 wheel 位于 `dist/`：

```bash
python -m pip install ./dist/hvgdecision-0.9.1-py3-none-any.whl
```

检查是否成功：

```bash
python -c "import hvgdecision as hd; print(hd.__version__); print(hd.VALID_MODES)"
```

正常应看到：

```text
0.9.1
('within_domain', 'cross_domain')
```

升级同一环境中的旧版本：

```bash
python -m pip install --upgrade ./hvgdecision-0.9.1-py3-none-any.whl
```

`pip install -e .` 是 **editable/development install**，用于开发者直接修改源码，不是普通用户的主安装方式。开发者如需使用：

```bash
python -m pip install -e ".[dev]"
```

---

# 2. 第一步：不要手动猜 raw counts，先让工具自动找

HVGDecision 默认：

```python
counts_layer="auto"
```

它会自动审计 AnnData 中可能的原始整数 counts，并按优先级寻找：

```text
常见 counts layer 名
    counts
    raw_counts
    raw.counts
    rawcounts
    umi_counts
    count
        ↓
adata.raw.X
        ↓
adata.X
        ↓
其它名称的 adata.layers[...]（也会逐一审计）
```

工具不是只看 layer 名；还会检查矩阵是否：

- 非负；
- 基本为整数型 count；
- 不是已经 Normalize/Log1p 后的小数表达矩阵。

## 2.1 推荐先检查一次

```python
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("your_data.h5ad")

count_check = hd.find_raw_counts(adata)

print("valid:", count_check.valid)
print("selected:", count_check.location)
print(count_check.audit)
```

例如可能输出：

```text
valid: True
selected: adata.layers['counts']
```

正常使用时无需再手动指定：

```python
study = hd.setup_reference_query(
    adata,
    mode="within_domain",
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2"],
    query=["D3"],
    counts_layer="auto",   # 默认就是 auto，可省略
)
```

每次 `setup_reference_query()` 都会把完整审计保存到：

```text
HVGDecision_results/.../raw_count_source_audit.csv
```

## 2.2 只有自动检测失败时才手动指定

如果你明确知道 raw counts 在：

```python
adata.layers["RNA_counts"]
```

可以写：

```python
counts_layer="RNA_counts"
```

如果 raw counts 在：

```python
adata.raw.X
```

写：

```python
counts_layer="raw"
```

如果你确认 `adata.X` 就是未经归一化的整数 counts：

```python
counts_layer=None
```

---

# 3. 必须先选择 mode

## 3.1 `within_domain`

适合：

- 同一组织；
- 同一或相近测序体系；
- 多 donor / batch；
- 目标是找出在控制 biological identity 后仍携带 donor/batch 风险的 HVG。

最终 raw risk：

$$
R_g^{within}=Z(L_g)+0.75Z(I_g)-Z(B_g)
$$

最终 harmful gate：

$$
H_g^{within}=
\mathbf{1}\left[
q_g\le0.05\ \land\
Z(L_g)\ge1\ \land\
Z_R(g)\ge1\ \land\
Z(B_g)\le0\ \land\
F_g^{boot}\ge0.80\ \land\
\neg P_g
\right]
$$

其中：

- $L_g$：biology residualization 后的 donor leakage；
- $I_g$：donor × biology interaction instability；
- $B_g$：biological explained variance；
- $q_g$：biology strata 内 donor-label permutation 后的 BH-FDR；
- $F_g^{boot}$：bootstrap recurrence；
- $P_g$：marker / explicit biological protection。

---

## 3.2 `cross_domain`

适合：

- cross-dataset；
- cross-technology；
- Reference 与 Query 不是同一个实验 domain。

Cross-domain 模式修正的是 **Query HVG panel**。

Reference risk percentile：

$$
R_g\in[0,1]
$$

Reference→Query mean-expression shift：

$$
E_g=
\frac{|\mu_{Q,g}-\mu_{R,g}|}
{\sqrt{(\mathrm{Var}_{R,g}+\mathrm{Var}_{Q,g})/2}+\epsilon}
$$

Detection-rate shift：

$$
D_g=|\pi_{Q,g}-\pi_{R,g}|
$$

Dataset/domain shift score：

$$
S_g=0.70E_g^*+0.30D_g^*
$$

Reference biology protection：

$$
B_g=
0.60B_g^{celltype}
+0.25M_g
+0.15P_g
$$

最终 Cross-domain Rule V3：

$$
R_g^{cross}
=
\min(R_g,S_g)
\left[
0.75+0.25\left(1-\left|R_g-S_g\right|\right)
\right]
\left(1-0.75B_g\right)
$$

然后从 Query HVG panel 中删除最高风险、且未被 hard protection 的 top-$k$ genes：

$$
\mathcal{H}_{final}^{Q}
=
\mathcal{H}_{base}^{Q}
\setminus
\operatorname{TopK}(R_g^{cross})
$$

Query 的真实 cell-type labels **不会用于 Cross-domain feature selection**。

---

# 4. Within-domain 完整教程

假设数据里：

```text
adata.obs['donor']
adata.obs['cell_type']
```

Reference donors = D1–D6，Query donors = D7–D8。

```python
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("PBMC.h5ad")

# 可选：先看自动找到哪一个 raw-count source
check = hd.find_raw_counts(adata)
print(check.location)
print(check.audit)

study = hd.setup_reference_query(
    adata,
    mode="within_domain",          # 必须显式选择
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2", "D3", "D4", "D5", "D6"],
    query=["D7", "D8"],
    counts_layer="auto",          # 推荐；自动寻找 raw counts
    dataset_name="PBMC",
    output_dir="HVGDecision_results/PBMC",
)

result = study.run(
    hvg_method="seurat_v3",
    return_details=True,
)

print("recommended base HVG N:", result.recommended_n_hvg)
print("harmful genes:", result.harmful_genes)
print("final HVG N:", result.final_n_hvg)
print("first final genes:", result.final_hvg_genes[:20])
```

`within_domain` 的 `study.run()` 会先执行 Reference-only minimum-sufficient HVG budget search，再执行 harmful-gene refinement。

如只需要预算推荐：

```python
result = study.find_best_hvg(return_details=True)
```

`find_best_hvg()` 只用于 `within_domain`。

---

# 5. Cross-domain：两个独立 h5ad 的完整教程

这是最推荐的跨数据集用法。

假设：

```text
Reference = reference_10x.h5ad
Query     = query_smartseq2.h5ad
```

两边 donor/batch 列最好统一成同一个名字，例如：

```text
donor
```

Reference 需要真实 biological label，例如：

```text
cell_type
```

Query 可以没有真实 cell-type label；如果有，也只用于你之后自己的 benchmark，HVGDecision 的 Cross-domain risk 不读取 Query true labels。

## 5.1 分别读取两个数据集

```python
import scanpy as sc
import hvgdecision as hd

ref = sc.read_h5ad("reference_10x.h5ad")
qry = sc.read_h5ad("query_smartseq2.h5ad")
```

## 5.2 工具分别自动寻找两个数据集的 raw counts

```python
print(hd.find_raw_counts(ref).location)
print(hd.find_raw_counts(qry).location)
```

Reference 和 Query 的 raw counts **不要求必须放在相同 layer 名**。

例如可以是：

```text
Reference -> adata.layers['counts']
Query     -> adata.raw.X
```

## 5.3 用工具自动准备 cross-domain 输入

```python
combined, input_audit = hd.prepare_cross_domain_inputs(
    ref,
    qry,
    reference_counts="auto",
    query_counts="auto",
    domain_key="hvgdecision_domain",
    reference_name="reference",
    query_name="query",
    return_audit=True,
)

print(input_audit)
print(combined)
```

这个函数会自动：

```text
Reference raw-count audit
Query raw-count audit
        ↓
分别读取真实 raw counts
        ↓
匹配共同 gene IDs
        ↓
只保留 Reference / Query 都存在的共同基因
        ↓
合并成一个 AnnData
        ↓
combined.layers['counts'] = raw counts
combined.obs['hvgdecision_domain'] = reference / query
```

因此不需要用户自己手工拼 counts matrix。

## 5.4 建立 Cross-domain study

```python
study = hd.setup_reference_query(
    combined,
    mode="cross_domain",             # 必须显式选择
    batch_key="donor",
    label_key="cell_type",          # Reference label
    split_key="hvgdecision_domain",
    reference=["reference"],
    query=["query"],
    counts_layer="auto",
    dataset_name="Skin_10X_to_SS2",
    output_dir="HVGDecision_results/Skin_cross_domain",
)
```

## 5.5 运行 Cross-domain Rule V3

```python
result = study.run(
    n_hvg=2000,
    hvg_method="seurat_v3",
    cross_domain_delete_budget=5,
    return_details=True,
)

print("removed genes:", result.harmful_genes)
print("final HVG N:", result.final_n_hvg)
print(result.final_hvg_genes[:20])
```

`cross_domain_delete_budget` 是删除深度，不是 mode：

```text
mode='cross_domain'
    决定使用 Cross-domain Rule V3

cross_domain_delete_budget=5
    决定从 Query HVG 中删除 top-5 eligible risk genes
```

论文主分析可使用保守的 `k=5`；软件仍允许显式指定其它非负整数。

---

# 6. 已经有本地 HVG 表：直接导入，不重新选择 HVG

HVGDecision 支持把 Seurat、Scanpy 或其它工具已经得到的 HVG 当作冻结 base panel。

最简单方式：**直接把本地文件路径传给 `study.refine_hvg()`**。

## 6.1 CSV：推荐格式

`my_hvg.csv`：

```text
gene
IL7R
LTB
MALAT1
CCR7
...
```

运行：

```python
result = study.refine_hvg(
    "my_hvg.csv",
    method_name="Seurat_v3",
    initial_n_hvg=2000,
    return_details=True,
)
```

## 6.2 TSV

`my_hvg.tsv`：

```text
gene
IL7R
LTB
MALAT1
CCR7
```

同样直接：

```python
result = study.refine_hvg(
    "my_hvg.tsv",
    method_name="Seurat_v3",
    initial_n_hvg=2000,
    return_details=True,
)
```

## 6.3 TXT：一行一个 gene，可以没有表头

`my_hvg.txt`：

```text
IL7R
LTB
MALAT1
CCR7
```

0.9.1 会自动按单列 gene list 读取，不会把第一条 gene 误当成表头：

```python
result = study.refine_hvg(
    "my_hvg.txt",
    method_name="external",
    initial_n_hvg=2000,
    return_details=True,
)
```

## 6.4 有 rank / budget / method 的完整表

也支持：

```text
method,n_hvg,gene_rank,gene
Seurat_v3,2000,1,IL7R
Seurat_v3,2000,2,LTB
Seurat_v3,2000,3,MALAT1
...
```

然后：

```python
result = study.refine_hvg(
    "hvg_full_table.csv",
    method_name="Seurat_v3",
    gene_column="gene",
    rank_column="gene_rank",
    budget_column="n_hvg",
    method_column="method",
    initial_n_hvg=2000,
    return_details=True,
)
```

### mode 对外部 HVG 的含义

```text
mode='within_domain'
    外部 HVG = 冻结的 within-domain base panel

mode='cross_domain'
    外部 HVG = 冻结的 Query base HVG panel
```

也就是说，Cross-domain 下本地上传的 HVG 应该是 **Query HVG**。

---

# 7. 最终结果怎么取

```python
adata_final = result.adata

final_hvgs = adata_final.var_names[
    adata_final.var["highly_variable"]
].tolist()

harmful = adata_final.var_names[
    adata_final.var["hvgdecision_harmful"]
].tolist()

print(len(final_hvgs))
print(harmful)
```

HVGDecision 不会把 AnnData 裁剪到只剩 HVG。最终对象仍保留 raw-count source 中的完整 gene space；`highly_variable=True` 标记最终面板。

核心 `var` 字段：

```text
highly_variable
hvgdecision_candidate
hvgdecision_risk_score
hvgdecision_risk_flagged
hvgdecision_marker_protected
hvgdecision_harmful
hvgdecision_removed
hvgdecision_final
hvgdecision_reason
```

---

# 8. 主要输出文件

所有模式都会输出 raw-count source 审计：

```text
raw_count_source_audit.csv
```

Within-domain 主要输出：

```text
00_budget_search/
01_hvg_decision/
    hvg_decision_table.csv
    harmful_genes.csv
    final_hvg_genes.csv
    decision.json
```

Cross-domain 主要输出包括：

```text
REFERENCE_RULE_V1_full_audit_for_crossdomain.csv
CROSSDOMAIN_reference_query_shift_audit.csv
CROSSDOMAIN_RULE_V3_ranking.csv
CROSSDOMAIN_RULE_V3_budget_gene_membership.csv
CROSSDOMAIN_RULE_V3_selection_manifest.csv
CROSSDOMAIN_RULE_V3_removal_audit.csv
```

---

# 9. CLI

CLI 使用已经合并好的 h5ad。

Within-domain：

```bash
hvgdecision refine \
  --config configs/within_domain_example.yml \
  --mode within_domain
```

Cross-domain：

```bash
hvgdecision refine \
  --config configs/cross_domain_example.yml \
  --mode cross_domain \
  --delete-budget 5
```

YAML 中也可以写 `mode:`。如果 CLI 和 YAML 同时提供，CLI 覆盖 YAML。

---

# 10. 可选：论文 benchmark 的三维综合评分

该评分只用于 benchmark，不参与 feature selection。

三个维度：

```text
Transfer
    Macro-F1
    Balanced accuracy
    Rare-cell Macro-F1

Batch correction
    Dataset mixing
    Donor mixing within cell type

Biological fidelity
    Leiden ARI
    Leiden NMI
    Cell-type silhouette
```

固定 `dataset × integration_method` 内先对 feature panels 做 min-max：

$$
\tilde{x}=\frac{x-x_{min}}{x_{max}-x_{min}}
$$

然后：

$$
S_{Transfer}
=\frac{
\widetilde{MacroF1}
+\widetilde{BalancedAccuracy}
+\widetilde{RareCellMacroF1}
}{3}
$$

$$
S_{Batch}
=\frac{
\widetilde{DatasetMixing}
+\widetilde{DonorMixingWithinCellType}
}{2}
$$

$$
S_{Bio}
=\frac{
\widetilde{ARI}
+\widetilde{NMI}
+\widetilde{CellTypeSilhouette}
}{3}
$$

默认 scIB-inspired overall：

$$
S_{Overall}
=0.30S_{Transfer}
+0.40S_{Batch}
+0.30S_{Bio}
$$

正式论文中建议称为：

```text
scIB-inspired three-domain composite score
```

而不是声称它是原始 scIB score 的完全复现。

---

# 11. Examples

ZIP 中附带：

```text
examples/00_raw_counts_auto.py
examples/01_within_domain.py
examples/02_cross_domain_two_h5ad.py
examples/03_local_hvg_table.py
```

建议第一次使用按 00 → 01/02 → 03 的顺序阅读。
