# HVGDecision 0.9.1 中文使用说明

**用于单细胞整合之前的 integration-aware HVG refinement。**

HVGDecision 是一个上游特征选择/特征修正工具。它不是 Seurat CCA、Harmony、BBKNN、Scanorama 或 scVI 的替代品，而是在这些整合方法之前，对将要用于整合的 HVG panel 进行筛选与修正。

0.9.1 版本要求用户必须显式选择一种分析模式：

```text
within_domain
cross_domain
```

- `within_domain`：Reference 与 Query 属于同一或相近实验体系，例如同一组织、多 donor、多 batch。
- `cross_domain`：Reference 与 Query 来自不同 dataset、不同技术或明显不同的实验 domain。

---

# 1. 普通用户安装：直接安装 WHL

下载并解压 HVGDecision ZIP 后，进入包含 wheel 文件的目录。

普通用户直接运行：

```bash
python -m pip install ./hvgdecision-0.9.1-py3-none-any.whl
```

如果 wheel 位于 `dist/` 目录：

```bash
python -m pip install ./dist/hvgdecision-0.9.1-py3-none-any.whl
```

如果环境中已经安装旧版本：

```bash
python -m pip install --upgrade ./hvgdecision-0.9.1-py3-none-any.whl
```

检查是否安装成功：

```bash
python -c "import hvgdecision as hd; print(hd.__version__); print(hd.VALID_MODES)"
```

正常应看到：

```text
0.9.1
('within_domain', 'cross_domain')
```

`pip install -e .` 是 editable/development install，主要用于开发者修改源码后立即生效，**不是普通用户的推荐安装方式**。

---

# 2. 第一步：让工具自动寻找 raw counts

HVGDecision 默认使用：

```python
counts_layer="auto"
```

因此用户不需要一开始就自己猜 raw counts 在哪里。

工具会优先审计常见 counts layer，例如：

```text
counts
raw_counts
raw.counts
rawcounts
umi_counts
count
```

同时还会检查：

```text
adata.raw.X
adata.X
其它 adata.layers[...]
```

HVGDecision 不只是检查 layer 名字，还会检查候选矩阵是否具有 raw-count 特征，例如：

- 非负；
- 基本为整数型 count；
- 不像已经 Normalize/Log1p 后的小数表达矩阵。

## 2.1 推荐先查看一次自动检测结果

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

正常分析时通常不需要再手工指定：

```python
study = hd.setup_reference_query(
    adata,
    mode="within_domain",
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2"],
    query=["D3"],
    counts_layer="auto",
)
```

完整的 raw-count source 审计会保存为：

```text
raw_count_source_audit.csv
```

## 2.2 只有自动检测失败或你明确想指定某个来源时才手工设置

如果 raw counts 明确在：

```python
adata.layers["RNA_counts"]
```

则写：

```python
counts_layer="RNA_counts"
```

如果 raw counts 在：

```python
adata.raw.X
```

则写：

```python
counts_layer="raw"
```

如果你确认：

```python
adata.X
```

本身就是未经归一化的整数 counts，则写：

```python
counts_layer=None
```

---

# 3. 必须先选择 mode

每次分析都必须显式写 `mode`。

## 3.1 `within_domain`

适合：

- 同一组织；
- 同一或相近测序体系；
- 多 donor / batch；
- donor-held-out 或 batch-held-out integration 场景。

目标是在控制 biological identity 后，找出仍然携带 donor/batch 风险的 HVG。

## 3.2 `cross_domain`

适合：

- cross-dataset；
- cross-technology；
- 10X Reference → Smart-seq2 Query；
- 不同实验平台或协议。

Cross-domain 模式修正的是 **Query HVG panel**。

---

# 4. Within-domain 方法

## 4.1 Biological explained variance

对于基因 `g`：

```math
B_g
=
\frac{
\sum_c n_c\left(\bar{x}_{cg}-\bar{x}_g\right)^2
}{
\sum_i\left(x_{ig}-\bar{x}_g\right)^2+\epsilon
}
```

## 4.2 在 biological strata 内进行残差化

```math
r_{ig}
=
x_{ig}-\bar{x}_{c_i g}
```

## 4.3 Donor leakage

```math
L_g
=
\frac{
\sum_d n_d\bar{r}_{dg}^{2}
}{
\sum_i r_{ig}^{2}+\epsilon
}
```

## 4.4 Donor × biology interaction instability

对于每个满足条件的 biological stratum `c`：

```math
I_{gc}
=
\frac{
\max_d \bar{x}_{cdg}-\min_d \bar{x}_{cdg}
}{
s_{cg}+\epsilon
}
```

随后得到 gene-level interaction score：

```math
I_g
=
\mathrm{median}_c\left(I_{gc}\right)
```

## 4.5 Within-domain risk

对三个分量分别在基因间进行 robust standardization：

```math
R_g^{\mathrm{within}}
=
Z\left(L_g\right)
+
0.75 Z\left(I_g\right)
-
Z\left(B_g\right)
```

随后再次对 raw risk 做 robust standardization，得到 `Z_R(g)`。

## 4.6 Conditional permutation

在 biological strata 内打乱 donor labels：

```math
p_g
=
\frac{
1+\sum_{p=1}^{P}
\mathbf{1}\left(L_g^{(p)}\ge L_g\right)
}{
P+1
}
```

随后进行 Benjamini-Hochberg FDR 校正得到 `q_g`。

默认：

```text
P = 100
FDR threshold = 0.05
```

## 4.7 Bootstrap recurrence

每次 bootstrap 中，必须同时满足：

```math
Z\left(L_g^{(b)}\right)\ge 1
```

```math
Z_R^{(b)}(g)\ge 1
```

```math
Z\left(B_g^{(b)}\right)\le 0
```

Bootstrap recurrence fraction：

```math
F_g^{\mathrm{boot}}
=
\frac{1}{N_{\mathrm{boot}}}
\sum_b
\mathbf{1}\left[\mathrm{pass}_b\right]
```

默认：

```text
N_boot = 20
minimum recurrence = 0.80
```

## 4.8 最终 Within-domain harmful gate

```math
H_g^{\mathrm{within}}
=
\mathbf{1}\left[
q_g\le 0.05
\;\land\;
Z(L_g)\ge 1
\;\land\;
Z_R(g)\ge 1
\;\land\;
Z(B_g)\le 0
\;\land\;
F_g^{\mathrm{boot}}\ge 0.80
\;\land\;
\neg P_g
\right]
```

其中 `P_g` 表示 replicated marker protection 或 explicit hard protection。

最终 panel：

```math
\mathcal{H}_{\mathrm{final}}
=
\mathcal{H}_{\mathrm{base}}
\setminus
\left\{g:H_g^{\mathrm{within}}=1\right\}
```

**删除 0 个基因也是合法结果。**

---

# 5. Within-domain 完整教程

假设：

```text
adata.obs['donor']
adata.obs['cell_type']
```

Reference donors = D1-D6，Query donors = D7-D8。

```python
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("PBMC.h5ad")

# 可选：先查看自动找到的 raw-count source
check = hd.find_raw_counts(adata)
print(check.location)
print(check.audit)

study = hd.setup_reference_query(
    adata,
    mode="within_domain",
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2", "D3", "D4", "D5", "D6"],
    query=["D7", "D8"],
    counts_layer="auto",
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

`within_domain` 下的 `study.run()` 会先执行 Reference-only minimum-sufficient HVG budget search，再执行 harmful-gene refinement。

如果只想显式运行预算搜索：

```python
budget_result = study.find_best_hvg(return_details=True)
```

`find_best_hvg()` 只适用于 `within_domain`。

---

# 6. Cross-domain 方法

Cross-domain 模式修正的是 **Query HVG panel**。

它使用：

- Reference technical risk；
- 无需 Query true labels 的 Reference→Query distribution shift；
- Reference biological protection。

**Query 的真实 cell-type labels 不会用于 Cross-domain feature selection。**

## 6.1 Reference technical-risk percentile

Reference Rule risk 转换为 percentile：

```math
R_g\in[0,1]
```

## 6.2 Reference→Query mean-expression shift

```math
E_g
=
\frac{
\left|\mu_{Q,g}-\mu_{R,g}\right|
}{
\sqrt{
\left(\mathrm{Var}_{R,g}+\mathrm{Var}_{Q,g}\right)/2
}
+
\epsilon
}
```

## 6.3 Detection-rate shift

```math
D_g
=
\left|\pi_{Q,g}-\pi_{R,g}\right|
```

随后分别转成 percentile：

```math
E_g^{*}
=
\mathrm{Pct}\left(E_g\right)
```

```math
D_g^{*}
=
\mathrm{Pct}\left(D_g\right)
```

## 6.4 Dataset/domain shift score

```math
S_g
=
0.70 E_g^{*}
+
0.30 D_g^{*}
```

## 6.5 Reference biology protection

```math
B_g
=
0.60 B_g^{\mathrm{celltype}}
+
0.25 M_g
+
0.15 P_g
```

其中：

- `B_g^{celltype}`：Reference cell-type biological effect 的 percentile score；
- `M_g`：eligible Reference donors 中 marker replication fraction；
- `P_g`：explicit/hard protection indicator。

Hard-protected genes 不进入可删除候选集合。

## 6.6 Dual technical consensus

```math
T_g
=
\min\left(R_g,S_g\right)
```

两个 technical evidence source 的一致性：

```math
A_g
=
1-
\left|R_g-S_g\right|
```

## 6.7 最终 Cross-domain Rule V3

```math
R_g^{\mathrm{cross}}
=
\min\left(R_g,S_g\right)
\left[
0.75
+
0.25
\left(
1-
\left|R_g-S_g\right|
\right)
\right]
\left(
1-0.75B_g
\right)
```

这个公式要求 Reference technical risk 和 Reference→Query shift 两条技术证据同时较高；两者越一致，risk 越高；如果一个基因具有较强 Reference biological support，则删除优先级会下降。

## 6.8 最终 Query HVG panel

原始 Query HVG panel：

```math
\mathcal{H}_{\mathrm{base}}^{Q}
```

排除 hard-protected genes 后，根据 `R_g^{cross}` 排序，删除最高风险的 top-`k` eligible genes：

```math
\mathcal{H}_{\mathrm{final}}^{Q}
=
\mathcal{H}_{\mathrm{base}}^{Q}
\setminus
\mathrm{TopK}
\left(
R_g^{\mathrm{cross}}
\right)
```

删除深度通过：

```python
cross_domain_delete_budget=5
```

进行控制。

`mode` 和 `delete_budget` 是两个不同概念：

```text
mode="cross_domain"
    决定使用 Cross-domain Rule V3

cross_domain_delete_budget=5
    决定从 Query HVG panel 中删除 top-5 eligible risk genes
```

---

# 7. Cross-domain：两个独立 h5ad 的完整教程

这是跨数据集/跨技术分析最推荐的使用方式。

例如：

```text
Reference = reference_10x.h5ad
Query     = query_smartseq2.h5ad
```

## 7.1 分别读取两个 AnnData

```python
import scanpy as sc
import hvgdecision as hd

ref = sc.read_h5ad("reference_10x.h5ad")
qry = sc.read_h5ad("query_smartseq2.h5ad")
```

Reference 至少应该具有：

```text
batch/donor column
biological label column
```

例如：

```text
donor
cell_type
```

Query 至少需要 batch/donor column。

Query 可以没有真实 cell-type label；即使存在，Cross-domain feature selection 也不会读取 Query true labels。

如果 Query 的 batch 列名字不同，可以先统一：

```python
qry.obs["donor"] = qry.obs["sample_id"].astype(str)
```

## 7.2 分别自动检测 Reference 和 Query 的 raw counts

```python
ref_check = hd.find_raw_counts(ref)
qry_check = hd.find_raw_counts(qry)

print("Reference:", ref_check.location)
print("Query:", qry_check.location)
```

Reference 和 Query 的 raw counts **不要求位于相同位置**。

例如完全可以是：

```text
Reference -> adata.layers['counts']
Query     -> adata.raw.X
```

## 7.3 自动对齐共同基因并构建 cross-domain 输入

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

该函数会自动：

1. 审计 Reference raw counts；
2. 审计 Query raw counts；
3. 读取两边通过验证的 raw-count matrix；
4. 寻找两边共有 gene IDs；
5. 按完全相同的 gene order 对齐 Reference/Query；
6. 合并 cells；
7. 将 raw counts 写入 `combined.layers['counts']`；
8. 创建 `combined.obs['hvgdecision_domain']`。

因此用户不需要自己手工拼接两个 counts matrix。

## 7.4 创建 Cross-domain study

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
    dataset_name="Skin_10X_to_SS2",
    output_dir="HVGDecision_results/Skin_cross_domain",
)
```

这里：

```text
label_key="cell_type"
```

指的是 **Reference biological label**。

Query cells 的 `cell_type` 可以为空；Cross-domain risk 不使用 Query true labels。

## 7.5 让 HVGDecision 自己在 Query 中计算 HVG2000

```python
result = study.run(
    n_hvg=2000,
    hvg_method="seurat_v3",
    cross_domain_delete_budget=5,
    return_details=True,
)

print("removed genes:", result.harmful_genes)
print("final HVG N:", result.final_n_hvg)
print("first final genes:", result.final_hvg_genes[:20])
```

---

# 8. 已经有本地 HVG 表：直接导入，不重新计算 HVG

HVGDecision 支持将 Seurat、Scanpy 或其他软件已经得到的 HVG panel 作为冻结 base panel 进行 refinement。

支持：

```text
.csv
.tsv
.txt
```

## 8.1 CSV：推荐包含 `gene` 列

例如 `my_hvg.csv`：

```text
gene
IL7R
LTB
MALAT1
CCR7
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

## 8.2 TSV

```python
result = study.refine_hvg(
    "my_hvg.tsv",
    method_name="Seurat_v3",
    initial_n_hvg=2000,
    return_details=True,
)
```

## 8.3 TXT：一行一个 gene，可以没有表头

例如 `my_hvg.txt`：

```text
IL7R
LTB
MALAT1
CCR7
```

直接运行：

```python
result = study.refine_hvg(
    "my_hvg.txt",
    method_name="external",
    initial_n_hvg=2000,
    return_details=True,
)
```

0.9.1 会按照单列 gene list 读取，第一条 gene 不会被误当成表头。

## 8.4 有 rank / budget / method 的完整表

也支持：

```text
method,n_hvg,gene_rank,gene
Seurat_v3,2000,1,IL7R
Seurat_v3,2000,2,LTB
Seurat_v3,2000,3,MALAT1
```

运行：

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

## 8.5 不同 mode 下，本地 HVG 的含义

```text
mode="within_domain"
    本地 HVG = 冻结的 within-domain base panel

mode="cross_domain"
    本地 HVG = 冻结的 Query base HVG panel
```

因此，Cross-domain 模式下上传/导入的 HVG 应该是 **Query HVG panel**。

---

# 9. 如何读取最终结果

```python
adata_final = result.adata

final_hvgs = adata_final.var_names[
    adata_final.var["highly_variable"]
].tolist()

harmful_genes = adata_final.var_names[
    adata_final.var["hvgdecision_harmful"]
].tolist()

print("final HVGs:", len(final_hvgs))
print("removed genes:", harmful_genes)
```

HVGDecision 不会把 AnnData 裁剪成只有 HVG。

最终对象仍然保留完整 gene space，最终 HVG panel 通过：

```text
adata.var['highly_variable'] == True
```

进行标记。

重要的 `var` 字段包括：

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

# 10. 主要输出文件

所有模式都会输出：

```text
raw_count_source_audit.csv
```

Within-domain 常见输出：

```text
00_budget_search/
01_hvg_decision/
    hvg_decision_table.csv
    harmful_genes.csv
    final_hvg_genes.csv
    decision.json
```

Cross-domain 常见输出：

```text
REFERENCE_RULE_V1_full_audit_for_crossdomain.csv
CROSSDOMAIN_reference_query_shift_audit.csv
CROSSDOMAIN_RULE_V3_ranking.csv
CROSSDOMAIN_RULE_V3_budget_gene_membership.csv
CROSSDOMAIN_RULE_V3_selection_manifest.csv
CROSSDOMAIN_RULE_V3_removal_audit.csv
```

---

# 11. CLI 使用

CLI 适用于已经准备/合并好的 h5ad。

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

YAML 中也可以写：

```yaml
mode: within_domain
```

或：

```yaml
mode: cross_domain
```

如果 CLI 和 YAML 同时指定 mode，则 CLI 参数覆盖 YAML。

---

# 12. 可选：论文 benchmark 的三维综合评分

这个评分只用于 benchmark evaluation，**不参与 feature selection**。

固定在同一个：

```text
dataset × integration_method
```

任务内，对不同 feature panels 的每个 metric 分别做 min-max normalization。

## 12.1 Transfer

```math
S_{\mathrm{Transfer}}
=
\frac{
\widetilde{\mathrm{MacroF1}}
+
\widetilde{\mathrm{BalancedAccuracy}}
+
\widetilde{\mathrm{RareCellMacroF1}}
}{3}
```

## 12.2 Batch correction

```math
S_{\mathrm{Batch}}
=
\frac{
\widetilde{\mathrm{DatasetMixing}}
+
\widetilde{\mathrm{DonorMixingWithinCellType}}
}{2}
```

## 12.3 Biological fidelity

```math
S_{\mathrm{Bio}}
=
\frac{
\widetilde{\mathrm{ARI}}
+
\widetilde{\mathrm{NMI}}
+
\widetilde{\mathrm{CellTypeSilhouette}}
}{3}
```

## 12.4 Overall

```math
S_{\mathrm{Overall}}
=
0.30 S_{\mathrm{Transfer}}
+
0.40 S_{\mathrm{Batch}}
+
0.30 S_{\mathrm{Bio}}
```

论文中应将其称为：

> **scIB-inspired three-domain composite score**

而不是 exact scIB score。

---

# 13. Query leakage 规则

## Within-domain

Feature-selection risk 使用：

- Reference expression；
- Reference donor/batch IDs；
- Reference biological labels。

Within-domain risk decision 不使用 Query expression 或 Query labels。

## Cross-domain

Feature selection 使用：

- Reference raw counts；
- Reference donor/batch IDs；
- Reference biological labels；
- Query raw counts；
- Query batch IDs，用于 Query HVG 构建和 Reference→Query shift estimation。

**Query 的真实 biological labels 不会用于 Cross-domain Rule V3 feature selection。**

Query true labels 只可以在 feature selection 和 downstream integration 全部完成后，用于独立 benchmark/evaluation。

---

# 14. 最小 API 总览

```python
import hvgdecision as hd

print(hd.__version__)
print(hd.VALID_MODES)

hd.find_raw_counts(...)
hd.prepare_cross_domain_inputs(...)
hd.setup_reference_query(...)
hd.three_domain_scores(...)
```

英文说明见 `README.md`。
