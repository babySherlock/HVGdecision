# 输入教程：raw counts 与本地 HVG

## Raw counts

默认使用：

```python
counts_layer="auto"
```

检查：

```python
check = hd.find_raw_counts(adata)
print(check.valid)
print(check.location)
print(check.audit)
```

只有 auto 失败时才建议人工指定 layer / raw / X。

## 本地 HVG

支持：

```text
.csv   推荐包含 gene 表头
.tsv   推荐包含 gene 表头
.txt   一行一个 gene，可没有表头
```

直接传文件路径：

```python
result = study.refine_hvg(
    "my_hvg.txt",
    initial_n_hvg=2000,
    return_details=True,
)
```

Cross-domain 中，外部 HVG 表始终代表 Query base HVG panel。
