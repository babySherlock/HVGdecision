# Examples

建议按顺序：

1. `00_raw_counts_auto.py`：查看工具自动选择的 raw-count source。
2. `01_within_domain.py`：Within-domain 完整最小示例。
3. `02_cross_domain_two_h5ad.py`：从两个独立 h5ad 开始的 Cross-domain 示例。
4. `03_local_hvg_table.py`：导入本地 CSV/TSV/TXT HVG 表，不重新筛 HVG。

所有正式分析必须显式选择：

```text
mode="within_domain"
```

或：

```text
mode="cross_domain"
```
