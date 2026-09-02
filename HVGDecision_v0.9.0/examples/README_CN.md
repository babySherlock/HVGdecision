# Examples

- `quickstart_in_memory.py`：`within_domain` 最小示例。
- `example_1_python_hvg.py`：`cross_domain`，内部选择 Query HVG，然后运行 Cross-domain Rule V3。
- `example_2_external_hvg.py`：对外部冻结 HVG 列表运行所选 mode。

0.9.0 中 `mode` 不再隐式推断，必须明确为：

```text
within_domain
cross_domain
```
