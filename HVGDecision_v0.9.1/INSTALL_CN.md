# HVGDecision 0.9.1 安装

普通用户：

```bash
python -m pip install ./hvgdecision-0.9.1-py3-none-any.whl
```

验证：

```bash
python -c "import hvgdecision as hd; print(hd.__version__)"
```

开发者需要修改源码时才使用：

```bash
python -m pip install -e ".[dev]"
```
