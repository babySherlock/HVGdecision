# 服务器安装

先激活你运行 notebook 的环境，例如已经使用的 SCVI：

```bash
conda activate SCVI
python -m pip install /path/to/hvgdecision-0.10.0-py3-none-any.whl
python -c "import hvgdecision as hd; print(hd.__version__); print(hd.__file__)"
```

把 `/path/to/` 改成 wheel 实际上传位置。源码安装则进入解压后的文件夹运行
`python -m pip install .`。不要仅把 wheel 放到服务器而不安装。
安装后重启 notebook kernel。Python 环境与 kernel 必须一致。

如需单独新建环境：

```bash
conda create -n hvgdecision python=3.12 pip -y
conda activate hvgdecision
python -m pip install /path/to/hvgdecision-0.10.0-py3-none-any.whl
```

不需要先安装五种整合方法或 R。用法从 `README_CN.md` 或
`examples/01_current_workflow.ipynb` 开始。
