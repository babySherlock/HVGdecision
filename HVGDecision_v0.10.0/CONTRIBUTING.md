# Contributing

Please open an issue before changing statistical defaults. Pull requests should
include a test, an explanation of leakage implications, and an update to the
method documentation when behavior changes.

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
```

