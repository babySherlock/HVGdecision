# HVGDecision v0.14.0 test report

## Build validation performed for this release bundle

The package source was syntax-checked and the wheel was built successfully as:

```text
dist/hvgdecision-0.14.0-py3-none-any.whl
```

The built wheel was installed into an isolated target directory with `--no-deps` and the lightweight public API imported successfully:

```text
hvgdecision.__version__ == 0.14.0
REFINEMENT_MODES == ('dual_evidence',)
```

The pure mathematical/audit test suite was run both from source and against the built wheel:

```text
7 passed
```

The full source test collection in the build environment reported:

```text
7 passed, 1 skipped
```

The skipped file is the end-to-end AnnData/Scanpy workflow suite because the artifact-build environment does not contain AnnData/Scanpy. Those tests are included in the repository and run once the declared package dependencies are installed.

## Core properties tested

- centered log-mean-exp symmetry;
- exact preservation when donor-aware and cell-level margins are equal;
- single-component disagreement penalty;
- invalid fusion temperatures are rejected;
- replication adequacy is reported by the design audit but is not used as a fusion weight;
- public v0.14 fusion tables do not expose obsolete confidence-weighting fields;
- component-specific protection sets the corresponding effective margin to `-inf`;
- built-wheel import and version metadata.

## Full workflow tests included

When AnnData, Scanpy and scikit-misc are installed, `tests/test_workflow.py` additionally checks:

- output roundtrip and no mutation of the input object;
- run-manifest fusion metadata;
- Reference/Query separation;
- Query labels are not used for Reference risk discovery;
- explicit protection and valid zero-removal results;
- deterministic repeated runs with the same seed.

## Run locally

```bash
python -m pip install '.[dev]'
python -m pytest -q
```
