# Release verification: 0.13.0

Verified locally on 2026-09-17. This report distinguishes implementation checks from historical benchmark reproduction.

- Source tests: 12 passed.
- Installed-wheel tests: the same 12 tests passed after installation into an isolated target directory with existing dependencies.
- Both dense and sparse counts; raw-count AnnData roundtrip; no mutation of input; refusal to overwrite results.
- Explicit Reference/Query isolation: changing Query labels and counts with frozen panels does not change Reference evidence or A.
- Both evidence engines, endpoint equivalence at weights 0 and 1, explicit protection, zero-removal handling, signed boundary rules.
- External CSV count alignment; explicit invalid counts rejected; numerical Seurat fits compared with independent Scanpy calls.
- A reproducibility test uses the default 100 permutations and 20 bootstraps. Most remaining tests use smaller synthetic runs.
- All 17 vendored function/class definitions match the supplied continuous-v0.3 prototype by AST comparison, except an explicitly ignored copy=True array-ownership compatibility change in percentile_rank.
- Human Lung metadata-only audit: 60,993 cells, 3 donors, 46 cell types; weighted A=1.7034374235369627. No expression-based risk discovery was rerun on this full dataset during release verification.

The prototype SHA256 is 082146b7cb476eb8979796c88e29f4771b3096fde3e4e10a2947f95e93753ecb.

## Local environment

Python 3.12.14; NumPy 2.5.2; pandas 3.0.5; SciPy 1.18.1; AnnData 0.13.3.post0; Scanpy 1.12.4; scikit-misc 0.5.2; scikit-learn 1.9.0; pytest 9.1.1.

These are release-test versions, not the historical server benchmark environment. The wheel installation used --no-deps with these installed dependencies, not a fresh dependency-resolution test across all supported Python versions.

## Not established by these checks

The full PBMC and Lung gene-discovery plus five-integration-method benchmark was not rerun. The two server examples were syntax-checked but not executed against server data. Exact historical deletion sets and downstream metrics must be validated from frozen inputs and run manifests. The package does not enforce three or ten deletions.
