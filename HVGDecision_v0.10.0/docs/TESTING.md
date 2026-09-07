# Release testing

The test suite exercises the installed package on small synthetic count matrices.
It checks design routing, disconnected/insufficient support, counts and external
cell-ID alignment, reference-only selection isolation, both risk engines,
explicit gene protection, zero-removal behavior, default Seurat v3 selection,
AnnData export/reload, downstream normalization/PCA/neighbors and the command line.

Regression fixtures were generated independently from the original Generic
cell-level V1 and scenario-specific donor-replicate notebook functions. Tests
compare all original gene evidence columns using the original default risk
parameters (100 permutations and 20 bootstraps for cell-level). The source
provenance and packaged hashes are in ENGINE_PROVENANCE.json. The pandas 3
writable-array adjustment is applied to the reference donor function solely
to execute the same calculation.

Additional tests use small permutation/bootstrap counts for interface checks.
These settings do not replace the defaults. The synthetic positive-removal
test checks that removal and explicit protection execute; it does not measure
the method's biological accuracy.

Local dependency versions and completed test commands are recorded in
TEST_REPORT.json. The wheel is installed independently of the editable source
checkout and the suite is rerun. The source archive and wheel contain the same
Python modules. CI is configured for Python 3.10, 3.11 and 3.12 on Ubuntu;
configured CI platforms are not described as locally tested platforms.

These are software and numerical-regression tests. The release has not rerun
the full PBMC/Human Lung datasets or trained all five integration methods.
It does not certify historical deletion counts or downstream improvements.
Use archived counts, metadata, gene panels and seeds for that validation.
