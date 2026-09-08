# How to merge this into the main HVGDecision v0.10.0 source

This companion wheel is intentionally core-preserving. For a true main-package v0.10.1 release:

1. Start from the exact current `HVGDecision_v0.10.0` source tree.
2. Add the Seurat-v3 retry sequence to the internal HVG selector: `0.3 -> 0.5 -> 0.7 -> 1.0`, recreating the fit object for each attempt.
3. Record attempted spans and the successful span in the audit.
4. Add a public `refine_reference_query()` wrapper with the exact semantics implemented in this patch:
   - Reference HVG selection starts at 0.3;
   - call the existing v0.10.0 donor-aware `refine()` on frozen Reference HVGs;
   - Query HVG selection is independent and starts again at 0.3;
   - delete only Reference-risk ∩ Query-HVG.
5. Export the wrapper from `hvgdecision.__init__`.
6. Add the Chinese/English README tutorial.
7. Run the Kang2018 regression test and require:
   - Reference span = 0.5;
   - Query span = 0.3;
   - Reference risk = 19;
   - Query deletion = CSF2, GJB2, AC147651.3;
   - final Query panel = 1997.
