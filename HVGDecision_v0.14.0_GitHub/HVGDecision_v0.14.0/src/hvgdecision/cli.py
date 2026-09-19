"""Command-line interface for HVGDecision v0.14."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd

from . import __version__
from .design import audit_design
from .fusion import UnifiedFusionConfig
from .refinement_workflow import refine


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        prog="hvgdecision",
        description=(
            "Integration-oriented HVG refinement with dual-evidence signed-margin "
            "fusion and dataset-level design audits"
        ),
    )
    command.add_argument("--version", action="version", version=__version__)
    command.add_argument("command", choices=["audit", "refine"])
    command.add_argument("--input", required=True, help="QC-filtered .h5ad")
    command.add_argument("--batch-key", required=True)
    command.add_argument("--label-key", required=True)
    command.add_argument("--reference", nargs="+")
    command.add_argument("--query", nargs="+", help="Explicit Query donors; disjoint from Reference")
    command.add_argument("--hvg-span", default="auto", help="auto or fixed Reference span, e.g. 0.5")
    command.add_argument("--query-hvg-span", help="auto or fixed Query span, e.g. 0.3")
    command.add_argument("--query-hvg-table", help="Ordered Query CSV with gene column")
    command.add_argument("--counts", default="auto", help="auto, X, raw, layer name or counts.csv/tsv")
    command.add_argument("--n-hvg", type=int, default=2000)
    command.add_argument("--hvg-table", help="Ordered Reference/pooled CSV with gene column")
    command.add_argument("--protected-table", help="CSV with gene column")
    command.add_argument("--tau-f", type=float, default=0.10, help="Fusion temperature; default 0.10")
    command.add_argument("--output", required=True, help="NEW/empty output directory")
    command.add_argument("--seed", type=int, default=20260829)
    return command


def main() -> int:
    args = parser().parse_args()
    output = Path(args.output).expanduser()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("Choose a new/empty output directory; existing results preserved")

    adata = ad.read_h5ad(args.input)
    common = dict(batch_key=args.batch_key, label_key=args.label_key, reference=args.reference)

    if args.command == "audit":
        result = audit_design(adata, **common)
        output.mkdir(parents=True, exist_ok=True)
        result.audit.to_csv(output / "design_audit.csv", index=False)
        result.coverage.to_csv(output / "donor_celltype_coverage.csv", index=False)
        print(json.dumps(result.summary, indent=2))
    else:
        panel = pd.read_csv(args.hvg_table)["gene"].astype(str).tolist() if args.hvg_table else None
        protected = (
            pd.read_csv(args.protected_table)["gene"].astype(str).tolist()
            if args.protected_table
            else None
        )
        query_panel = (
            pd.read_csv(args.query_hvg_table)["gene"].astype(str).tolist()
            if args.query_hvg_table
            else None
        )
        span = args.hvg_span if args.hvg_span == "auto" else float(args.hvg_span)
        query_span = (
            args.query_hvg_span
            if args.query_hvg_span in (None, "auto")
            else float(args.query_hvg_span)
        )
        result = refine(
            adata,
            **common,
            counts=args.counts,
            n_hvg=args.n_hvg,
            hvg_genes=panel,
            protected_genes=protected,
            seed=args.seed,
            output_dir=output,
            return_details=True,
            query=args.query,
            query_hvg_genes=query_panel,
            hvg_span=span,
            query_hvg_span=query_span,
            fusion_config=UnifiedFusionConfig(tau=args.tau_f),
        )
        print(result)
        print("Reference/pooled discovery-risk genes:", result.reference_risk_genes)
        print("Final removed genes:", result.removed_genes)

    print("Saved:", output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
