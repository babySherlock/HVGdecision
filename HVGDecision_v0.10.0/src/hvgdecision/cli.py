"""Current multi-donor CLI. Historical YAML: hvgdecision-legacy."""

from __future__ import annotations

import argparse
import json

from . import __version__
from pathlib import Path
import anndata as ad
import pandas as pd
from .routing import audit_design
from .workflow import refine


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        prog="hvgdecision",
        description="Auditable multi-donor HVG refinement with automatic routing",
    )
    command.add_argument("--version", action="version", version=__version__)
    command.add_argument('command', choices=['audit', 'refine'])
    command.add_argument('--input', required=True, help='QC-filtered .h5ad')
    command.add_argument('--batch-key', required=True)
    command.add_argument('--label-key', required=True)
    command.add_argument('--reference', nargs='+')
    command.add_argument('--counts', default='auto', help='auto, X, raw, layer name or counts.csv/tsv')
    command.add_argument('--n-hvg', type=int, default=2000)
    command.add_argument('--hvg-table', help='Ordered CSV with gene column')
    command.add_argument('--protected-table', help='CSV with gene column')
    command.add_argument('--output', required=True, help='NEW/empty output directory')
    command.add_argument('--seed', type=int, default=20260829)
    return command


def main() -> int:
    args = parser().parse_args()
    output = Path(args.output).expanduser()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError('Choose a new/empty output directory; existing results preserved')
    adata = ad.read_h5ad(args.input)
    common = dict(batch_key=args.batch_key, label_key=args.label_key, reference=args.reference)
    if args.command == 'audit':
        result = audit_design(adata, **common)
        output.mkdir(parents=True, exist_ok=True)
        result.audit.to_csv(output / 'routing_summary.csv', index=False)
        result.coverage.to_csv(output / 'donor_celltype_coverage.csv', index=False)
        print(json.dumps(result.summary, indent=2))
    else:
        panel = pd.read_csv(args.hvg_table)['gene'].tolist() if args.hvg_table else None
        protected = pd.read_csv(args.protected_table)['gene'].tolist() if args.protected_table else None
        result = refine(adata, **common, counts=args.counts, n_hvg=args.n_hvg, hvg_genes=panel,
                        protected_genes=protected, seed=args.seed, output_dir=output, return_details=True)
        print(result)
        print('Removed genes:', result.removed_genes)
    print('Saved:', output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
