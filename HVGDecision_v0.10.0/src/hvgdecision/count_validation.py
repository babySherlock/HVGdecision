"""Exhaustive numeric checks; integer-like values cannot prove raw provenance."""
import numpy as np
from scipy import sparse


def audit_counts(matrix, name, tolerance=1e-6):
    if not hasattr(matrix, 'shape') or len(matrix.shape) != 2:
        raise TypeError('Counts must be a two-dimensional matrix')
    shape = matrix.shape
    is_sparse = sparse.issparse(matrix)
    total = finite_n = integer_n = 0
    minimum, maximum = np.inf, -np.inf
    row_step = max(1, 1_000_000 // max(1, shape[1]))
    if is_sparse:
        chunks = (matrix.data[i:i+1_000_000] for i in range(0, len(matrix.data), 1_000_000))
    else:
        chunks = (np.asarray(matrix[i:i+row_step, :]).ravel() for i in range(0, shape[0], row_step))
    error = ''
    for values in chunks:
        if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
            error = f'{name}: counts must be real numeric values'
            break
        total += values.size
        finite = np.isfinite(values)
        finite_n += int(finite.sum())
        v = values[finite]
        if v.size:
            minimum = min(minimum, float(v.min()))
            maximum = max(maximum, float(v.max()))
            integer_n += int(np.isclose(v, np.round(v), atol=tolerance, rtol=0).sum())
    if not error:
        if finite_n != total:
            error = f'{name}: counts contain non-finite values'
        elif total == 0 or maximum <= 0:
            error = f'{name}: count matrix is empty or contains no positive signal'
        elif minimum < 0:
            error = f'{name}: count matrix contains negative values'
        elif integer_n != total:
            error = f'{name}: counts are not integer-like; normalized expression is not valid raw counts'
    return dict(source=name, shape=f'{shape[0]} x {shape[1]}', n_cells=int(shape[0]),
                n_genes=int(shape[1]), storage='sparse' if is_sparse else 'dense',
                dtype=str(getattr(matrix, 'dtype', 'unknown')), sampled_values=total,
                checked_values=total, validation_scope='all_stored_values',
                minimum=float(minimum) if np.isfinite(minimum) else np.nan,
                maximum=float(maximum) if np.isfinite(maximum) else np.nan,
                integer_like_fraction=integer_n / total if total else np.nan,
                valid=not error, error=error)
