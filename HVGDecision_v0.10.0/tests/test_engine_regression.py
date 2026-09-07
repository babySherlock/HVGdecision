from pathlib import Path
import pandas as pd
from test_current_workflow import example
from hvgdecision._cell_level import find_harmful_v1_celllevel
from hvgdecision._donor_aware import find_harmful_within_protocol_donor_replicate


def test_frozen_cell_engine_default_parameters():
    a = example()
    table = find_harmful_v1_celllevel(a, a.X, a.var_names, a.var_names.tolist(), 'donor', 'type')['table']
    expected = pd.read_csv(Path(__file__).parent / 'fixtures/cell_level_original.csv', keep_default_na=False)
    table = table.sort_values('gene').reset_index(drop=True)
    pd.testing.assert_frame_equal(table, expected, check_dtype=False, rtol=1e-7, atol=1e-9)


def test_frozen_donor_engine_default_parameters():
    a = example()
    table = find_harmful_within_protocol_donor_replicate(a, a.X, a.var_names, a.var_names.tolist(), 'donor', 'type')['table']
    expected = pd.read_csv(Path(__file__).parent / 'fixtures/donor_aware_original.csv', keep_default_na=False)
    table = table.sort_values('gene').reset_index(drop=True)
    pd.testing.assert_frame_equal(table, expected, check_dtype=False, rtol=1e-7, atol=1e-9)
