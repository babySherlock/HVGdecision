import numpy as np

from hvgdecision.statistics import benjamini_hochberg, eta_squared, holm_adjust, robust_z


def test_robust_z_is_centered_for_simple_vector():
    values = robust_z(np.array([1.0, 2.0, 3.0]))
    assert np.isclose(np.median(values), 0.0)


def test_eta_squared_detects_group_signal():
    matrix = np.array([[0.0, 0.9], [0.0, 1.1], [4.0, 1.1], [4.0, 0.9]])
    labels = np.array(["a", "a", "b", "b"])
    eta = eta_squared(matrix, labels)
    assert eta[0] > 0.99
    assert eta[1] < 0.1


def test_multiple_testing_adjustments_are_bounded():
    values = np.array([0.001, 0.01, 0.2, 0.9])
    for adjusted in (benjamini_hochberg(values), holm_adjust(values)):
        assert np.all(adjusted >= 0)
        assert np.all(adjusted <= 1)
