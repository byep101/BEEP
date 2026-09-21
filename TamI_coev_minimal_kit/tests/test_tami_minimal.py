import importlib.util
import math
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "tami_minimal.py"
SPEC = importlib.util.spec_from_file_location("tami_minimal", MODULE_PATH)
assert SPEC and SPEC.loader
tami = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tami)


class WeightTests(unittest.TestCase):
    def test_identical_sequences_share_weight(self):
        weights, neighbors, minimum_overlap = tami.compute_weights(
            ["AAAA", "AAAA", "CCCC"], identity_threshold=0.8, overlap_fraction=0.7
        )
        self.assertEqual(minimum_overlap, 3)
        self.assertEqual(neighbors, [2, 2, 1])
        self.assertEqual(weights, [0.5, 0.5, 1.0])
        self.assertEqual(sum(weights), 2.0)

    def test_low_overlap_sequences_are_not_neighbors(self):
        weights, neighbors, _ = tami.compute_weights(
            ["AAA-", "AAA-"], identity_threshold=0.8, overlap_fraction=1.0
        )
        self.assertEqual(neighbors, [1, 1])
        self.assertEqual(weights, [1.0, 1.0])


class MutualInformationTests(unittest.TestCase):
    def test_independent_columns_have_zero_mi(self):
        score, observable = tami.weighted_mi(
            list("AACC"), list("ACAC"), [1.0, 1.0, 1.0, 1.0]
        )
        self.assertAlmostEqual(score, 0.0, places=12)
        self.assertEqual(observable, 4.0)

    def test_perfectly_correlated_binary_columns(self):
        score, observable = tami.weighted_mi(
            list("AACC"), list("AACC"), [1.0, 1.0, 1.0, 1.0]
        )
        self.assertAlmostEqual(score, math.log(2.0), places=12)
        self.assertEqual(observable, 4.0)

    def test_missing_values_are_excluded_pairwise(self):
        score, observable = tami.weighted_mi(
            list("AA-C"), list("AAAC"), [1.0, 1.0, 5.0, 1.0]
        )
        self.assertTrue(math.isfinite(score))
        self.assertEqual(observable, 3.0)


if __name__ == "__main__":
    unittest.main()
