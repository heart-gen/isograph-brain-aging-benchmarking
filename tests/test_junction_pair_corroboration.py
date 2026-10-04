"""Guard against plausible scientific bookkeeping errors in Figure 6c."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.real_data.junction_pair_corroboration import (
    aggregate_counts, canonical_pair, donor_grid, junction_match, nominated_pairs,
    rank_concordance, COVARIATES,
)


class JunctionCorroborationTests(unittest.TestCase):
    def test_nomination_is_an_edge_not_a_clique(self):
        pairs = pd.DataFrame({"gene_id": ["g.1", "g.1", "other.1"],
                              "transcript_id_1": ["a.1", "b.1", "a.1"],
                              "transcript_id_2": ["b.1", "c.1", "c.1"]})
        event = {"gene": "g", "junction_transcripts": "a.1", "switch_pair": "a.1 | b.1 | c.1"}
        got = nominated_pairs(event, pairs)
        self.assertEqual(list(zip(got.transcript_id_1, got.transcript_id_2)), [("a.1", "b.1")])
        self.assertNotEqual(canonical_pair("a.1", "b.1"), canonical_pair("a.2", "b.1"))
        self.assertEqual(canonical_pair("a.1", "b.1"), canonical_pair("b.1", "a.1"))

    def test_coordinate_conversion_rejects_cross_junction_endpoint_match(self):
        j = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [101, 151],
                          "end": [149, 199], "specific": [True, True]})
        self.assertEqual(junction_match("chr1:100-200(+)", {"strand": "+"}, j), (False, False))
        self.assertEqual(junction_match("chr1:100-150(+)", {"strand": "+"}, j), (True, True))
        self.assertEqual(junction_match("chr1:100-150(-)", {"strand": "+"}, j), (False, False))

    def test_streaming_pools_hap0_and_excludes_unselected_donors(self):
        rows = pd.DataFrame({"sample_id": ["s", "s", "s", "case", "s"],
                             "donor_id": ["d", "d", "d", "case_d", "d"],
                             "pair_id": ["p", "p", "p", "p", "other"],
                             "isoform": [1, 1, 2, 1, 1], "hap": [0, 1, 2, 0, 0],
                             "n_frag": [10, 2, 3, 1000, 2000]})
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "counts.parquet"
            rows.to_parquet(p)
            a = aggregate_counts(p, {"s"}, {"p"}, batch_size=1)
            b = aggregate_counts(p, {"s"}, {"p"}, batch_size=100)
        pd.testing.assert_frame_equal(a, b)
        self.assertEqual(a.n_frag.tolist(), [12, 3])

    def test_zero_counts_only_for_completed_samples(self):
        pairs = pd.DataFrame({"count_pair_id": ["p"], "pair_key": ["a|b"],
                              "t1": ["a"], "t2": ["b"], "gene": ["g"], "gene_name": ["G"]})
        completed = pd.DataFrame({"sample_id": ["s1", "s2"], "donor_id": ["d1", "d2"]})
        counts = pd.DataFrame({"sample_id": ["s1"], "donor_id": ["d1"],
                               "pair_id": ["p"], "isoform": [1], "n_frag": [10]})
        d = donor_grid(counts, completed, pairs)
        self.assertEqual(len(d), 2)
        self.assertEqual(d.total.tolist(), [10, 0])
        self.assertEqual(d.n2.tolist(), [0, 0])
        self.assertTrue(np.isnan(d.junction_fraction_1.iloc[1]))

    def test_concordance_missing_covariates_does_not_silently_change_model(self):
        f = pd.DataFrame({"x": range(40), "y": range(40)})
        r = rank_concordance(f, "x", "y")
        self.assertTrue(r["concordance_status"].startswith("missing_covariates:"))
        self.assertTrue(np.isnan(r["rho_adjusted"]))

    def test_covariate_only_agreement_is_not_residual_evidence(self):
        f = pd.DataFrame({c: np.zeros(60) for c in COVARIATES})
        f["Age"] = np.arange(60.)
        f["x"] = f.Age
        f["y"] = 2 * f.Age
        r = rank_concordance(f, "x", "y")
        self.assertEqual(r["concordance_status"], "constant_residual")
        self.assertTrue(np.isnan(r["rho_adjusted"]))


if __name__ == "__main__":
    unittest.main()
