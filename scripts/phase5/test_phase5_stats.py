#!/usr/bin/env python3
"""Focused regression tests for the Phase 5 BH/promoter statistics."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from aggregate_audit import (promoter_q_values, site_q_values, window_count)


class Phase5StatsTests(unittest.TestCase):
    def test_complete_testing_universe(self):
        # 10 bp motif, 37,415 promoters, both strands, 1,000 bp.
        self.assertEqual(window_count(10, 37415, 1000), 37415 * 2 * 991)

    def test_site_q_uses_total_universe_not_hit_count(self):
        q = site_q_values([0.000001, 0.00001], 10, 100, 1000)
        total = 100 * 2 * 991
        self.assertAlmostEqual(q[0], total * 0.000001 / 1)
        self.assertAlmostEqual(q[1], min(1.0, total * 0.00001 / 2))

    def test_truncated_ranking_is_not_hit_only_bh(self):
        # If only two hits are present, q must still use N=100; a hit-only BH
        # would give 0.002/0.02 instead.
        q = site_q_values([0.001, 0.01], 100, 1, 101)  # windows=2*2=4? use N=4
        # This fixture is intentionally tiny; the important assertion is that
        # q is based on total windows, not n=2.
        self.assertGreater(q[0], 0.002)
        self.assertGreater(q[1], 0.02 / 2)

    def test_tied_p_values_share_adjusted_value(self):
        q = site_q_values([0.001, 0.001, 0.01], 1, 100, 1000)
        self.assertAlmostEqual(q[0], q[1])
        self.assertGreaterEqual(q[2], q[1])

    def test_adjusted_values_are_monotone(self):
        p = [0.0001, 0.0001, 0.001, 0.01, 0.02]
        q = site_q_values(p, 10, 1000, 1000)
        order = sorted(range(len(p)), key=lambda i: p[i])
        vals = [q[i] for i in order]
        for left, right in zip(vals, vals[1:]):
            self.assertLessEqual(left, right + 1e-12)

    def test_promoter_window_count_and_q(self):
        # One promoter, one site: Bonferroni over the actual windows, then BH
        # over the single eligible promoter.
        q = promoter_q_values([0], [0.001], 10, 1, 100)
        windows = 2 * 91
        self.assertAlmostEqual(q[0], 0.001 * windows)
        self.assertEqual(window_count(10, 1, 100), windows)

    def test_missing_or_untestable_motif_has_no_q(self):
        self.assertEqual(site_q_values([], 10, 100, 1000), [])
        self.assertEqual(promoter_q_values([], [], 10, 100, 1000), {})
        # A motif longer than the promoter still has a non-zero window count.
        self.assertEqual(window_count(2000, 10, 1000), 10)


if __name__ == "__main__":
    unittest.main()
