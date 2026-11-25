#!/usr/bin/env python3
"""Unit tests for the Phase 3 independent scanner validation."""

import math
import unittest

import synthetic_validation as sv


class SyntheticValidationTests(unittest.TestCase):
    def test_probability_normalization_sums_to_one(self):
        matrix = sv.normalize_probability_matrix(
            [[0.7, 0.1, 0.1, 0.1], [0.25, 0.25, 0.25, 0.25]]
        )
        for column in matrix:
            self.assertAlmostEqual(sum(column), 1.0, places=12)

    def test_null_matches_bruteforce_two_positions(self):
        matrix = sv.normalize_probability_matrix(
            [[0.7, 0.1, 0.1, 0.1], [0.4, 0.3, 0.2, 0.1]]
        )
        background = [0.25] * 4
        odds = sv.log_odds(matrix, background)
        rows, tail = sv.null_distribution(odds, background)
        dist = {}
        for a in range(4):
            for b in range(4):
                score = rows[0][a] + rows[1][b]
                dist[score] = dist.get(score, 0.0) + 0.25 * 0.25
        keys, values = sv.tail_arrays(tail)
        for score in dist:
            expected = sum(p for other, p in dist.items() if other >= score)
            got = sv.tail_pvalue_sorted(keys, values, score)
            self.assertAlmostEqual(got, expected, places=12)

    def test_tail_pvalue_is_monotone(self):
        matrix = sv.normalize_probability_matrix(
            [[0.7, 0.1, 0.1, 0.1], [0.4, 0.3, 0.2, 0.1], [0.1, 0.1, 0.7, 0.1]]
        )
        background = [0.25] * 4
        rows, tail = sv.null_distribution(sv.log_odds(matrix, background), background)
        keys, values = sv.tail_arrays(tail)
        best = sum(max(row) for row in rows)
        worst = sum(min(row) for row in rows)
        self.assertGreaterEqual(
            sv.tail_pvalue_sorted(keys, values, worst),
            sv.tail_pvalue_sorted(keys, values, best),
        )

    def test_binomial_tail_known_value(self):
        # P(X >= 2) for Binomial(4, 0.5) = 11/16
        self.assertAlmostEqual(sv.binomial_tail(4, 2, 0.5), 11 / 16, places=12)
        self.assertEqual(sv.binomial_tail(4, 0, 0.5), 1.0)
        self.assertEqual(sv.binomial_tail(4, 5, 0.5), 0.0)

    def test_gc_frequencies(self):
        frequencies = sv.base_frequencies_gc(0.8)
        self.assertAlmostEqual(sum(frequencies), 1.0, places=12)
        self.assertAlmostEqual(frequencies[0], frequencies[3], places=12)
        self.assertAlmostEqual(frequencies[1], frequencies[2], places=12)
        self.assertAlmostEqual(frequencies[1] + frequencies[2], 0.8, places=12)

    def test_ks_uniform_grid(self):
        values = [(index + 0.5) / 1000 for index in range(1000)]
        d, p = sv.ks_uniform(values)
        self.assertLess(d, 0.002)
        self.assertGreater(p, 0.99)


if __name__ == "__main__":
    unittest.main()
