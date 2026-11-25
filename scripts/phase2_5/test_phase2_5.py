#!/usr/bin/env python3
"""Synthetic checks for Phase 2.5 PWM and calibration helpers."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_pairwise import auc  # noqa: E402
from build_specificity_groups import cluster_motifs  # noqa: E402
from calibrate import loso_predictions, motif_match  # noqa: E402
from common import (  # noqa: E402
    information_content,
    matrix_hash,
    motif_similarity,
    reverse_complement_matrix,
)


class Phase25Tests(unittest.TestCase):
    def setUp(self):
        self.matrix = [
            [0.91, 0.03, 0.03, 0.03],
            [0.03, 0.91, 0.03, 0.03],
            [0.03, 0.03, 0.91, 0.03],
            [0.03, 0.03, 0.03, 0.91],
            [0.25, 0.25, 0.25, 0.25],
        ]

    def test_identity_and_reverse_complement_similarity(self):
        self.assertAlmostEqual(motif_similarity(self.matrix, self.matrix)["correlation"], 1.0, places=6)
        reverse = reverse_complement_matrix(self.matrix)
        self.assertAlmostEqual(motif_similarity(self.matrix, reverse)["correlation"], 1.0, places=6)
        self.assertIn(motif_similarity(self.matrix, reverse)["orientation"], {"forward", "reverse"})

    def test_exact_hash_is_stable_and_orientation_specific(self):
        self.assertEqual(matrix_hash(self.matrix), matrix_hash([list(row) for row in self.matrix]))
        self.assertNotEqual(matrix_hash(self.matrix), matrix_hash(reverse_complement_matrix(self.matrix)))

    def test_information_content(self):
        self.assertGreater(information_content(self.matrix), 1.0)
        self.assertAlmostEqual(information_content([[0.25] * 4]), 0.0, places=6)

    def test_direct_motif_match_requires_overlap(self):
        similar = {(("M1", "M2")): (-1.0, 0.0)}
        self.assertTrue(motif_match(["M1"], ["M1"], similar))
        self.assertFalse(motif_match(["M1"], ["M2"], similar))
        similar[(("M1", "M2"))] = (0.95, 0.60)
        self.assertTrue(motif_match(["M1"], ["M2"], similar))
        similar[(("M1", "M2"))] = (0.95, 0.20)
        self.assertFalse(motif_match(["M1"], ["M2"], similar))


    def test_loso_excludes_held_out_donors(self):
        donors = ["D1", "D2", "D3", "D4"]
        species = {"D1": "sp1", "D2": "sp1", "D3": "sp2", "D4": "sp2"}
        motifs = {donor: ["M1"] for donor in donors}
        alignments = {}
        reverse = {}
        for query in donors:
            alignments[query] = [
                {
                    "subject": subject,
                    "bitscore": 100.0 if subject == query else 50.0,
                    "identity": 90.0,
                    "coverage": 0.9,
                    "evalue": 1e-10,
                }
                for subject in donors
            ]
            for subject in donors:
                reverse.setdefault(subject, []).append((query, 100.0 if subject == query else 50.0))
        predictions = loso_predictions(
            donors,
            species,
            alignments,
            reverse,
            motifs,
            {},
            (30.0, 0.5),
        )
        self.assertIsNotNone(predictions)
        self.assertTrue(all(row["top_donor"] != row["query"] for row in predictions if row["covered"]))

    def test_complete_linkage_keeps_exact_groups(self):
        similarities = {
            tuple(sorted(("M1", "M2"))): 0.99,
            tuple(sorted(("M2", "M3"))): 0.99,
            tuple(sorted(("M1", "M3"))): 0.50,
        }
        clusters = cluster_motifs(["M1", "M2", "M3"], similarities, 0.90)
        # Complete linkage must not chain M1-M2-M3 through a weak M1-M3 edge.
        self.assertEqual(len(clusters), 2)

    def test_auc_orders_scores(self):
        rows = [
            {"score": 0.9, "label": True},
            {"score": 0.8, "label": True},
            {"score": 0.2, "label": False},
            {"score": 0.1, "label": False},
        ]
        self.assertAlmostEqual(auc(rows, "score"), 1.0)
        for row in rows:
            row["score"] = 1.0 - row["score"]
        self.assertAlmostEqual(auc(rows, "score"), 0.0)


if __name__ == "__main__":
    unittest.main()
