#!/usr/bin/env python3
"""Synthetic checks for Phase 2.6 ablation helpers."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ablation import global_species_folds, predict_one, top_donors  # noqa: E402
from audit import homology_cluster_map, motif_match_at  # noqa: E402


class Phase26Tests(unittest.TestCase):
    def setUp(self):
        self.donors = ["D1", "D2", "D3", "D4"]
        self.species = {"D1": "sp1", "D2": "sp1", "D3": "sp2", "D4": "sp2"}
        self.family = {donor: "A" for donor in self.donors}
        self.motifs = {"D1": ["M1"], "D2": ["M2"], "D3": ["M3"], "D4": ["M4"]}
        self.groups = {"M1": "G1", "M2": "G1", "M3": "G2", "M4": "G2"}
        self.alignments = {
            "D1": [
                {"subject": "D2", "bitscore": 100.0, "identity": 90.0, "coverage": 0.9, "evalue": 1e-10},
                {"subject": "D3", "bitscore": 50.0, "identity": 60.0, "coverage": 0.8, "evalue": 1e-5},
            ]
        }
        self.reverse = {"D2": [("D1", 100.0), ("D3", 90.0)], "D3": [("D1", 100.0)]}

    def test_species_folds_are_disjoint(self):
        folds = global_species_folds(self.donors, self.species, k=2)
        self.assertEqual(len(folds), 2)
        seen = [donor for fold in folds for donor in fold["donors"]]
        self.assertEqual(sorted(seen), sorted(self.donors))
        self.assertEqual(len(seen), len(set(seen)))

    def test_reciprocal_top_donor_requires_reverse_best(self):
        donors, _, reason = top_donors(
            "D1",
            set(self.donors) - {"D1"},
            self.alignments,
            self.family,
            (30.0, 0.5),
            True,
            self.reverse,
        )
        self.assertEqual(donors, ["D2"])
        self.assertEqual(reason, "ok")
        reverse = {"D2": [("D3", 100.0)]}
        donors, _, reason = top_donors(
            "D1",
            set(self.donors) - {"D1"},
            self.alignments,
            self.family,
            (30.0, 0.5),
            True,
            reverse,
        )
        self.assertEqual(donors, [])
        self.assertEqual(reason, "not_reciprocal")

    def test_nonreciprocal_can_use_multiple_queries(self):
        donors, _, reason = top_donors(
            "D1",
            set(self.donors) - {"D1"},
            self.alignments,
            self.family,
            (30.0, 0.5),
            False,
            {},
        )
        self.assertEqual(donors, ["D2"])
        self.assertEqual(reason, "ok")

    def test_group_aware_prediction_retains_ambiguity(self):
        alignments = {
            "D1": [
                {"subject": "D2", "bitscore": 100.0, "identity": 90.0, "coverage": 0.9, "evalue": 1e-10},
                {"subject": "D3", "bitscore": 90.0, "identity": 80.0, "coverage": 0.9, "evalue": 1e-10},
            ]
        }
        row = predict_one(
            "E",
            "D1",
            {"D2", "D3"},
            alignments,
            {},
            self.family,
            self.motifs,
            self.groups,
            {},
            (30.0, 0.5),
            False,
        )
        self.assertEqual(row["status"], "ambiguous")
        groups = {"M2": "G1", "M3": "G1"}
        row = predict_one(
            "E",
            "D1",
            {"D2", "D3"},
            alignments,
            {},
            self.family,
            self.motifs,
            groups,
            {},
            (30.0, 0.5),
            False,
        )
        self.assertEqual(row["status"], "called")
        self.assertEqual(row["predicted_motifs"], "M2;M3")

    def test_homology_cluster_map_and_overlap_match(self):
        mapping = homology_cluster_map(
            ["D1", "D2"],
            {"D1": [{"subject": "D2", "identity": 90.0, "coverage": 0.9}]},
        )
        self.assertEqual(mapping["D1"], mapping["D2"])
        similar = {("M1", "M2"): (0.95, 0.6)}
        self.assertTrue(motif_match_at(["M1"], ["M2"], similar, 0.80, 0.50))
        self.assertFalse(motif_match_at(["M1"], ["M2"], similar, 0.97, 0.50))
        self.assertFalse(motif_match_at(["M1"], ["M2"], similar, 0.80, 0.70))


if __name__ == "__main__":
    unittest.main()
