#!/usr/bin/env python3
"""Small unit checks for Phase 5 coordinate grouping and policy helpers."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from aggregate_audit import dedup_exact, group_rows, overlap_coeff, select_policies


def row(motif, start, end, status="dbd_exploratory", rank=(1.0, 0.9, 80.0, 40.0)):
    return {"motif_id": motif, "mapping_status": status, "query_full_dbd_start": start,
            "query_full_dbd_end": end, "query_dbd_coverage": 0.9, "donor_dbd_coverage": 0.9,
            "dbd_identity": 80.0, "dbd_bitscore": 40.0, "_qs": start, "_qe": end,
            "_qcov": 0.9, "_dcov": 0.9, "_id": 80.0, "_bits": 40.0, "_rank": rank,
            "extraction_family": "bHLH", "domain_mode": "full_architecture", "domain_count": "1"}


class Phase5Tests(unittest.TestCase):
    def test_overlap(self):
        self.assertAlmostEqual(overlap_coeff((10, 50), (20, 60)), 31 / 41)
        self.assertEqual(overlap_coeff((10, 20), (30, 40)), 0.0)

    def test_one_region_groups(self):
        rows = [row("M1", 10, 50), row("M2", 11, 50), row("M3", 12, 49)]
        self.assertEqual(len(group_rows(rows, 0.8)), 1)

    def test_disjoint_regions_split(self):
        rows = [row("M1", 10, 50), row("M2", 200, 240)]
        self.assertEqual(len(group_rows(rows, 0.8)), 2)

    def test_exact_pwm_dedup(self):
        rows = [row("M1", 10, 50), row("M2", 10, 50)]
        hashes = {"M1": "p1", "M2": "p1"}
        self.assertEqual(len(dedup_exact(rows, hashes)), 1)

    def test_best_is_one_per_group(self):
        rows = [row("M1", 10, 50, rank=(1.0, 0.9, 80.0, 40.0)),
                row("M2", 11, 50, rank=(1.0, 0.8, 79.0, 39.0))]
        policies = select_policies(rows, {}, {})
        self.assertEqual(len(policies["all"]), 2)
        self.assertLessEqual(len(policies["best@0.8"]), 1)
        self.assertLessEqual(len(policies["reps@0.8"]), 2)


if __name__ == "__main__":
    unittest.main()
