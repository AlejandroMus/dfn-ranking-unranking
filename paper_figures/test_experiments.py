import random
import unittest

import paper_algorithm as paper
import benchmark_experiments as experiments


class ExperimentTests(unittest.TestCase):
    def test_four_orders_round_trip(self):
        n, m = 4, 5
        total = paper.dcru.total_dfns(n, m)
        rng = random.Random(20260803)
        for order in experiments.ORDERS:
            paper.clear_preprocessing_caches()
            paper.preprocess_interval_order(order, n, m)
            for index in [0, total - 1] + [rng.randrange(total) for _ in range(20)]:
                levels = experiments._proposed_unrank(order, n, m, index)
                self.assertEqual(
                    experiments._proposed_rank(order, n, m, levels),
                    index,
                    (order, index, levels),
                )

    def test_log_collection_is_opt_in(self):
        sequence, _, log = paper.unrank_tinc_by_cuts(5, 6, 49, collect_log=False)
        self.assertEqual(log, "")
        index, log = paper.rank_tinc_by_cuts_from_levels(
            list(paper.seq_levels(sequence, 6)), 5, 6, collect_log=False
        )
        self.assertEqual(index, 49)
        self.assertEqual(log, "")
        _, _, verbose_log = paper.unrank_tinc_by_cuts(5, 6, 49, collect_log=True)
        self.assertIn("t-inc UNRANK", verbose_log)

    def test_complete_enumeration_count_and_all_orders(self):
        n, m = 3, 4
        for order in experiments.ORDERS:
            brute = experiments.enumerate_and_sort_dfns(n, m, order)
            self.assertEqual(len(brute), paper.dcru.total_dfns(n, m))
            for index, expected in enumerate(brute):
                self.assertEqual(
                    tuple(experiments._proposed_unrank(order, n, m, index)),
                    expected,
                    (order, index),
                )


    def test_peer_review_application_positions(self):
        n, m = 4, 4
        assessments = {
            "A1": [0, 1, 3, 1, 0],
            "A2": [0, 0, 3, 3, 2],
            "A3": [0, 0, 2, 3, 1],
            "A4": [0, 0, 2, 3, 2],
            "A5": [0, 0, 1, 3, 3],
            "A6": [0, 0, 0, 2, 3],
        }
        expected = {
            "lex1": [117, 148, 169, 172, 193, 204],
            "lex2": [92, 144, 156, 170, 193, 204],
            "xy": [97, 145, 161, 174, 193, 204],
            "t-inc": [150, 112, 181, 176, 163, 204],
        }
        for order, positions in expected.items():
            self.assertEqual(
                [
                    experiments._proposed_rank(order, n, m, levels)
                    for levels in assessments.values()
                ],
                positions,
            )


if __name__ == "__main__":
    unittest.main()
