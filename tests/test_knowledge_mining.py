import importlib.util
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "06_knowledge_mining.py"
)
SPEC = importlib.util.spec_from_file_location("knowledge_mining", MODULE_PATH)
knowledge_mining = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(knowledge_mining)


class OraTests(unittest.TestCase):
    def test_bh_includes_eligible_zero_overlap_gene_sets(self):
        universe = {f"G{i}" for i in range(1, 21)}
        top_genes = {"G1", "G2"}
        gene_sets = {
            "two_hits": {f"G{i}" for i in range(1, 11)},
            "one_hit": {"G1", *{f"G{i}" for i in range(3, 12)}},
            "zero_hits": {f"G{i}" for i in range(3, 13)},
        }

        result = knowledge_mining.ora(top_genes, universe, gene_sets)
        by_term = result.set_index("term")

        p_two_hits = 45 / 190
        p_one_hit = 1 - p_two_hits
        self.assertAlmostEqual(by_term.loc["two_hits", "p"], p_two_hits)
        self.assertAlmostEqual(by_term.loc["one_hit", "p"], p_one_hit)
        self.assertEqual(by_term.loc["zero_hits", "overlap"], 0)
        self.assertEqual(by_term.loc["zero_hits", "genes"], "")
        self.assertEqual(by_term.loc["zero_hits", "p"], 1.0)

        # BH is calculated over all three eligible tests, including zero_hits.
        self.assertAlmostEqual(by_term.loc["two_hits", "fdr"], 3 * p_two_hits)
        self.assertEqual(by_term.loc["one_hit", "fdr"], 1.0)
        self.assertEqual(by_term.loc["zero_hits", "fdr"], 1.0)

    def test_returns_empty_result_when_no_gene_sets_are_eligible(self):
        universe = {f"G{i}" for i in range(1, 21)}
        result = knowledge_mining.ora(
            {"G1", "G2"},
            universe,
            {"too_small": {"G1", "G2"}},
        )

        self.assertTrue(result.empty)

    def test_zero_overlap_only_family_has_unit_p_values_and_fdr(self):
        universe = {f"G{i}" for i in range(1, 21)}
        result = knowledge_mining.ora(
            {"G1", "G2"},
            universe,
            {"zero_hits": {f"G{i}" for i in range(3, 13)}},
        )

        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "p"], 1.0)
        self.assertEqual(result.loc[0, "fdr"], 1.0)


if __name__ == "__main__":
    unittest.main()
