import unittest
import sys
import os
from pathlib import Path

# Current file's directory
current_dir = Path(__file__).resolve().parent

# Go up one level
parent_dir = current_dir.parent

# Add to sys.path
sys.path.append(str(parent_dir))

# Add path to module
# Add path
from utils import *

class TestUtils(unittest.TestCase):
    def setUp(self):
        self.power_of_two_test_cases = [
            {
                "n": 1,
                "is_power_of_two": True
            },
            {
                "n": 2,
                "is_power_of_two": True
            },
            {
                "n": 3,
                "is_power_of_two": False
            }
        ]
        self.connection_test_cases = [
            {
                "nodes": ((0, 1), (2, 3), (4, 5)),
                "is_connected": False
            },
            {
                "nodes": ((0, 1), (1, 2), (2, 3)),
                "is_connected": True
            }
        ]
        self.cost_test_cases = [
            {
                "Xs": [0b1111, 0b0001, 0b1101],
                "Zs": [(1, 0), (1, 12)],
                "expected_best_Xs_reduced": [1, 2, 12],
                "expected_best_cost": 12
            }
        ]

    def test_is_power_of_two(self):
        for case in self.power_of_two_test_cases:
            with self.subTest(case=case):
                n = case["n"]
                result = is_power_of_two(n)
                self.assertEqual(result, case["is_power_of_two"])
        
    def test_is_connected(self):
        for case in self.connection_test_cases:
            with self.subTest(case=case):
                nodes = case["nodes"]
                result = is_connected(nodes)
                self.assertEqual(result, case["is_connected"])

    def test_find_best_cost(self):
        for case in self.cost_test_cases:
            with self.subTest(case=case):
                best_Xs_reduced, best_cost = find_best_cost(case["Xs"], case["Zs"])
                self.assertEqual(set(best_Xs_reduced), set(case["expected_best_Xs_reduced"]))
                self.assertEqual(best_cost, case["expected_best_cost"])


if __name__ == '__main__':
    unittest.main()