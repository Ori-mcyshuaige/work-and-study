#!/usr/bin/env python3
###############################################################
## Copyright: GeneGenieDx Corp 2022
## Author: whgu
# Date of creation: 04/25/2022
# Date of revision: 04/26/2022
#
## AIM
## Description: Unit test for functions of calculating feature types' importance and
# selecting the best-performing feature types' combination.
#
###############################################################
import unittest
import os
import sys

cur_folder = os.path.dirname(os.path.realpath(__file__))
last_folder = os.path.dirname(cur_folder)
sys.path.append(last_folder)


from select_feature_types import *
from generate_feature_matrix import *


class TestSelectFeature(unittest.TestCase):
    """Test select_feature_types.py"""

    def setUp(self) -> None:
        feature_matrix, feature_type_maps = generate_mock_feature_matrix(200, 3)
        self.feature_matrix = feature_matrix
        self.feature_type_maps = feature_type_maps

    def test_get_feature_type_combination(self) -> None:
        """
        Test function of combining feature types.
        """
        with self.assertRaises(TypeError):
            feature_types_combinations = get_feature_type_combination("t0")

        feature_types = ["t0", "t1", "t2"]
        feature_types_combinations = get_feature_type_combination(feature_types)
        self.assertEqual(
            feature_types_combinations.sort(),
            [
                ("t0",),
                ("t1",),
                ("t2",),
                ("t0", "t1"),
                ("t0", "t2"),
                ("t1", "t2"),
                ("t0", "t1", "t2"),
            ].sort(),
        )

        # Test if the output is correct when feature_types == [].
        feature_types = []
        feature_types_combinations = get_feature_type_combination(feature_types)
        self.assertEqual(
            feature_types_combinations, [],
        )

        # Test if the output is correct when len(feature_types) == 1.
        feature_types = ["t0"]
        feature_types_combinations = get_feature_type_combination(feature_types)
        self.assertEqual(
            feature_types_combinations, [("t0",)],
        )

        # Test if the output is correct when feature_types has duplicates.
        feature_types = ["t0", "t1", "t1"]
        feature_types_combinations = get_feature_type_combination(feature_types)
        self.assertEqual(
            feature_types_combinations.sort(), [("t0",), ("t1",), ("t0", "t1")].sort(),
        )

    def test_feature_type_importance_input(self) -> None:
        """
        Test if the input of calculate_feature_type_importance() are valid.
        """
        calculate_feature_type_importance(self.feature_matrix, self.feature_type_maps)
        calculate_feature_type_importance(
            self.feature_matrix, self.feature_type_maps, model_type=ModelType.DT
        )
        calculate_feature_type_importance(
            self.feature_matrix, self.feature_type_maps, metrics=Metrics.Fscore
        )
        calculate_feature_type_importance(
            self.feature_matrix, self.feature_type_maps, cv=3
        )
        calculate_feature_type_importance(
            self.feature_matrix, self.feature_type_maps, n_jobs=2
        )

        with self.assertRaises(AssertionError):
            calculate_feature_type_importance(self.feature_matrix, {"f-1": "t0"})
        with self.assertRaises(ValueError):
            calculate_feature_type_importance(self.feature_matrix, {})
        with self.assertRaises(TypeError):
            calculate_feature_type_importance(
                self.feature_matrix, self.feature_type_maps, model_type="SVM"
            )
        with self.assertRaises(TypeError):
            calculate_feature_type_importance(
                self.feature_matrix, self.feature_type_maps, metrics="accuracy"
            )
        with self.assertRaises(ValueError):
            calculate_feature_type_importance(
                self.feature_matrix, self.feature_type_maps, cv=1
            )
        with self.assertRaises(ValueError):
            calculate_feature_type_importance(
                self.feature_matrix, self.feature_type_maps, cv=2.0
            )
        with self.assertRaises(ValueError):
            calculate_feature_type_importance(
                self.feature_matrix, self.feature_type_maps, cv=1000
            )
        with self.assertRaises(ValueError):
            calculate_feature_type_importance(
                self.feature_matrix, self.feature_type_maps, n_jobs=0
            )
        with self.assertRaises(ValueError):
            calculate_feature_type_importance(
                self.feature_matrix, self.feature_type_maps, n_jobs=8.0
            )

    def test_feature_type_importance_output(self) -> None:
        """
        Test if the output of calculate_feature_type_importance() is correct.
        """
        feature_type_importance = calculate_feature_type_importance(
            self.feature_matrix, self.feature_type_maps
        )
        for feature_type, importance in feature_type_importance:
            self.assertLessEqual(importance, 1)

        # Test if the output is correct when len(feature_types) == 1.
        feature_matrix, feature_type_maps = generate_mock_feature_matrix(20, 3, [100])
        feature_type_importance = calculate_feature_type_importance(
            feature_matrix, feature_type_maps
        )
        self.assertEqual(feature_type_importance, [("t0", 1.0)])

        # Test if the output is correct when add an extra feature types.
        feature_matrix, feature_type_maps = generate_mock_feature_matrix(20, 3, [100])
        feature_matrix["f101"] = 0
        feature_type_maps["f101"] = "t1"
        feature_type_importance = calculate_feature_type_importance(
            feature_matrix, feature_type_maps
        )
        feature_type_importance.sort(key=lambda t: t[0])
        self.assertEqual(feature_type_importance[1], ("t1", 0))

    def test_select_feature_types_input(self) -> None:
        """
        Test if the input of select_feature_types() are valid.
        """
        select_best_feature_type_combination(
            self.feature_matrix, self.feature_type_maps
        )
        select_best_feature_type_combination(self.feature_matrix, {})
        select_best_feature_type_combination(
            self.feature_matrix, self.feature_type_maps, model_type=ModelType.DT
        )
        select_best_feature_type_combination(
            self.feature_matrix, self.feature_type_maps, metrics=Metrics.Fscore
        )
        select_best_feature_type_combination(
            self.feature_matrix, self.feature_type_maps, cv=3
        )
        select_best_feature_type_combination(
            self.feature_matrix, self.feature_type_maps, n_jobs=2
        )

        with self.assertRaises(AssertionError):
            select_best_feature_type_combination(self.feature_matrix, {"f-1": "t0"})
        with self.assertRaises(TypeError):
            select_best_feature_type_combination(
                self.feature_matrix, self.feature_type_maps, model_type="SVM"
            )
        with self.assertRaises(TypeError):
            select_best_feature_type_combination(
                self.feature_matrix, self.feature_type_maps, metrics="accuracy"
            )
        with self.assertRaises(ValueError):
            select_best_feature_type_combination(
                self.feature_matrix, self.feature_type_maps, cv=1
            )
        with self.assertRaises(ValueError):
            select_best_feature_type_combination(
                self.feature_matrix, self.feature_type_maps, cv=2.0
            )
        with self.assertRaises(ValueError):
            select_best_feature_type_combination(
                self.feature_matrix, self.feature_type_maps, cv=1000
            )
        with self.assertRaises(ValueError):
            select_best_feature_type_combination(
                self.feature_matrix, self.feature_type_maps, n_jobs=0
            )
        with self.assertRaises(ValueError):
            select_best_feature_type_combination(
                self.feature_matrix, self.feature_type_maps, n_jobs=8.0
            )

    def test_select_best_feature_type_combination_output(self) -> None:
        """
        Test if the output of select_best_feature_type_combination() are correct.
        """
        (
            best_feature_types_combination,
            best_performance,
        ) = select_best_feature_type_combination(
            self.feature_matrix, self.feature_type_maps
        )
        self.assertLessEqual(best_performance, 1)
        self.assertGreaterEqual(best_performance, 0)

        # Test if the output is correct when len(feature_types) == 1.
        feature_matrix, feature_type_maps = generate_mock_feature_matrix(20, 3, [100])
        (
            best_feature_types_combination,
            best_performance,
        ) = select_best_feature_type_combination(feature_matrix, feature_type_maps)
        self.assertEqual(best_feature_types_combination, ("t0",))

        # Test if the output is correct when add an extra feature types.
        feature_matrix, feature_type_maps = generate_mock_feature_matrix(20, 3, [100])
        feature_matrix["f101"] = 0
        feature_type_maps["f101"] = "t1"
        (
            best_feature_types_combination,
            best_performance,
        ) = select_best_feature_type_combination(feature_matrix, feature_type_maps)
        self.assertEqual(best_feature_types_combination, ("t0",))


def main() -> None:
    unittest.main()


if __name__ == "__main__":
    main()
