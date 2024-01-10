#!/usr/bin/env python3
###############################################################
# Copyright: GeneGenieDx Corp 2022
# Author: whgu
# Date of creation: 04/03/2022
# Date of revision: 04/19/2022
## AIM
## Description: Unit test for generate_mock_feature_matrix().
#
## usage:
#   Test_generate_matrix.py
###############################################################
import unittest
import os
import sys
import numpy as np

cur_folder = os.path.dirname(os.path.realpath(__file__))
last_folder = os.path.dirname(cur_folder)
sys.path.append(last_folder)

from generate_feature_matrix import *


class TestGenerateMatrix(unittest.TestCase):
    """Test generate_feature_matrix.py"""

    def test_generate_mock_feature_matrix_input(self) -> None:
        """
        Test whether the input of generate_mock_feature_matrix() are valid.
        """
        with self.assertRaises(ValueError):
            generate_mock_feature_matrix(n_samples=-1)

        with self.assertRaises(TypeError):
            generate_mock_feature_matrix(n_samples=20.0)

        with self.assertRaises(AssertionError):
            generate_mock_feature_matrix(n_classes=0)

        with self.assertRaises(AssertionError):
            generate_mock_feature_matrix(n_classes=-1)

        with self.assertRaises(TypeError):
            generate_mock_feature_matrix(n_classes=3.0)

        with self.assertRaises(ValueError):
            generate_mock_feature_matrix(feature_count_per_type=[-10])

        with self.assertRaises(TypeError):
            generate_mock_feature_matrix(feature_count_per_type=[10.0])

        with self.assertRaises(ValueError):
            generate_mock_feature_matrix(cluster_std=-1)

        with self.assertRaises(TypeError):
            generate_mock_feature_matrix(random_state=2.0)

    def test_generate_mock_feature_matrix_output(self) -> None:
        """
        Test if the output of generate_mock_feature_matrix() are correct.
        """
        feature_matrix, feature_to_type = generate_mock_feature_matrix(
            n_samples=10, feature_count_per_type=[1, 2, 3]
        )

        self.assertIsInstance(feature_matrix, pd.DataFrame)
        self.assertEqual(feature_matrix.shape, (10, 1 + 2 + 3 + 1))
        self.assertEqual(feature_matrix.columns.values.tolist()[-1], "label")
        for label in feature_matrix["label"].values.tolist():
            self.assertIn(label, [0, 1, 2])
        self.assertDictEqual(
            feature_to_type,
            {"f0": "t0", "f1": "t1", "f2": "t1", "f3": "t2", "f4": "t2", "f5": "t2"},
        )

        # Test if the output are correct when n_sample == 0.
        feature_matrix, _ = generate_mock_feature_matrix(
            n_samples=0, n_classes=3, feature_count_per_type=[100, 150, 200]
        )
        self.assertEqual(feature_matrix.shape, (0, 100 + 150 + 200 + 1))

        # Test if the output are correct when feature_count_per_type == [].
        feature_matrix, feature_to_type = generate_mock_feature_matrix(
            n_samples=150, n_classes=3, feature_count_per_type=[]
        )
        self.assertEqual(feature_matrix.shape, (150, 1))
        self.assertDictEqual(feature_to_type, {})

        # Test if the output are correct when n_classes == 1.
        feature_matrix, feature_to_type = generate_mock_feature_matrix(n_classes=1)
        for label in feature_matrix["label"].values.tolist():
            self.assertEqual(label, 0)

        # Test whether the generated samples are reproducible.
        feature_matrix_1, _ = generate_mock_feature_matrix(random_state=123)
        feature_matrix_2, _ = generate_mock_feature_matrix(random_state=123)
        assert np.array_equal(feature_matrix_1.values, feature_matrix_2.values)


def main() -> None:
    unittest.main()


if __name__ == "__main__":
    main()
