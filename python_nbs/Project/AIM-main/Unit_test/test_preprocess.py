#!/usr/bin/env python3
###############################################################
## Copyright: GeneGenieDx Corp 2022
# Author: whgu
# Date of creation: 04/04/2022
# Date of revision: 05/24/2022
#
## AIM
## Description: Unit test for functions of preprocessing the features,
# prepare_train_test_data and k-fold split.
#
###############################################################
import unittest
import os
import sys

cur_folder = os.path.dirname(os.path.realpath(__file__))
last_folder = os.path.dirname(cur_folder)
sys.path.append(last_folder)

from preprocess import *

from generate_feature_matrix import *


class TestPreprocess(unittest.TestCase):
    """Test preprocess.py"""

    def setUp(self) -> None:
        """
        Generate feature_matrix and feature_map.
        """
        feature_matrix, feature_map = generate_mock_feature_matrix()
        self.feature_matrix = feature_matrix
        self.feature_map = feature_map

    def test_preprocess_input(self) -> None:
        """
        Test whether the input of proprecess() are valid.
        """
        # Test if the input features of proprecess() is valid.
        preprocess(features=np.array([[1], [2]]))
        preprocess(features=np.array([[1, 2, 3]]))

        with self.assertRaises(TypeError):
            preprocess(features=[[1], [2]])

        with self.assertRaises(ValueError):
            preprocess(features=np.array([]))

        with self.assertRaises(ValueError):
            preprocess(features=np.array([1, 2, 3]))

        # Test if the input impute_missing_value of proprecess() is valid.
        preprocess(features=np.array([[1], [2]]), impute_missing_value=True)
        preprocess(features=np.array([[1], [2]]), impute_missing_value=False)
        preprocess(features=np.array([[1], [2]]), impute_missing_value=1)
        preprocess(features=np.array([[1], [2]]), impute_missing_value=0)

        # Test if the input scale_by_minmax of proprecess() is valid.
        preprocess(features=np.array([[1], [2]]), scale_by_minmax=True)
        preprocess(features=np.array([[1], [2]]), scale_by_minmax=False)
        preprocess(features=np.array([[1], [2]]), scale_by_minmax=1)
        preprocess(features=np.array([[1], [2]]), scale_by_minmax=0)

        # Testif the input impute_strategy of proprecess() is valid.
        preprocess(features=np.array([[1], [2]]), impute_strategy="mean")
        preprocess(features=np.array([[1], [2]]), impute_strategy="median")
        preprocess(features=np.array([[1], [2]]), impute_strategy="most_frequent")

        with self.assertRaises(AssertionError):
            preprocess(features=np.array([[1], [2]]), impute_strategy="common")

        # Test if the input scale_range of proprecess() is valid.
        preprocess(features=np.array([[1], [2]]), scale_range=(0, 1))
        preprocess(features=np.array([[1], [2]]), scale_range=(-1, 1))

        with self.assertRaises(ValueError):
            preprocess(features=np.array([[1], [2]]), scale_range=(1, 0))

        with self.assertRaises(ValueError):
            preprocess(features=np.array([[1], [2]]), scale_range=(0, 0))

    def test_imputation(self) -> None:
        """
        Test function of imputing missing values with different impute strategies in preprocess().
        """
        # Test mean imputation.
        X = np.array([[1, 2, 3], [4, np.nan, 6], [7, 8, 9]])
        preprocessed_features = preprocess(X)

        self.assertIsInstance(preprocessed_features, np.ndarray)

        X_ = np.array([[1, 2, 3], [4, np.mean([2, 8]), 6], [7, 8, 9]], dtype=float)
        scale = 1 / (X_.max(axis=0) - X_.min(axis=0))
        X_scaled = scale * X_ - X_.min(axis=0) * scale
        assert np.array_equal(preprocessed_features, X_scaled)

        # Test median imputation.
        X = np.array([[1, 2, 3], [4, np.nan, 6], [7, 8, 9], [10, 11, 12]])
        preprocessed_features = preprocess(X, impute_strategy="median")

        X_ = np.array([[1, 2, 3], [4, 8, 6], [7, 8, 9], [10, 11, 12]], dtype=float)
        scale = 1 / (X_.max(axis=0) - X_.min(axis=0))
        X_scaled = scale * X_ - X_.min(axis=0) * scale
        assert np.array_equal(preprocessed_features, X_scaled)

        # Test most_frequent imputation.
        X = np.array([[1, 2, 3], [4, np.nan, 6], [7, 8, 9], [10, 2, 12]])
        preprocessed_features = preprocess(X, impute_strategy="most_frequent")

        X_ = np.array([[1, 2, 3], [4, 2, 6], [7, 8, 9], [10, 2, 12]], dtype=float)
        scale = 1 / (X_.max(axis=0) - X_.min(axis=0))
        X_scaled = scale * X_ - X_.min(axis=0) * scale
        assert np.array_equal(preprocessed_features, X_scaled)

    def test_drop_column(self) -> None:
        """
        Test function of dropping columns that contain missing values in preprocess().
        """
        # Test dropping columns containing missing values.
        X = np.array([[1, 2, 3], [4, np.nan, 6], [7, 8, 9]])
        preprocessed_features = preprocess(X, impute_missing_value=False)

        X_ = np.array([[1, 3], [4, 6], [7, 9]], dtype=float)
        scale = 1 / (X_.max(axis=0) - X_.min(axis=0))
        X_scaled = scale * X_ - X_.min(axis=0) * scale
        assert np.array_equal(preprocessed_features, X_scaled)

    def test_MinMaxScaler(self) -> None:
        """
        Test function of transforming features by scaling each feature to a given range in preprocess().
        """
        # Test MinMaxScaler with (0,1).
        X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        preprocessed_features = preprocess(X)

        scale = 1 / (X.max(axis=0) - X.min(axis=0))
        X_scaled = scale * X - X.min(axis=0) * scale
        assert np.array_equal(preprocessed_features, X_scaled)

        # Test MinMaxScaler with (-1,1).
        X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        preprocessed_features = preprocess(X, scale_range=(-1, 1))

        scale = 2 / (X.max(axis=0) - X.min(axis=0))
        min_ = -1 - X.min(axis=0) * scale
        X_scaled = scale * X + min_
        assert np.array_equal(preprocessed_features, X_scaled)

    def test_StandardScaler(self) -> None:
        """
        Test function of converting features into a standard normal distribution in preprocess().
        """
        X = np.array([[0, 0], [0, 0], [1, 1], [1, 1]])
        preprocessed_features = preprocess(X, scale_by_minmax=False)

        target = np.array([[-1, -1], [-1, -1], [1, 1], [1, 1]])
        assert np.array_equal(preprocessed_features, target)

    def test_prepare_train_test_data_input(self) -> None:
        """
        Test whether the input of prepare_train_test_data() are valid.
        """
        # Test if the input data of prepare_train_test_data() is valid.
        prepare_train_test_data(
            self.feature_matrix, features=["f0"], label_column_name="label"
        )
        with self.assertRaises(TypeError):
            prepare_train_test_data(
                np.array([[1, 2, 3]]), features=["f0"], label_column_name="label"
            )
        prepare_train_test_data(
            self.feature_matrix, features=["f0"], label_column_name="f1"
        )

        # Test if the input features of prepare_train_test_data() is valid.
        prepare_train_test_data(
            self.feature_matrix,
            features=["f0", "f15", "f20"],
            label_column_name="label",
        )
        with self.assertRaises(TypeError):
            prepare_train_test_data(
                self.feature_matrix, features="f0", label_column_name="label"
            )
        with self.assertRaises(AssertionError):
            prepare_train_test_data(
                self.feature_matrix, features=["f500"], label_column_name="label"
            )

        # Test if the input label_column_name of prepare_train_test_data() is valid.
        with self.assertRaises(AssertionError):
            prepare_train_test_data(
                self.feature_matrix, features=["f0"], label_column_name="class"
            )

        # Test if the input test_ratio of prepare_train_test_data() is valid.
        prepare_train_test_data(self.feature_matrix, features=["f0"], test_ratio=0)
        prepare_train_test_data(self.feature_matrix, features=["f0"], test_ratio=1)
        with self.assertRaises(ValueError):
            prepare_train_test_data(
                self.feature_matrix, features=["f0"], test_ratio=1.1
            )
        with self.assertRaises(ValueError):
            prepare_train_test_data(
                self.feature_matrix, features=["f0"], test_ratio=-0.1
            )

        # Test if the input shuffle_data of prepare_train_test_data() is valid.
        prepare_train_test_data(self.feature_matrix, features=["f0"], shuffle_data=True)
        prepare_train_test_data(
            self.feature_matrix, features=["f0"], shuffle_data=False
        )

        # Test if the input random_state of prepare_train_test_data() is valid.
        prepare_train_test_data(
            self.feature_matrix, features=["f0"], shuffle_data=True, random_state=None
        )
        prepare_train_test_data(
            self.feature_matrix, features=["f0"], shuffle_data=True, random_state=123
        )
        with self.assertRaises(TypeError):
            prepare_train_test_data(
                self.feature_matrix,
                features=["f0"],
                label_column_name="label",
                random_state=1.0,
            )

    def test_prepare_train_test_data_output(self) -> None:
        """
        Test whether the output of prepare_train_test_data() are correct.
        """
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            self.feature_matrix,
            ["f0", "f2", "f100", "f200"],
            label_column_name="label",
            shuffle_data=False,
        )

        total_X = self.feature_matrix.iloc[:, :-1].values
        total_y = self.feature_matrix.iloc[:, -1].values
        selected_feature_columns = [
            i
            for i, feature_type in enumerate(self.feature_matrix.columns.values)
            if feature_type in ["f0", "f2", "f100", "f200"]
        ]

        assert np.array_equal(
            test_X.values,
            total_X[0 : int(len(self.feature_matrix) * 0.2), selected_feature_columns],
        )
        assert np.array_equal(
            train_X.values,
            total_X[int(len(self.feature_matrix) * 0.2) :, selected_feature_columns],
        )

        assert np.array_equal(test_y, total_y[0 : int(len(self.feature_matrix) * 0.2)])
        assert np.array_equal(train_y, total_y[int(len(self.feature_matrix) * 0.2) :])

        # Test if the output are correct when input features == [].
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            self.feature_matrix, [], label_column_name="label",
        )
        self.assertEqual(train_X.shape[1], 0)
        self.assertEqual(test_X.shape[1], 0)

        # Test if the output are correct when input features == all features.
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            self.feature_matrix,
            list(self.feature_map.keys()),
            label_column_name="label",
        )
        self.assertEqual(train_X.shape[1], self.feature_matrix.shape[1] - 1)
        self.assertEqual(test_X.shape[1], self.feature_matrix.shape[1] - 1)

        # Test if the output are correct when input test_ratio == 0.
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            self.feature_matrix, ["f0", "f2", "f100", "f200"], test_ratio=0
        )
        self.assertEqual(train_X.shape[0], len(self.feature_matrix))
        self.assertEqual(train_y.shape[0], len(self.feature_matrix))
        self.assertEqual(test_X.shape[0], 0)
        self.assertEqual(test_y.shape[0], 0)

        # Test if the output are correct when input test_ratio == 1.
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            self.feature_matrix, ["f0", "f2", "f100", "f200"], test_ratio=1
        )
        self.assertEqual(train_X.shape[0], 0)
        self.assertEqual(train_y.shape[0], 0)
        self.assertEqual(test_X.shape[0], len(self.feature_matrix))
        self.assertEqual(test_y.shape[0], len(self.feature_matrix))

    def test_prepare_train_test_data_shuffle(self) -> None:
        """
        Test whether shuffling data of prepare_train_test_data() works.
        """
        # Test if the output are correct when shuffe_data == True.
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            self.feature_matrix,
            ["f0", "f2", "f100", "f200"],
            label_column_name="label",
            shuffle_data=True,
        )

        total_X = self.feature_matrix.iloc[:, :-1].values
        total_y = self.feature_matrix.iloc[:, -1].values
        selected_feature_columns = [
            i
            for i, feature_type in enumerate(self.feature_matrix.columns.values)
            if feature_type in ["f0", "f2", "f100", "f200"]
        ]

        self.assertEqual(
            train_X.shape,
            total_X[
                int(len(self.feature_matrix) * 0.2) :, selected_feature_columns
            ].shape,
        )
        self.assertEqual(
            test_X.shape,
            total_X[
                0 : int(len(self.feature_matrix) * 0.2), selected_feature_columns
            ].shape,
        )
        self.assertEqual(
            train_y.shape, total_y[int(len(self.feature_matrix) * 0.2) :].shape,
        )
        self.assertEqual(
            test_y.shape, total_y[0 : int(len(self.feature_matrix) * 0.2) :].shape,
        )

        assert not np.array_equal(
            train_X.values,
            total_X[int(len(self.feature_matrix) * 0.2) :, selected_feature_columns],
        )
        assert not np.array_equal(
            test_X.values,
            total_X[0 : int(len(self.feature_matrix) * 0.2), selected_feature_columns],
        )
        assert not np.array_equal(
            train_y, total_y[int(len(self.feature_matrix) * 0.2) :]
        )
        assert not np.array_equal(
            test_y, total_y[0 : int(len(self.feature_matrix) * 0.2)]
        )

        # Test whether the shuffed results are reproducible.
        train_X_1, train_y_1, test_X_1, test_y_1 = prepare_train_test_data(
            self.feature_matrix,
            ["f0", "f2", "f100", "f200"],
            label_column_name="label",
            shuffle_data=True,
            random_state=123,
        )
        train_X_2, train_y_2, test_X_2, test_y_2 = prepare_train_test_data(
            self.feature_matrix,
            ["f0", "f2", "f100", "f200"],
            label_column_name="label",
            shuffle_data=True,
            random_state=123,
        )
        assert np.array_equal(train_X_1, train_X_2)
        assert np.array_equal(train_y_1, train_y_2)
        assert np.array_equal(test_X_1, test_X_2)
        assert np.array_equal(test_y_1, test_y_2)

    def test_k_fold_split_input(self) -> None:
        """
        Test whether the input of k_fold_split() are valid.
        """
        # Test if the input data of k_fold_split() is valid.
        k_fold_split(
            self.feature_matrix,
            features=["f0", "f10", "f100"],
            label_column_name="label",
        )
        with self.assertRaises(TypeError):
            k_fold_split(
                np.array([[1, 2, 3]]),
                features=["f0", "f10", "f100"],
                label_column_name="label",
            )

        # Test if the input features of k_fold_split() is valid.
        k_fold_split(
            self.feature_matrix, features=["f0"], label_column_name="label",
        )
        with self.assertRaises(TypeError):
            k_fold_split(self.feature_matrix, features="f0", label_column_name="label")
        with self.assertRaises(AssertionError):
            k_fold_split(
                self.feature_matrix, features=["f500"], label_column_name="label"
            )

        # Test if the input label_column_name of k_fold_split() is valid.
        k_fold_split(
            self.feature_matrix, features=["f0", "f10", "f100"], label_column_name="f1"
        )
        with self.assertRaises(AssertionError):
            k_fold_split(
                self.feature_matrix, features=["f0"], label_column_name="class"
            )

        # Test if the input n_splits of k_fold_split() is valid.
        k_fold_split(self.feature_matrix, features=["f0"], n_splits=10)
        k_fold_split(self.feature_matrix, features=["f0"], n_splits=2)
        with self.assertRaises(ValueError):
            k_fold_split(self.feature_matrix, features=["f0"], n_splits=1)
        with self.assertRaises(ValueError):
            k_fold_split(self.feature_matrix, features=["f0"], n_splits=10000)
        with self.assertRaises(ValueError):
            k_fold_split(self.feature_matrix, features=["f0"], n_splits=3.0)

        # Test if the input shuffle of k_fold_split() is valid.
        k_fold_split(self.feature_matrix, features=["f0"], shuffle_data=True)
        k_fold_split(self.feature_matrix, features=["f0"], shuffle_data=False)

        # Test if the input random_state of k_fold_split() is valid.
        k_fold_split(
            self.feature_matrix, features=["f0"], shuffle_data=True, random_state=None
        )
        k_fold_split(
            self.feature_matrix, features=["f0"], shuffle_data=True, random_state=123
        )
        with self.assertRaises(TypeError):
            k_fold_split(
                self.feature_matrix,
                features=["f0"],
                label_column_name="label",
                random_state=1.0,
            )

    def test_k_fold_split_output(self) -> None:
        """
        Test whether the output of k_fold_split() are correct.
        """
        # Test if the output are correct when n_samples % n_split ==0
        X = pd.DataFrame(
            {
                "f0": [6, 5, 4, 3, 2, 1],
                "f1": [1, 2, 3, 4, 5, 6],
                "label": [1, 0, 1, 0, 1, 1],
            }
        )
        folds = k_fold_split(X, features=["f0", "f1"], n_splits=3, shuffle_data=False)
        folds = [
            (
                fold[0].values.tolist(),
                fold[1].tolist(),
                fold[2].values.tolist(),
                fold[3].tolist(),
            )
            for fold in folds
        ]
        self.assertEqual(
            folds,
            [
                (
                    [[4, 3], [3, 4], [2, 5], [1, 6]],
                    [1, 0, 1, 1],
                    [[6, 1], [5, 2]],
                    [1, 0],
                ),
                (
                    [[6, 1], [5, 2], [2, 5], [1, 6]],
                    [1, 0, 1, 1],
                    [[4, 3], [3, 4]],
                    [1, 0],
                ),
                (
                    [[6, 1], [5, 2], [4, 3], [3, 4]],
                    [1, 0, 1, 0],
                    [[2, 5], [1, 6]],
                    [1, 1],
                ),
            ],
        )

        # Test if the output are correct when n_samples % n_split !=0
        X = pd.DataFrame(
            {
                "f0": [7, 6, 5, 4, 3, 2, 1],
                "f1": [0, 1, 2, 3, 4, 5, 6],
                "label": [0, 1, 0, 1, 0, 1, 1],
            }
        )
        folds = k_fold_split(X, features=["f0", "f1"], n_splits=3, shuffle_data=False)
        folds = [
            (
                fold[0].values.tolist(),
                fold[1].tolist(),
                fold[2].values.tolist(),
                fold[3].tolist(),
            )
            for fold in folds
        ]
        self.assertEqual(
            folds,
            [
                (
                    [[4, 3], [3, 4], [2, 5], [1, 6]],
                    [1, 0, 1, 1],
                    [[7, 0], [6, 1], [5, 2]],
                    [0, 1, 0],
                ),
                (
                    [[7, 0], [6, 1], [5, 2], [2, 5], [1, 6]],
                    [0, 1, 0, 1, 1],
                    [[4, 3], [3, 4]],
                    [1, 0],
                ),
                (
                    [[7, 0], [6, 1], [5, 2], [4, 3], [3, 4]],
                    [0, 1, 0, 1, 0],
                    [[2, 5], [1, 6]],
                    [1, 1],
                ),
            ],
        )

        # Test if the output are correct when shuffe == True.
        folds_unshuffled = k_fold_split(
            X, features=["f0", "f1"], n_splits=3, shuffle_data=False
        )
        folds_shuffled = k_fold_split(
            X, features=["f0", "f1"], n_splits=3, shuffle_data=True
        )
        self.assertNotEqual(
            [
                (
                    fold[0].values.tolist(),
                    fold[1].tolist(),
                    fold[2].values.tolist(),
                    fold[3].tolist(),
                )
                for fold in folds_unshuffled
            ],
            [
                (
                    fold[0].values.tolist(),
                    fold[1].tolist(),
                    fold[2].values.tolist(),
                    fold[3].tolist(),
                )
                for fold in folds_shuffled
            ],
        )

        # Test whether the shuffed results are reproducible.
        folds_1 = k_fold_split(
            X, features=["f0", "f1"], n_splits=3, shuffle_data=True, random_state=123
        )
        folds_2 = k_fold_split(
            X, features=["f0", "f1"], n_splits=3, shuffle_data=True, random_state=123
        )
        self.assertEqual(
            [
                (
                    fold[0].values.tolist(),
                    fold[1].tolist(),
                    fold[2].values.tolist(),
                    fold[3].tolist(),
                )
                for fold in folds_1
            ],
            [
                (
                    fold[0].values.tolist(),
                    fold[1].tolist(),
                    fold[2].values.tolist(),
                    fold[3].tolist(),
                )
                for fold in folds_2
            ],
        )


def main() -> None:
    unittest.main()


if __name__ == "__main__":
    main()
