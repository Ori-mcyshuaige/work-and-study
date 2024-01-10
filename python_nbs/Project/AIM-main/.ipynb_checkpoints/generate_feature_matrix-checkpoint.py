#!/usr/bin/env python3
###############################################################
# Copyright: GeneGenieDx Corp 2022
# Author: whgu, xyang
# Date of creation: 04/03/2022
# Date of revision: 04/19/2022
## AIM
## Description: This file is to prepare mock data (feature matrix)
# for classification.
#
###############################################################
from sklearn.datasets import make_blobs
from typing import List, Dict
import pandas as pd


def generate_mock_feature_matrix(
    n_samples: int = 150,
    n_classes: int = 3,
    feature_count_per_type: List[int] = [100, 150, 200],
    cluster_std: float = 1.0,
    random_state: int = 2,
) -> (pd.DataFrame, Dict[str, str]):
    """
    Generate a mock 2d feature matrix.

    Params:
        n_samlpes (int): The total number of samples (rows) to generate, default 150.
        n_classes (int): Number of classes default 3.
        feature_count_per_type (List[int]): count of features per each type,
            this is to mock how many features to generate per type. Features are given names of `f{#}` and feature types are given names of `t{#}`
        cluster_std (float): Variance for each class, default = 1.0.
        random_state (int): Random number seed, default = 2.
    Returns:
        feature_matrix (DataFrame): contains features, feature column names, and sample label, e.g.
                f0        f1        f2        f3        f4        f5  label
        0 -6.245849  2.149235 -4.644562 -5.851067  1.001459  0.429347      1
        1 -2.027973 -9.472450  0.115142 -1.449986 -1.336074 -4.382083      0

        and feature name to feature type map
    """
    # Check if input are valid
    assert n_classes > 0, "n_class should be greater than 0!"
    if not isinstance(n_classes, int):
        raise TypeError("n_class should be of type int!")
    if not isinstance(random_state, int):
        raise TypeError("random_state should be of type int!")

    n_features = sum(feature_count_per_type)
    X, y = make_blobs(
        n_samples=n_samples,
        centers=n_classes,
        n_features=n_features,
        cluster_std=cluster_std,
        random_state=random_state,
    )

    feature_names: List[str] = []
    feature_to_type: Dict[str, str] = {}
    cnt = 0
    for i, c in enumerate(feature_count_per_type):
        for j in range(c):
            feature_name = f"f{cnt}"
            cnt += 1
            feature_names.append(feature_name)
            feature_to_type[feature_name] = f"t{i}"

    feature_matrix = pd.DataFrame(columns=feature_names, data=X)
    feature_matrix["label"] = y

    return feature_matrix, feature_to_type


def main() -> None:
    df, feature_to_type = generate_mock_feature_matrix(
        n_samples=2, feature_count_per_type=[1, 2, 3]
    )
    print(f"feature to type: {feature_to_type}")
    print(f"Mock matrix:\n{df}")


if __name__ == "__main__":
    main()
