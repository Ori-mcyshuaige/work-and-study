###############################################################
## Copyright: GeneGenieDx Corp 2022
# Author: whgu
# Date of creation: 04/26/2022
# Date of revision: 05/14/2022
#
## AIM
## Description: features can be of different types, e.g. protein,
# RNA-expression. This file calculates the importance
# of each feature type and selecting the best-performing combination of feature
# types
###############################################################
from itertools import combinations
from model import *


def get_feature_type_combination(feature_types: List[str],) -> List[Tuple[str,]]:
    """
    Enumerate all combinations of feature types.

    Param:
        feature_types (List[str]): All feature types.
    Return:
        feature_types_combinations (List[Tuple[str,]]): All combinations of feature types.
    """
    # Check if input are valid
    if not isinstance(feature_types, list):
        raise TypeError("feature_types should be of type list!")

    # Deduplication
    feature_types = list(set(feature_types))

    feature_types_combinations = []
    for i in range(1, len(feature_types) + 1):
        for c in combinations(feature_types, i):
            feature_types_combinations.append(c)

    return feature_types_combinations


def calculate_feature_type_importance(
    feature_matrix: pd.DataFrame,
    feature_type_map: Dict[str, str],
    model_type: ModelType = ModelType.SVM,
    metrics: Metrics = Metrics.Accuracy,
    test_ratio: float = 0.2,
    cv: int = 5,
    n_jobs: int = 8,
) -> List[Tuple[str, float]]:
    """
    Calculate the importance of each feature type.
    Detailed steps:
    1. Split feature_matrix into training set and test set.
    2. Select all feature types of training set.
    3. Train model (including parameter optimization via CV on training set)
        and test the performance on test set.
    4. Then remove one feature type at a time to see the impact on performance.
        The importance of feature_type_i is defined as:
        tanh(1 - performance_without_feature_type_i/performance_with_all_features).
    The model is specified by "model_type" and the performance is measured by "metrics".

    Params:
        feature_matrix (pd.DataFrame): DataFrame contains features, feature column names, and sample label.
        feature_type_map (Dict[str, str]): feature name to feature type map.
        model_type (ModelType): defined in utils.ModelType,
            one of [SVM, LR, RF, KNN, DT, GBDT].
        metrics (Metrics = Metrics.Accuracy): Metrics to measure model's performance.
        test_ratio (float = 0.2): Proportion of the dataset to include in the test split.
        cv (int = 5): Specify the number of folds in a KFold when searching the hyper-parameter space for the best cross validation.
        n_job (int = 8): Number of jobs to run in parallel when search the hyper-parameter space for the best cross validation score.
    Return:
        feature_type_importance (List[Tuple[str, float]]): The importance of each feature type.
    """
    # Check if all feature_type_map's features exist in the column name of feature_matrix.
    for feature in feature_type_map.keys():
        assert (
            feature in feature_matrix.columns.values
        ), f"{feature} is not feature_matrix!"
    if feature_type_map == {}:
        raise ValueError("feature_type_map cannot be empty!")

    feature_types = list(set(feature_type_map.values()))
    # If feature_types contains only one feature type, return the result directly
    if len(feature_types) == 1:
        return [(feature_types[0], 1.0)]

    # Train model with all features.
    train_X, train_y, test_X, test_y = prepare_train_test_data(
        feature_matrix,
        list(feature_type_map.keys()),
        test_ratio=test_ratio,
        shuffle_data=False,
    )
    classifier, performance_with_all_features, feature_importance = train(
        train_X,
        train_y,
        test_X,
        test_y,
        model_type=model_type,
        metrics=metrics,
        cv=cv,
        n_jobs=n_jobs,
    )

    feature_type_importance = []
    for feature_type in feature_types:
        # Remove one feature type at a time to see the impact on performance.
        features = [
            key for key in feature_type_map if feature_type_map[key] != feature_type
        ]
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            feature_matrix, features, shuffle_data=False,
        )
        classifier, performance, feature_importance = train(
            train_X,
            train_y,
            test_X,
            test_y,
            model_type=model_type,
            metrics=metrics,
            cv=cv,
            n_jobs=n_jobs,
        )

        importance = 1 - performance / performance_with_all_features
        feature_type_importance.append((feature_type, importance))

    return feature_type_importance


def select_best_feature_type_combination(
    feature_matrix: pd.DataFrame,
    feature_type_map: Dict[str, str],
    model_type: ModelType = ModelType.SVM,
    metrics: Metrics = Metrics.Accuracy,
    cv: int = 5,
    n_jobs: int = 8,
):
    """
    Select the best-performing feature type combination by enumerating all
    combinations of feature types for a given `model_type`

    Detailed steps:
    1. Enumerate all combinations of feature types.
    2. Split feature_matrix into training set and test set.
    3. For each feature types' combination:
        - Select features of training set belonging to feature types' combination.
        - Train model (including parameter optimization via CV on training set) and test the performance on test set.
    The model is specified by "model_type" and the performance is measured by "metrics". Note that if the performance of two combinations is the same,
    we will choose the one with less number of feature types.

    Params:
        feature_matrix (pd.DataFrame): DataFrame contains features, feature column names, and sample label.
        feature_type_map (Dict[str, str]): feature name to feature type map.
        model_type (ModelType): defined in utils.ModelType,
            one of [SVM, LR, RF, KNN, DT, GBDT].
        metrics (Metrics = Metrics.Accuracy): Metrics to measure model's performance.
        cv (int = 5): Specify the number of folds in a KFold when searching the
        hyper-parameter space for the best cross validation.
        n_job (int = 8): Number of jobs to run in parallel when search the
        hyper-parameter space for the best cross validation score.
    Returns:
        best_feature_types_combination (Tuple[str]): Selected feature types with best performance on specified model.
        best_performance (float): The performance according to the given metrics on the selected feature types.
    """
    # Check if all feature_type_map's features exist in the column name of feature_matrix.
    for feature in feature_type_map.keys():
        assert (
            feature in feature_matrix.columns.values
        ), f"{feature} is not feature_matrix!"

    feature_types = list(set(feature_type_map.values()))
    best_performance = 0
    best_feature_types_combination = tuple()

    # Enumerate all combinations of feature types and test the performance.
    feature_types_combinations = get_feature_type_combination(feature_types)
    for combination in feature_types_combinations:
        # Select features belong to the feature types.
        features = [
            key for key in feature_type_map if feature_type_map[key] in combination
        ]
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            feature_matrix, features, shuffle_data=False,
        )
        classifier, performance, feature_importance = train(
            train_X,
            train_y,
            test_X,
            test_y,
            model_type=model_type,
            metrics=metrics,
            cv=cv,
            n_jobs=n_jobs,
        )

        if performance > best_performance:
            best_performance = performance
            best_feature_types_combination = combination
        elif performance == best_performance and len(combination) < len(
            best_feature_types_combination
        ):
            best_feature_types_combination = combination

    return best_feature_types_combination, best_performance
