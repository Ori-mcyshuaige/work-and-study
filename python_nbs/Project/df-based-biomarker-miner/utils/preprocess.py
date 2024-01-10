###############################################################
## Copyright: GeneGenieDx Corp 2022
# Author: whgu
# Date of creation: 04/04/2022
# Date of revision: 05/24/2022
#
## AIM
## Description: Functions to preprocess the feature matrix,
# prepare train and test data and k-fold split.
#
###############################################################
from typing import Tuple, List
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.utils import shuffle
from sklearn.model_selection import KFold

import numpy as np
import pandas as pd


def preprocess(
    features: np.ndarray,
    impute_missing_value: bool = True,
    impute_strategy: str = "mean",
    scale_by_minmax: bool = False,
    scale_range: Tuple[float, float] = (0, 1),
) -> np.ndarray:
    """
    Preprocess missing values and apply feature normalization.

    There are two ways address missing feature (column) values:
        a. Impute by the mean, median or mode.
        b. Drop features that contain missing values.

    There are also two ways for feature (column) normalization:
        a. MinMaxScaler: scale each feature to a predefined range.
        b. StandardScaler: convert feature to N(0, 1). The premise is that the feature distribution obey a normal distribution.

    Params:
        features (numpy matrix of shape (n_samlpes, n_features)): Raw features.
        impute_missing_value (bool = True):
            True: impute the missing values with impute_strategy.
            False: then drop columns with missing values.
        scale_by_minmax (bool = True): Scaler type.
            True: scale the value interval of the features to a given range.
            False: convert the features into a N(0, 1).
        scale_range (Tuple[float, float] = (0, 1)):
            specify desired range of transformed data for MinMaxScaler.
        impute_strategy (str = "mean"): by column
            "mean": replace missing values using the mean.
            "median": replace missing values using the median.
            "most_frequent": replace missing using the most frequent value.

    Return:
        preprocessed_features (numpy matrix of shape (n_samlpes, n_features)):
            Feature matrix after preprocessing.
    """
    # Check if input are valid
    if not isinstance(features, np.ndarray):
        raise TypeError("features should be of type np.ndarray!")
    assert impute_strategy in [
        "mean",
        "median",
        "most_frequent",
    ], "Invalid imputation strategy!"

    if impute_missing_value:
        # Impute missing values with specified strategy.
        imp = SimpleImputer(missing_values=np.nan, strategy=impute_strategy)
        preprocessed_features = imp.fit_transform(features)
    else:
        # Drop columns with missing values
        preprocessed_features = pd.DataFrame(features).dropna(axis=1).values

    if scale_by_minmax:
        # Min-max normalization.
        preprocessed_features = MinMaxScaler(feature_range=scale_range).fit_transform(
            preprocessed_features
        )
    else:
        # Standardization.
        preprocessed_features = StandardScaler().fit_transform(preprocessed_features)

    return preprocessed_features


def prepare_train_test_data(
    data: pd.DataFrame,
    features: List[str],
    label_column_name: str = "label",
    test_ratio: float = 0.2,
    shuffle_data: bool = True,
    random_state: int = None,
) -> Tuple[pd.DataFrame, np.ndarray, pd.DataFrame, np.ndarray]:
    """
    Split DataFrame into train and test subsets containing the specified
    features. By default, shuffle the data (shuffle_data = True) before
    splitting. Note that the shuffled data are copies of the input. The
    input `data` is not modified.

    Params:
        data (pd.Dataframe): input data.
        features (List[str]): features (columns) to be considered.
        label_column_name (str = "label"): The name of the column specifying the class labels.
        test_ratio (float = 0.2): Proportion of the dataset to include in the test split.
        shuffle_data (bool = True):  Shuffle the data before splitting or not.
        random_state (int = None): Pass an int for reproducible results across multiple function calls.
    Return:
        train_X (DataFrame of shape (n_samlpes, n_features)): Features of train set.
        train_y (numpy matrix of shape (n_samlpes, )): Labels of train set.
        test_X (DataFrame of shape (n_samlpes, n_features)): Features of test set.
        test_y (numpy matrix of shape (n_samlpes, )): Labels of test set.
    """
    # Check if input are valid
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data should be of type pd.DataFrame!")

    if not isinstance(features, List):
        raise TypeError("features should be of type List!")

    if test_ratio < 0 or test_ratio > 1:
        raise ValueError("test_ratio should be in [0,1]!")

    if random_state and not isinstance(random_state, int):
        raise TypeError("random_state should be of type int!")

    # Check if passed label column name and features exist in the column name
    # of data.
    column_names = data.columns.values.tolist()
    assert (
        label_column_name in column_names
    ), f"{label_column_name} doesn't exist in the column name of data."
    for feature in features:
        assert (
            feature in column_names
        ), f"{feature} doesn't exist in the column name of data."

    if shuffle_data:
        # Shuffle the data before splitting.
        data = shuffle(data, random_state=random_state)

    # Split into train and test subsets.
    train_X = data.iloc[int(len(data) * test_ratio) :].loc[:, features]
    train_y = data.iloc[int(len(data) * test_ratio) :].loc[:, label_column_name].values
    test_X = data.iloc[: int(len(data) * test_ratio)].loc[:, features]
    test_y = data.iloc[: int(len(data) * test_ratio)].loc[:, label_column_name].values

    return train_X, train_y, test_X, test_y


def k_fold_split(
    data: pd.DataFrame,
    features: List[str],
    label_column_name: str = "label",
    n_splits: int = 3,
    shuffle_data: bool = True,
    random_state: int = None,
) -> List[Tuple[pd.DataFrame, np.ndarray, pd.DataFrame, np.ndarray]]:
    """
    Split dataset into k == n_split consecutive folds containing the specified features.
    Each fold is then used once as a validation while the k - 1 remaining folds form the training set.
    If shuffle_data == True, only generated indices are shuffled and the input `data` is not modified.
    Note that the first n_samples % n_splits folds have size n_samples // n_splits + 1,
    other folds have size n_samples // n_splits, where n_samples is the number of samples.

    Params:
        data (pd.Dataframe): input data.
        features (List[str]): features (columns) to be considered.
        label_column_name (str = "label"): The name of the column specifying the class labels.
        n_split (int = 3): Number of folds. Must be at least 2.
        shuffle_data (bool = True):  Whether to shuffle the data before splitting into batches.
        random_state (int = None): Pass an int for reproducible results across multiple function calls.
    Returns:
        folds (List[Tuple[pd.DataFrame, np.ndarray, pd.DataFrame, np.ndarray]]): A list of (train_X, train_y, test_X and test_y).
    """
    # Check if input are valid
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data should be of type pd.DataFrame!")

    if not isinstance(features, List):
        raise TypeError("features should be of type List!")

    if random_state and not isinstance(random_state, int):
        raise TypeError("random_state should be of type int!")

    # Check if passed label column name and features exist in the column name
    # of data.
    column_names = data.columns.values.tolist()
    assert (
        label_column_name in column_names
    ), f"{label_column_name} doesn't exist in the column name of data."
    for feature in features:
        assert (
            feature in column_names
        ), f"{feature} doesn't exist in the column name of data."

    X = data.loc[:, features]
    y = data.loc[:, label_column_name]
    kf = KFold(n_splits=n_splits, shuffle=shuffle_data, random_state=random_state)

    folds = []
    # Generate indices to split data into training and test set.
    for train_index, test_index in kf.split(X):
        train_X, test_X = X.iloc[train_index], X.iloc[test_index]
        train_y, test_y = y.iloc[train_index], y.iloc[test_index]
        folds.append((train_X, train_y.values, test_X, test_y.values))

    return folds