# AI-Module (AIM)
```
Copyright: GeneGenieDx Corp 2021

Author: Wenhao Gu

Date of creation: 05/09/2022
```

The AIM mainly consists of feature engineering layer and machine learning layer. The goal of the AIM is to automatically apply a set of supervised machine learning methods on the input feature matrix data and generate the best model and features to explain the biological data. The main functions of this library include: 

- `generate_mock_feature_matrix()`
- `prepare_train_test_data()`
- `train()`
- `calculate_feature_type_importance()`
- `select_best_feature_type_combination()`

Some possible applications of this library include:

- Identification of the best model
- Select best feature type combination

For more information, please refer to the design document API: <https://docs.google.com/document/d/1P94eI1n1TU4psuINpQFrBPAb32hLfAv8>

# Examples


```python
from generate_feature_matrix import *
from select_feature_types import *
```

## Generate a mock 2d feature matrix.


```python
feature_matrix, feature_to_type = generate_mock_feature_matrix(
    feature_count_per_type=[4, 5, 6]
)

print(feature_matrix)
```

               f0         f1        f2        f3        f4        f5        f6  \
    0   -1.228331  -4.722796  3.987628  0.295844  0.541162 -3.080539  6.114187   
    1    6.422858   6.809397  1.317255  6.063358 -8.812726  1.057632 -7.839650   
    2    5.478865   6.311221  0.291409  7.578688 -8.612007 -0.202509 -7.576588   
    3   -0.238763  -8.506257  1.745384  0.546023 -2.063949 -3.932844 -6.042063   
    4   -2.121215  -8.747930  1.162462 -0.772695 -1.043783 -3.341773 -4.519976   
    ..        ...        ...       ...       ...       ...       ...       ...   
    145 -2.250645 -10.413647 -0.280318 -1.023522 -2.981393 -1.303268 -6.723011   
    146  6.402470   6.919573 -0.498770  8.590848 -6.672449 -0.267820 -7.551561   
    147  0.289417  -6.975000  0.847864 -2.273653 -1.108150 -4.615332 -4.447003   
    148 -0.025410  -6.010656  2.254331 -1.806333 -0.061683 -2.610419  6.413356   
    149  0.679093  -6.767254  3.387889 -0.085110  2.010616 -2.096499  5.390148   
    
               f7        f8        f9        f10       f11       f12       f13  \
    0    1.209640 -7.981336  2.414359  10.043919 -0.842776  6.561039 -4.103456   
    1   -1.068744 -8.386696 -7.108036   1.575386 -5.136713 -8.756422 -6.705829   
    2   -2.923908 -8.186970 -5.786687   1.742523 -4.882675 -7.763267 -4.217732   
    3    1.209023 -2.573917 -6.338006   2.801335 -0.068471 -8.254790  0.212586   
    4    4.880607 -3.946544 -5.255814   1.695118  1.933839 -6.646939  2.121843   
    ..        ...       ...       ...        ...       ...       ...       ...   
    145  2.364564 -4.115199 -4.569754   1.571396 -1.571885 -7.429315 -0.704815   
    146 -1.331050 -8.425862 -8.769228   1.668455 -4.557114 -6.171032 -5.397968   
    147  2.093372 -4.152274 -5.033793   1.416254  0.504104 -6.784677  2.432565   
    148  2.214427 -6.573489  4.390073   9.527225  1.085975  9.093339 -1.517997   
    149  1.453707 -6.752029  4.305121   8.951805  0.854286  6.814927 -2.230908   
    
              f14  label  
    0   -1.203397      2  
    1   -2.943550      1  
    2   -2.257172      1  
    3   -5.981733      0  
    4   -6.081805      0  
    ..        ...    ...  
    145 -7.905434      0  
    146 -3.521183      1  
    147 -6.533884      0  
    148  1.721513      2  
    149  0.101926      2  
    
    [150 rows x 16 columns]


## Split dataset


```python
selected_features = [
    feature
    for feature in feature_to_type.keys()
    if feature_to_type[feature] in ["t1", "t2"]
]
train_X, train_y, test_X, test_y = prepare_train_test_data(
    feature_matrix, selected_features, test_ratio=0.2
)

print("Size of train set:", len(train_X))
print("Size of test set:", len(test_X))
```

    Size of train set: 120
    Size of test set: 30


## Train and test model
Include parameter optimization via CV on training set. 
The performance is measured on test set.


```python
accuracy, feature_importance = train(
    train_X, train_y, test_X, test_y, model_type=ModelType.SVM, metrics=Metrics.Accuracy
)
print(accuracy, feature_importance)
```

    0.7666666666666667 [-0.05333333  0.02       -0.06        0.04666667  0.05333333 -0.1
      0.         -0.02       -0.06666667  0.05333333  0.03333333]


##  Identify the best model


```python
for model_type in ModelType:
    accuracy, feature_importance = train(
        train_X,
        train_y,
        test_X,
        test_y,
        model_type=model_type,
        metrics=Metrics.Accuracy,
        n_jobs=20,
    )
    print(model_type, accuracy)
```

    ModelType.SVM 0.7666666666666667
    ModelType.LR 0.7333333333333333
    ModelType.RF 0.7666666666666667
    ModelType.KNN 0.7
    ModelType.DT 0.4
    ModelType.GBDT 0.26666666666666666


## Select best feature type combination


```python
best_feature_types_combination, best_performance = select_best_feature_type_combination(
    feature_matrix, feature_to_type, model_type=ModelType.SVM, metrics=Metrics.Accuracy
)

print(best_feature_types_combination, best_performance)
```

    ('t0', 't2') 1.0


# API

###  `generate_mock_feature_matrix(n_samples, n_classes, feature_count_per_type, cluster_std, random_state)`

Generate a mock 2d feature matrix.

The parameters to the function are as follows:

- `n_samlpes` (int): The total number of samples (rows) to generate, default 150.
- `n_classes` (int): Number of classes default 3.
- `feature_count_per_type` (List[int]): count of features per each type,
- this is to mock how many features to generate per type. Features are given names of `f{#}` and feature types are given names of `t{#}`
- `cluster_std` (float): Variance for each class, default = 1.0.
- `random_state` (int): Random number seed, default = 2.

Returns:

- feature_matrix (DataFrame): contains features, feature column names, and sample label, e.g.

```
            f0        f1        f2        f3        f4        f5  label
    0 -6.245849  2.149235 -4.644562 -5.851067  1.001459  0.429347      1
    1 -2.027973 -9.472450  0.115142 -1.449986 -1.336074 -4.382083      0
```

- feature name to feature type map.

### `prepare_train_test_data(data, features, label_column_name, test_ratio, shuffle_data, random_state)`

Split DataFrame into train and test subsets containing the specified features. By default, shuffle the data (shuffle_data = True) before splitting. Note that the shuffled data are copies of the input. The input `data` is not modified.

The parameters to the function are as follows:

- data (pd.Dataframe): input data.
- features (List[str]): features (columns) to be considered.
- `label_column_name` (str = "label"): The name of the column specifying the class labels.
- `test_ratio` (float = 0.2): Proportion of the dataset to include in the test split.
- `shuffle_data` (bool = True):  Shuffle the data before splitting or not.
- `random_state` (int = None): Pass an int for reproducible results across multiple function calls.

Returns:

- `train_X` (numpy matrix of shape (`n_samlpes, n_features`)): Features of train set.
- `train_y` (numpy matrix of shape (`n_samlpes`, 1)): Labels of train set.
- `test_X` (numpy matrix of shape (`n_samlpes, n_features`)): Features of test set.
- `test_y` (numpy matrix of shape (`n_samlpes`, 1)): Labels of test set.

### `preprocess(features, impute_missing_value, impute_strategy, scale_by_minmax, scale_range)`

Preprocess missing values and apply feature normalization.

There are two ways address missing feature (column) values:

- Impute by the mean, median or mode.
- Drop features that contain missing values.

There are also two ways for feature (column) normalization:

- MinMaxScaler: scale each feature to a predefined range.
- StandardScaler: convert feature to N(0, 1). The premise is that the feature distribution obey a normal distribution.

The parameters to the function are as follows:

- features (numpy matrix of shape (n_samlpes, n_features)): Raw features.
- `impute_missing_value` (bool = True):
            True: impute the missing values with impute_strategy.
            False: then drop columns with missing values.
- `scale_by_minmax` (bool = True): Scaler type.
            True: scale the value interval of the features to a given range.
            False: convert the features into a N(0, 1).
- `scale_range` (Tuple[float, float] = (0, 1)):
            specify desired range of transformed data for MinMaxScaler.
- `impute_strategy` (str = "mean"): by column
            "mean": replace missing values using the mean.
            "median": replace missing values using the median.
            "most_frequent": replace missing using the most frequent value.

Return:

- `preprocessed_features` (numpy matrix of shape (`n_samlpes, n_features`)): Feature matrix after preprocessing.

### `train(train_X, train_y, test_X, test_y, model_type, metrics, average, cv, n_jobs)`

This function uses parameter optimization and obtain model with the best performance via CV on training set. Then return the results on test set and permutation feature importance.

The parameters to the function are as follows:

- train_X (np.ndarray of shape (`n_samlpes, n_features`)): Features of train set.
- train_y (np.ndarray of shape (`n_samlpes`, 1)): Labels of train set.
- test_X (np.ndarray of shape (`n_samlpes, n_features`)): Features of test set.
- test_y (np.ndarray of shape (`n_samlpes`, 1)): Labels of test set.
- model_type (ModelType): Specified model type, one of (ModelType.SVM, ModelType.LR, ModelType.RF, ModelType.KNN, ModelType.DT, ModelType.GBDT).
- metrics (Metrics = Metrics.Accuracy): Way to measure model's performance.
- average (Average = Average.Micro): When the metrics in (Metrics.Fscore, Metrics.Precision, Metrics.Recall), average determines the type of averaging performed on the data.
- cv (int = 5): Specify the number of folds in a KFold when search the hyper-parameter space for the best cross validation score.
- `n_job` (int = 8): Number of jobs to run in parallel when search the hyper-parameter space for the best cross validation score.

Returns:

- model's performance (float): Model's performance on test set measured by specified metrics.
- permutation feature importance (np.ndarray of shape (`n_features`, 1)): Permutation importance for feature evaluation.

### `k_fold_cross_validation(X, y, classifier, metrics, average, cv, n_jobs)`

This function evaluates model's performance by cross-validation and return performance of the model for each run of the cross validation. The performance is measured according to the given `metrics` and `average`.

The parameters to the function are as follows:

- X (np.ndarray of shape (`n_samlpes, n_features`)): Training vectors.
- y (np.ndarray of shape (`n_samlpes`, 1)): Class labels.
- classifier (estimator object): Scikit-learn classifier.
- metrics (Metrics = Metrics.Accuracy): Metrics to measure model's performance.
- average (Average = Average.Micro): When the metrics in (Metrics.Fscore, Metrics.Precision, Metrics.Recall), average determines the type of averaging performed on the data.
- cv (int = 5): Specify the number of folds in a KFold.
- `n_job` (int = 8): Number of jobs to run in parallel.

Return:

- `k_fold_performance` (np.ndarray of float of shape=cv): Array of performance of the model for each run of the cross validation.

### `select_best_feature_type_combination(feature_matrix, feature_type_map, model_type, metrics, cv, n_jobs)`

Select the best-performing feature type combination by enumerating all combinations of feature types for a given `model_type`

The parameters to the function are as follows:

- `feature_matrix` (pd.DataFrame): DataFrame contains features, feature column names, and sample label.
- `feature_type_map` (Dict[str, str]): feature name to feature type map.
- `model_type` (ModelType): defined in utils.ModelType, one of (SVM, LR, RF, KNN, DT, GBDT).
- metrics (Metrics = Metrics.Accuracy): Metrics to measure model's performance.
- cv (int = 5): Specify the number of folds in a KFold when searching the hyper-parameter space for the best cross validation.
- n_job (int = 8): Number of jobs to run in parallel when search the hyper-parameter space for the best cross validation score.

Returns:

- `best_feature_types_combination` (Tuple[str]): Selected feature types with best performance on specified model.
- `best_performance` (float): The performance according to the given metrics on the selected feature types.

### `calculate_feature_type_importance(feature_matrix, feature_type_map, model_type, metrics, test_ratio, cv, n_jobs)` 

Calculate the importance of each feature type.

The parameters to the function are as follows:

- `feature_matrix` (pd.DataFrame): DataFrame contains features, feature column names, and sample label.
- `feature_type_map` (Dict[str, str]): feature name to feature type map.
- `model_type` (ModelType): defined in utils.ModelType,one of (SVM, LR, RF, KNN, DT, GBDT).
- metrics (Metrics = Metrics.Accuracy): Metrics to measure model's performance.
- `test_ratio` (float = 0.2): Proportion of the dataset to include in the test split.
- cv (int = 5): Specify the number of folds in a KFold when searching the hyper-parameter space for the best cross validation.
- n_job (int = 8): Number of jobs to run in parallel when search the hyper-parameter space for the best cross validation score.

Return:

- `feature_type_importance` (List[Tuple[str, float]]): The importance of each feature type.

___

Confidentiality Notice: This document is confidential and contains proprietary information and intellectual property of GeneGenieDx, Inc. Except as required by federal, state or local laws or regulations, neither this document nor any of the information contained herein may be reproduced or disclosed under any circumstances without the express written permission of GeneGenieDx, Inc.