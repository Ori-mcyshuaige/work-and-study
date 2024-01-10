###############################################################
## Copyright: GeneGenieDx Corp 2022
# Author: whgu
# Date of creation: 04/04/2022
# Date of revision: 06/24/2022
## AIM
## Description: Functions of training and testing models.
#
## usage
#   import model
###############################################################
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, cross_val_score
from skopt import BayesSearchCV
from sklearn.inspection import permutation_importance

from preprocess import *
from utils import *


def calculate_prediction_performance(
    y_true,
    y_pred,
    metrics: Metrics = Metrics.Accuracy,
    average: Average = Average.Micro,
) -> float:
    """
    Implement accuracy, fscore, precision and recall to measure classification performance.
    These implementations are work in both binary classification and multilabel case.

    Parmas:
        y_true (1d array-like): Ground truth (correct) target values.
        y_pred (1d array-like): Estimated targets as returned by a classifier.
        metrics (Metrics = Metrics.Accuracy): Metrics to measure classification performance.
                If Metrics.Accuracy, then compute accuracy classification score.
                If Metrics.Fscore, then compute the F1 score.
                If Metrics.Precision, then compute the precision.
                If Metrics.Recall, then compute the recall.
        average (Average = Average.Micro): When the metrics in [Metrics.Fscore, Metrics.Precision, Metrics.Recall],
                average determines the type of averaging performed on the data.
                If Average.Binary, then only report results for the class == 1.
                This is applicable only if targets (y_{true,pred}) are binary.
                If Average.Micro, then calculate metrics globally by counting the total true positives,
                false negatives and false positives.
                If Average.Macro, then calculate metrics for each label, and find their unweighted mean.
                This does not take label imbalance into account.
                If Average.Weighted, then calculate metrics for each label, and find their average weighted by
                support (the number of true instances for each label).
    Return:
        accuracy, fscore, precision or recall (float): Calculated result.
    """
    # Check if input are valid.
    if not isinstance(metrics, Metrics):
        raise TypeError("metrics should be of type Metrics!")
    if not isinstance(average, Average):
        raise TypeError("average should be of type Average!")

    if metrics == Metrics.Accuracy:
        return accuracy_score(y_true, y_pred)
    elif metrics == Metrics.Fscore:
        return f1_score(y_true, y_pred, average=average.value)
    elif metrics == Metrics.Precision:
        return precision_score(y_true, y_pred, average=average.value)
    elif metrics == Metrics.Recall:
        return recall_score(y_true, y_pred, average=average.value)


def train(
    train_X: pd.DataFrame,
    train_y: np.ndarray,
    test_X: pd.DataFrame,
    test_y: np.ndarray,
    model_type: ModelType,
    hpo_algorithm: HPOAlgorithm = HPOAlgorithm.GridSearch,
    metrics: Metrics = Metrics.Accuracy,
    average: Average = Average.Micro,
    cv: int = 5,
    n_jobs: int = 8,
    hpo_search_iter: int = 100,
    display_attribute: bool = True,
    display_cv_results: bool = False,
) -> Tuple[BaseEstimator, float, List[Tuple[str, float]]]:
    """
    This function uses parameter optimization and obtain model with the best
    performance via CV on training set. Then return the results on test set
    and permutation feature importance.

    Detailed steps:
    1. Train the specified model on training set via cross validation (cv=5)
        and return the parameters that achieved the best results. During
        training, hyper parameter optimization method (HPO) was used.
    2. Once the best hyper parameter set was obtained from the previous step,
        the model was initialized and fit again on the training set.
    3. Using the trained model from the previous step, obtain and return the
        the performance on the test set.
    4. Each feature column from the test set is permuted and the metric is evaluated again.
        The permutation importance is defined to be the difference between the baseline metric and metric
        from permutating the feature column (Mean value over n_repeats == 5).

    Parmas:
        train_X (DataFrame of shape (n_samlpes, n_features)): Features of train set.
        train_y (numpy matrix of shape (n_samlpes, )): Labels of train set.
        test_X (DataFrame of shape (n_samlpes, n_features)): Features of test set.
        test_y (numpy matrix of shape (n_samlpes, )): Labels of test set.
        model_type (ModelType): Specified model type, one of [ModelType.SVM, ModelType.LR, ModelType.RF,
        ModelType.KNN, ModelType.DT, ModelType.GBDT].
        hpo_algorithm (HPOAlgorithm = HPOAlgorithm.GridSearch): Approaches to parameter search,
            If HPOAlgorithm.GridSearch, then exhaustive search over specified parameter values for an estimator;
            If HPOAlgorithm.RandomSearch, then sample a given number of candidates from a parameter space with
                a specified distribution;
            If HPOAlgorithm.BayesianSearch, then apply Bayesian optimization over hyper parameters.
        hpo_search_iter（int = 100) : Number of parameter settings that are sampled when using RandomSearch or BayesianSearch
                hpo_search_iter trades off runtime vs quality of the solution.
        metrics (Metrics = Metrics.Accuracy): Way to measure model's performance.
        average (Average = Average.Micro): When the metrics in [Metrics.Fscore, Metrics.Precision, Metrics.Recall],
                average determines the type of averaging performed on the data.
        cv (int = 5): Specify the number of folds in a KFold when search the hyper-parameter space
                        for the best cross validation score.
        n_job (int = 8): Number of jobs to run in parallel when search the hyper-parameter space
                        for the best cross validation score.
        display_attribute (bool = True): Whether to display classifier's attributes.
        display_cv_results (bool = False): Whether to display the cross-validated score of different parameter settings.
    Return:
        model (BaseEstimator): Trained model.
        model's performance (float): Model's performance on test set measured by specified metrics.
        feature_importance_with_name (List[Tuple[str, float]]]): Permutation importance for feature evaluation.
    """
    # Check if input are valid.
    check_X_y(train_X, train_y)
    check_X_y(test_X, test_y)
    if not isinstance(model_type, ModelType):
        raise TypeError("model_type should be of type ModelType!")
    scoring = get_scoring_str(metrics, average)
    feature_names = train_X.columns.values.tolist()

    if model_type == ModelType.SVM:
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data.
        best_params, best_score, cv_results = HPO(
            train_X.values,
            train_y,
            SVC(),
            ParamConfig.get_param_config(model_type, hpo_algorithm),
            hpo_algorithm=hpo_algorithm,
            metrics=metrics,
            average=average,
            cv=cv,
            n_jobs=n_jobs,
            hpo_search_iter=hpo_search_iter,
        )
        # Initialize with best parameters.
        classifier = SVC(
            kernel=best_params["kernel"],
            C=best_params["C"],
            gamma=best_params["gamma"],
            probability=True,
        )

    elif model_type == ModelType.LR:
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data.
        best_params, best_score, cv_results = HPO(
            train_X.values,
            train_y,
            LogisticRegression(),
            ParamConfig.get_param_config(model_type, hpo_algorithm),
            hpo_algorithm=hpo_algorithm,
            metrics=metrics,
            average=average,
            cv=cv,
            n_jobs=n_jobs,
            hpo_search_iter=hpo_search_iter,
        )
        # Initialize with best parameters.
        classifier = LogisticRegression(
            penalty=best_params["penalty"],
            C=best_params["C"],
            multi_class=best_params["multi_class"],
            solver=best_params["solver"],
            max_iter=best_params["max_iter"],
            l1_ratio=best_params["l1_ratio"],
        )

    elif model_type == ModelType.RF:
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data.
        best_params, best_score, cv_results = HPO(
            train_X.values,
            train_y,
            RandomForestClassifier(),
            ParamConfig.get_param_config(model_type, hpo_algorithm),
            hpo_algorithm=hpo_algorithm,
            metrics=metrics,
            average=average,
            cv=cv,
            n_jobs=n_jobs,
            hpo_search_iter=hpo_search_iter,
        )
        # Initialize with best parameters.
        classifier = RandomForestClassifier(
            n_estimators=best_params["n_estimators"],
            bootstrap=best_params["bootstrap"],
            criterion=best_params["criterion"],
            max_depth=best_params["max_depth"],
            max_features=best_params["max_features"],
            min_samples_leaf=best_params["min_samples_leaf"],
        )

    elif model_type == ModelType.KNN:
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data.
        best_params, best_score, cv_results = HPO(
            train_X.values,
            train_y,
            KNeighborsClassifier(),
            ParamConfig.get_param_config(model_type, hpo_algorithm),
            hpo_algorithm=hpo_algorithm,
            metrics=metrics,
            average=average,
            cv=cv,
            n_jobs=n_jobs,
            hpo_search_iter=hpo_search_iter,
        )
        # Initialize with best parameters.
        classifier = KNeighborsClassifier(
            n_neighbors=best_params["n_neighbors"], weights=best_params["weights"]
        )

    elif model_type == ModelType.DT:
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data.
        best_params, best_score, cv_results = HPO(
            train_X.values,
            train_y,
            DecisionTreeClassifier(),
            ParamConfig.get_param_config(model_type, hpo_algorithm),
            hpo_algorithm=hpo_algorithm,
            metrics=metrics,
            average=average,
            cv=cv,
            n_jobs=n_jobs,
            hpo_search_iter=hpo_search_iter,
        )
        # Initialize with best parameters.
        classifier = DecisionTreeClassifier(
            criterion=best_params["criterion"],
            splitter=best_params["splitter"],
            max_depth=best_params["max_depth"],
            max_features=best_params["max_features"],
            min_samples_leaf=best_params["min_samples_leaf"],
        )

    elif model_type == ModelType.GBDT:
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data.
        best_params, best_score, cv_results = HPO(
            train_X.values,
            train_y,
            GradientBoostingClassifier(),
            ParamConfig.get_param_config(model_type, hpo_algorithm),
            hpo_algorithm=hpo_algorithm,
            metrics=metrics,
            average=average,
            cv=cv,
            n_jobs=n_jobs,
            hpo_search_iter=hpo_search_iter,
        )
        # Initialize with best parameters.
        classifier = GradientBoostingClassifier(
            n_estimators=best_params["n_estimators"],
            learning_rate=best_params["learning_rate"],
            loss=best_params["loss"],
            max_depth=best_params["max_depth"],
            max_features=best_params["max_features"],
            min_samples_leaf=best_params["min_samples_leaf"],
        )
    elif model_type == ModelType.MLP:
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data.
        best_params, best_score, cv_results = HPO(
            train_X.values,
            train_y,
            MLPClassifier(),
            ParamConfig.get_param_config(model_type, hpo_algorithm),
            hpo_algorithm=hpo_algorithm,
            metrics=metrics,
            average=average,
            cv=cv,
            n_jobs=n_jobs,
        )
        # Initialize with best parameters.
        classifier = MLPClassifier(
            hidden_layer_sizes=best_params["hidden_layer_sizes"],
            activation=best_params["activation"],
            alpha=best_params["alpha"],
            batch_size=best_params["batch_size"],
            learning_rate_init=best_params["learning_rate_init"],
            max_iter=best_params["max_iter"],
        )

    # Train the model.
    classifier.fit(
        train_X, train_y,
    )

    # Assess the performance.
    y_pred = classifier.predict(test_X)
    performance = calculate_prediction_performance(test_y, y_pred, metrics, average)

    # Permutation feature importance.
    feature_importance = permutation_importance(
        classifier, test_X, test_y, n_jobs=n_jobs, scoring=scoring
    ).importances_mean
    feature_importance_with_name = list(zip(feature_names, feature_importance))
    # Rank feature importance by decreasing order.
    feature_importance_with_name.sort(key=lambda t: t[1], reverse=True)

    # Display classifer's attributes and cv results.
    if display_attribute:
        display_classifier_attr(classifier, best_params)
    if display_cv_results:
        print(
            "Cross-validated score of different parameter settings:",
            pd.DataFrame(cv_results),
        )

    return classifier, performance, feature_importance_with_name


def HPO(
    X: np.ndarray,
    y: np.ndarray,
    classifier: BaseEstimator,
    param_config: List[Dict],
    hpo_algorithm: HPOAlgorithm = HPOAlgorithm.GridSearch,
    metrics: Metrics = Metrics.Accuracy,
    average: Average = Average.Micro,
    cv: int = 5,
    n_jobs: int = 8,
    hpo_search_iter: int = 100,
) -> Tuple[Dict, float, Dict]:
    """
    Hyperparameter optimization: search the hyper-parameter space for the best cross validation score.

    Three approache to parameter search are provided currently:
    GridSearchCV, exhaustively search over specified parameter values;
    RandomSearchCV, sample a given number of candidates from a parameter space;
    BayesianSearchCV, apply Bayesian optimization over hyperparameters.

    Params:
        X (np.ndarray of shape (n_samlpes, n_features)): Training vectors.
        y (np.ndarray of shape (n_samlpes, 1)): Class labels.
        classifier (estimator object): Scikit-learn classifier.
        param_config (List[Dict[str]]): A list of dictionary with parameters names as keys and lists
                                        of parameter settings to try as values.
        hpo_algorithm (HPOAlgorithm = HPOAlgorithm.GridSearch): Approaches to parameter search,
            If HPOAlgorithm.GridSearch, then exhaustive search over specified parameter values for an estimator;
            If HPOAlgorithm.RandomSearch, then sample a given number of candidates from a parameter space with
                a specified distribution;
            If HPOAlgorithm.BayesianSearch, then apply Bayesian optimization over hyper parameters.
        metrics (Metrics = Metrics.Accuracy): Metrics to measure model's performance.
        average (Average = Average.Micro): When the metrics in [Metrics.Fscore, Metrics.Precision, Metrics.Recall],
                average determines the type of averaging performed on the data.
        cv (int = 5): Specify the number of folds in a KFold.
        n_job (int = 8): Number of jobs to run in parallel.
        hpo_search_iter（int = 100) : Number of parameter settings that are sampled when using RandomSearch or BayesianSearch
                hpo_search_iter trades off runtime vs quality of the solution.
    Returns:
        best_params (Dict[str,]): Parameter setting that gave the best results on the hold out data.
        best_score (float): Mean cross-validated score of the best_estimator.
        cv_results (Dict[str,]): Cross-validated score of different parameter settings.
    """
    # Check if input are valid.
    if not isinstance(X, np.ndarray):
        raise TypeError("X should be of type np.ndarray!")
    if not isinstance(y, np.ndarray):
        raise TypeError("y should be of type np.ndarray!")
    if not isinstance(hpo_algorithm, HPOAlgorithm):
        raise TypeError("hpo_algorithm should be of type HPOAlgorithm!")
    if not isinstance(metrics, Metrics):
        raise TypeError("metrics should be of type Metrics!")
    if not isinstance(average, Average):
        raise TypeError("average should be of type Average!")
    assert hpo_search_iter > 0, f"{hpo_search_iter} must be greater than 0!"
    scoring = get_scoring_str(metrics, average)

    if hpo_algorithm == HPOAlgorithm.GridSearch:
        clf = GridSearchCV(
            classifier, param_config, cv=cv, scoring=scoring, n_jobs=n_jobs
        )
    elif hpo_algorithm == HPOAlgorithm.RandomSearch:
        clf = RandomizedSearchCV(
            classifier,
            param_config,
            cv=cv,
            scoring=scoring,
            n_jobs=n_jobs,
            n_iter=hpo_search_iter,
        )
    elif hpo_algorithm == HPOAlgorithm.BayesianSearch:
        clf = BayesSearchCV(
            classifier, param_config, scoring=scoring, n_jobs=n_jobs, n_iter=hpo_search_iter
        )

    clf.fit(X, y)
    best_params = clf.best_params_
    best_score = clf.best_score_
    cv_results = clf.cv_results_
    return best_params, best_score, cv_results


def k_fold_cross_validation(
    X: np.ndarray,
    y: np.ndarray,
    classifier: BaseEstimator,
    metrics: Metrics = Metrics.Accuracy,
    average: Average = Average.Micro,
    cv: int = 5,
    n_jobs: int = 8,
) -> np.ndarray:
    """
    This function evaluates model's performance by cross-validation and return
    performance of the model for each run of the cross validation.
    The performance is measured according to the given `metrics` and `average`.

    Detail steps:
    1. The input is split into k=cv smaller sets.
    2. A model is trained using k-1 of the folds as training data.
    3. The resulting model is validated on the remaining part of the data.

    Params:
        X (np.ndarray of shape (n_samlpes, n_features)): Training vectors.
        y (np.ndarray of shape (n_samlpes, 1)): Class labels.
        classifier (estimator object): Scikit-learn classifier.
        metrics (Metrics = Metrics.Accuracy): Metrics to measure model's performance.
        average (Average = Average.Micro): When the metrics in [Metrics.Fscore, Metrics.Precision, Metrics.Recall],
            average determines the type of averaging performed on the data.
        cv (int = 5): Specify the number of folds in a KFold.
        n_job (int = 8): Number of jobs to run in parallel.

    Return:
        k_fold_performance (np.ndarray of float of shape=cv): Array of performance of the model for
            each run of the cross validation.
    """
    # Check if input are valid.
    if not isinstance(X, np.ndarray):
        raise TypeError("X should be of type np.ndarray!")
    if not isinstance(y, np.ndarray):
        raise TypeError("y should be of type np.ndarray!")
    scoring = get_scoring_str(metrics, average)

    k_fold_performance = cross_val_score(
        classifier, X, y, scoring=scoring, cv=cv, n_jobs=n_jobs
    )
    return k_fold_performance
