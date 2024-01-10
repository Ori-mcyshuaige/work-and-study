###############################################################
## Copyright: GeneGenieDx Corp 2022
# Author: whgu
# Date of creation: 04/23/2022
# Date of revision: 06/24/2022
## AIM
## Description: Definition of Metrics, Average, ModelType, HPOAlgorithm and ParamConfig;
# Check if the input are valid;
# Display classifier's attributes and plot performance.
#
## usage:
#   import utils
###############################################################
from enum import Enum
from typing import List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes._axes import Axes
from sklearn.base import BaseEstimator
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.utils.validation import check_is_fitted
from sklearn.metrics import roc_curve, precision_recall_curve
import scipy.stats as stats
from scipy.stats import randint as sp_randint
from skopt.space import Real, Categorical, Integer


class Metrics(Enum):
    """
    Definition of metrics.
    """

    Accuracy = "accuracy"
    Fscore = "f1"
    Precision = "precision"
    Recall = "recall"


class Average(Enum):
    """
    Definition of average.
    """

    Binary = "binary"
    Micro = "micro"
    Macro = "macro"
    Weighted = "weighted"


class ModelType(Enum):
    """
    Definition of ModelType.
    """

    SVM = "SVM"
    LR = "LogisticRegression"
    RF = "RandomForest"
    KNN = "KNN"
    DT = "DecisionTree"
    GBDT = "GradientBoosting"
    MLP = "MLP"


class HPOAlgorithm(Enum):
    """
    Definition of HPOAlgorithm.
    """

    GridSearch = "grid_search"
    RandomSearch = "random_search"
    BayesianSearch = "bayesian_search"


class ParamConfig:
    """
    Define the hyper parameter configuration space.
    """

    @staticmethod
    def get_param_config(
        model_type: ModelType, hpo_algotithm: HPOAlgorithm = HPOAlgorithm.GridSearch
    ) -> List[Dict]:
        """
        Get param_config according to the ModelType and HPOAlgorithm.

        Params:
            model_type (ModelType): specified model type
            hpo_algotithm (HPOAlgorithm = HPOAlgorithm.GridSearch) : GridSearch, RandomizedSearch or BayesianSearch.
        Returns:
            param_config (List[Dict[str]]): A list of dictionary with parameter
            names as keys and lists of parameter settings to try as values.
        """
        # Check if input are valid.
        if not isinstance(model_type, ModelType):
            raise TypeError("model_type should be of type ModelType!")
        if not isinstance(hpo_algotithm, HPOAlgorithm):
            raise TypeError("hpo_algotithm should be of type HPOAlgorithm!")

        if model_type == ModelType.SVM:
            return ParamConfig.SVM_param_config(hpo_algotithm)
        elif model_type == ModelType.LR:
            return ParamConfig.LR_param_config(hpo_algotithm)
        elif model_type == ModelType.KNN:
            return ParamConfig.KNN_param_config(hpo_algotithm)
        elif model_type == ModelType.RF:
            return ParamConfig.RF_param_config(hpo_algotithm)
        elif model_type == ModelType.DT:
            return ParamConfig.DT_param_config(hpo_algotithm)
        elif model_type == ModelType.GBDT:
            return ParamConfig.GBDT_param_config(hpo_algotithm)
        elif model_type == ModelType.MLP:
            return ParamConfig.MLP_param_config(hpo_algotithm)

    @staticmethod
    def SVM_param_config(hpo_algotithm: HPOAlgorithm) -> List[Dict]:
        """
        Define the hyper parameter configuration space of SVM.

        Hyper parameters:
        kernel: SVM kernel type {"Linear Kernel", "Polynomial Kernel",
            "Gaussian Kernel", "Sigmoid Kernel"},
        C: float, Default = 1.0
           Penalty parameter C of the error term.
           The larger C is, the worse the generalization ability is, and the phenomenon of over-fitting is prone to occur;
           The smaller the C is, the better the generalization ability is, and the phenomenon of under-fitting is prone to occur.
        gamma: {"scale", "auto"} or float, default = "scale"
            Kernel coefficient for "rbf", "poly" and "sigmoid".
            If "scale", then use 1 / (n_features * X.var()) as value of gamma.
            If "auto", then use 1 / n_features.
            The larger γ is, the better the training set fits, the worse the generalization ability,
            and the phenomenon of over-fitting is prone to occur.

        Params:
            hpo_algotithm (HPOAlgorithm=HPOAlgorithm.GridSearch) : GridSearch, RandomizedSearch or BayesianSearch
                If GridSearch, then return param list;
                If RandomizedSearch, then return param distributions;
                If BayesianSearch, then return param range.
        Returns:
            param_config (List[Dict[str, List]]) :
                A list of dictionary with parameter names as keys and lists of parameter settings to try as values.
        """
        # gamma is kernel coefficient for "rbf", "poly" and "sigmoid",
        # which becomes irrelevant when kernel == "linear".
        if hpo_algotithm == HPOAlgorithm.GridSearch:
            param_config = [
                {"kernel": ["linear"], "C": [1, 10, 100, 1000], "gamma": ["scale"]},
                {
                    "kernel": ["rbf", "sigmoid", "poly"],
                    "C": [1, 10, 100, 1000],
                    "gamma": ["auto", "scale", 1, 0.1, 0.01, 0.001],
                },
            ]
        elif hpo_algotithm == HPOAlgorithm.RandomSearch:
            param_config = [
                {"kernel": ["linear"], "C": stats.uniform(1, 100), "gamma": ["scale"]},
                {
                    "kernel": ["rbf", "sigmoid", "poly"],
                    "C": stats.uniform(1, 100),
                    "gamma": ["auto", "scale", 1, 0.1, 0.01, 0.001],
                },
            ]
        elif hpo_algotithm == HPOAlgorithm.BayesianSearch:
            param_config = [
                {
                    "kernel": Categorical(["linear"]),
                    "C": Real(1, 100),
                    "gamma": Categorical(["scale"]),
                },
                {
                    "kernel": Categorical(["rbf", "sigmoid", "poly"]),
                    "C": Real(1, 100),
                    "gamma": Real(0.001, 1),
                },
            ]

        return param_config

    @staticmethod
    def LR_param_config(hpo_algotithm: HPOAlgorithm) -> List[Dict]:
        """
        Define the hyper parameter configuration space of LogisticRegression.

        Hyper parameters:
        penalty: {"l1", "l2", "elasticnet", "none"}, default = "l2"
            Specify the norm of the penalty.
            'none': no penalty is added;
            'l2': L2 penalty term and it is the default choice;
            'l1': L1 penalty term;
            'elasticnet': using both L1 and L2 penalty terms. Can only be used with saga solver
        C: float, default = 1.0
            Inverse of regularization strength.
            Smaller values specify stronger regularization.
        multi_class: {"auto", "ovr", "multinomial"}, default = "auto"
            If the option chosen is "ovr", then a binary problem is fit for each label.
            For "multinomial" the loss minimised is the multinomial loss fit across the entire probability distribution,
            even when the data is binary.
        solver: optimization algorithm -- {"newton-cg", "lbfgs", "liblinear",
            "sag", "saga"}, default = "lbfgs"
            For small datasets, 'liblinear' is a good choice, whereas 'sag' and 'saga' are faster for large ones;
            For multiclass problems, only 'newton-cg', 'sag', 'saga' and 'lbfgs' handle multinomial loss;
        max_iter: int, default = 100
            Maximum number of iterations taken for the solvers to converge.
        l1_ratio: float,
            Only used if penalty='elasticnet'.
            Setting l1_ratio=0 is equivalent to using penalty='l2', while setting l1_ratio=1 is equivalent to using penalty='l1'.
            For 0 < l1_ratio <1, the penalty is a combination of L1 and L2
        Params:
            hpo_algotithm (HPOAlgorithm=HPOAlgorithm.GridSearch) : GridSearch, RandomizedSearch or BayesianSearch
                If GridSearch, then return param list;
                If RandomizedSearch, then return param distributions;
                If BayesianSearch, then return param range.
        Returns:
            param_config (List[Dict[str]]):
                A list of dictionary with parameters names as keys and lists of
                parameter settings to try as values.
        """
        # The choice of the algorithm depends on the penalty chosen:
        # "newton-cg" - ["l2", "none"]
        # "lbfgs" - ["l2", "none"]
        # "liblinear" - ["l1", "l2"]
        # "sag" - ["l2", "none"]
        # "saga" - ["elasticnet", "l1", "l2", "none"]
        # "liblinear" is limited to one-versus-rest schemes.
        if hpo_algotithm == HPOAlgorithm.GridSearch:
            param_config = [
                {
                    "penalty": ["l2"],
                    "C": [0.01, 0.1, 1, 10, 100],
                    "multi_class": ["ovr", "multinomial"],
                    "solver": ["newton-cg", "lbfgs", "sag", "saga"],
                    "max_iter": [100, 500],
                    "l1_ratio": [None]
                },
                {
                    "penalty": ["l1"],
                    "C": [0.01, 0.1, 1, 10, 100],
                    "multi_class": ["ovr"],
                    "solver": ["liblinear", "saga"],
                    "max_iter": [100, 500],
                    "l1_ratio": [None]
                },
                {
                    "penalty": ["elasticnet"],
                    "C": [0.01, 0.1, 1, 10, 100],
                    "multi_class": ["ovr", "multinomial"],
                    "solver": ["saga"],
                    "max_iter": [100, 500],
                    "l1_ratio": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
                },
            ]
        elif hpo_algotithm == HPOAlgorithm.RandomSearch:
            param_config = [
                {
                    "penalty": ["l2"],
                    "C": stats.uniform(0.01, 100),
                    "multi_class": ["ovr", "multinomial"],
                    "solver": ["newton-cg", "lbfgs", "sag", "saga"],
                    "max_iter": [100, 500],
                    "l1_ratio": [None]
                },
                {
                    "penalty": ["l1"],
                    "C": stats.uniform(0.01, 100),
                    "multi_class": ["ovr"],
                    "solver": ["liblinear", "saga"],
                    "max_iter": [100, 500],
                    "l1_ratio": [None]
                },
                {
                    "penalty": ["elasticnet"],
                    "C": stats.uniform(0.01, 100),
                    "multi_class": ["ovr", "multinomial"],
                    "solver": ["saga"],
                    "max_iter": [100, 500],
                    "l1_ratio": stats.uniform(0.1, 0.9)
                },
            ]
        elif hpo_algotithm == HPOAlgorithm.BayesianSearch:
            param_config = [
                {
                    "penalty": Categorical(["l2"]),
                    "C": Real(0.01, 100),
                    "multi_class": Categorical(["ovr", "multinomial"]),
                    "solver": Categorical(["newton-cg", "lbfgs", "sag", "saga"]),
                    "max_iter": [100, 500],
                    "l1_ratio": [None]
                },
                {
                    "penalty": Categorical(["l1"]),
                    "C": Real(0.01, 100),
                    "multi_class": Categorical(["ovr"]),
                    "solver": Categorical(["liblinear", "saga"]),
                    "max_iter": [100, 500],
                },
                {
                    "penalty": Categorical(["elasticnet"]),
                    "C": Real(0.01, 100),
                    "multi_class": Categorical(["ovr", "multinomial"]),
                    "solver": Categorical(["saga"]),
                    "max_iter": [100, 500],
                    "l1_ratio": Real(0.1, 0.9)
                },
            ]

        return param_config

    @staticmethod
    def RF_param_config(hpo_algotithm: HPOAlgorithm) -> List[Dict]:
        """
        Define the hyper parameter configuration space of Random Forest.

        Hyper parameters:
        n_estimators: int, Default = 100
            The number of trees in the forest.
            The larger n_estimators is, the phenomenon of over-fitting is prone to occur;
            The smaller the n_estimators is, the phenomenon of under-fitting is prone to occur.
        bootstrap: bool, Default = True
            Whether bootstrap samples are used when building trees.
        criterion: {"gini", "entropy"}, default = "gini"
            The function to measure the quality of a split.
        max_depth: int, default = None
            Maximum depth of decision tree.
            The value can be set to 10-100 with large sample size or a large number of features.
        max_features: {"auto", "sqrt", "log2"}, default = "auto"
            The number of features to consider when looking for the best split.
            If None, then max_features = n_features.
        min_samples_leaf: int or float, default = 1
            The minimum number of samples required to be at a leaf node.
            If the number of leaf nodes is less than the number of samples, they will be pruned together with sibling nodes.

        Params:
            hpo_algotithm (HPOAlgorithm=HPOAlgorithm.GridSearch) : GridSearch, RandomizedSearch or BayesianSearch
                If GridSearch, then return param list;
                If RandomizedSearch, then return param distributions;
                If BayesianSearch, then return param range.
        Returns:
            param_config (List[Dict[str]]): A list of dictionary with parameters names as keys and lists
                                            of parameter settings to try as values.
        """
        if hpo_algotithm == HPOAlgorithm.GridSearch:
            param_config = [
                {
                    "n_estimators": [10, 50, 100, 150, 200],
                    "bootstrap": [True, False],
                    "criterion": ["gini", "entropy"],
                    "max_depth": [None, 10, 20, 30, 50],
                    "max_features": [None, "sqrt", "log2"],
                    "min_samples_leaf": [1, 2, 4, 8],
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.RandomSearch:
            param_config = [
                {
                    "n_estimators": sp_randint(10, 1000),
                    "bootstrap": [True, False],
                    "criterion": ["gini", "entropy"],
                    "max_depth": sp_randint(5, 50),
                    "max_features": [None, "sqrt", "log2"],
                    "min_samples_leaf": sp_randint(1, 10),
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.BayesianSearch:
            param_config = [
                {
                    "n_estimators": Integer(10, 1000),
                    "bootstrap": Categorical([True, False]),
                    "criterion": Categorical(["gini", "entropy"]),
                    "max_depth": Integer(5, 50),
                    "max_features": Categorical([None, "sqrt", "log2"]),
                    "min_samples_leaf": Integer(1, 10),
                }
            ]

        return param_config

    @staticmethod
    def KNN_param_config(hpo_algotithm: HPOAlgorithm) -> List[Dict]:
        """
        Define the hyper parameter configuration space of KNN.

        Hyper parameters:
        n_neighbors: int, default = 5
            Number of neighbors to use by default for kneighbors queries.
        weights: {"uniform", "distance"}, default='uniform'
            Weight function used in prediction.
            "uniform" : uniform weights. All points in each neighborhood are weighted equally.
            "distance" : weight points by the inverse of their distance.
            In this case, closer neighbors of a query point will have a greater influence than neighbors which are further away.

        Params:
            hpo_algotithm (HPOAlgorithm=HPOAlgorithm.GridSearch) : GridSearch, RandomizedSearch or BayesianSearch
                If GridSearch, then return param list;
                If RandomizedSearch, then return param distributions;
                If BayesianSearch, then return param range.
        Returns:
            param_config (List[Dict[str]]): A list of dictionary with parameters names as keys and lists of parameter settings to try as values.
        """
        if hpo_algotithm == HPOAlgorithm.GridSearch:
            param_config = [
                {
                    "n_neighbors": [i for i in range(1, 51)],
                    "weights": ["uniform", "distance"],
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.RandomSearch:
            param_config = [
                {"n_neighbors": range(1, 51), "weights": ["uniform", "distance"],}
            ]
        elif hpo_algotithm == HPOAlgorithm.BayesianSearch:
            param_config = [
                {
                    "n_neighbors": Integer(1, 50),
                    "weights": Categorical(["uniform", "distance"]),
                }
            ]

        return param_config

    @staticmethod
    def DT_param_config(hpo_algotithm: HPOAlgorithm) -> List[Dict]:
        """
        Define the hyper parameter configuration space of Decision Tree.

        Hyper parameters:
        criterion: {"gini", "entropy"}, default = "gini"
            The function to measure the quality of a split.
            "gini" for the Gini impurity and “entropy” for the information gain.
        splitter: {"best", "random"}, default="best"
            The strategy used to choose the split at each node.
            Supported strategies are "best" to choose the best split and "random" to choose the best random split.
        max_depth: int, default = None
            Maximum depth of decision tree.
            The value can be set to 10-100 with large sample size or a large number of features.
        max_features: {"auto", "sqrt", "log2"}, default = None
            The number of features to consider when looking for the best split.
            If None, then max_features = n_features.
        min_samples_leaf: int or float, default = 1
            The minimum number of samples required to be at a leaf node.
            If the number of leaf nodes is less than the number of samples, they will be pruned together with sibling nodes.

        Params:
            hpo_algotithm (HPOAlgorithm=HPOAlgorithm.GridSearch) : GridSearch, RandomizedSearch or BayesianSearch
                If GridSearch, then return param list;
                If RandomizedSearch, then return param distributions;
                If BayesianSearch, then return param range.
        Returns:
            param_config (List[Dict[str]]): A list of dictionary with parameters names as keys and lists of parameter settings to try as values.
        """
        if hpo_algotithm == HPOAlgorithm.GridSearch:
            param_config = [
                {
                    "criterion": ["gini", "entropy"],
                    "splitter": ["best", "random"],
                    "max_depth": [None, 10, 20, 30, 50],
                    "max_features": [None, "sqrt", "log2"],
                    "min_samples_leaf": [1, 2, 4, 8],
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.RandomSearch:
            param_config = [
                {
                    "criterion": ["gini", "entropy"],
                    "splitter": ["best", "random"],
                    "max_depth": sp_randint(5, 50),
                    "max_features": [None, "sqrt", "log2"],
                    "min_samples_leaf": sp_randint(1, 10),
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.BayesianSearch:
            param_config = [
                {
                    "criterion": Categorical(["gini", "entropy"]),
                    "splitter": Categorical(["best", "random"]),
                    "max_depth": Integer(5, 50),
                    "max_features": Categorical([None, "sqrt", "log2"]),
                    "min_samples_leaf": Integer(1, 10),
                }
            ]
        return param_config

    @staticmethod
    def GBDT_param_config(hpo_algotithm: HPOAlgorithm) -> List[Dict]:
        """
        Define the hyper parameter configuration space of Gradient Boosting Classifier.

        Hyper parameters:
        n_estimators: int, default=100,
            The number of boosting stages to perform.
            Gradient boosting is fairly robust to over-fitting so a large number usually results in better performance.
        learning_rate: float, default=0.1
            Learning rate shrinks the contribution of each tree by learning_rate.
            There is a trade-off between learning_rate and n_estimators.
        loss: {"deviance", "exponential"}, default = "deviance"
            The loss function to be optimized.
            "deviance" refers to deviance (= logistic regression) for classification with probabilistic outputs.
            For loss "exponential" gradient boosting recovers the AdaBoost algorithm
        max_depth: int, default = None
            Maximum depth of decision tree.
            The value can be set to 10-100 with large sample size or a large number of features.
        max_features: {"auto", "sqrt", "log2"}, default = None
            The number of features to consider when looking for the best split.
            If None, then max_features = n_features.
        min_samples_leaf: int or float, default = 1
            The minimum number of samples required to be at a leaf node.
            If the number of leaf nodes is less than the number of samples, they will be pruned together with sibling nodes.

        Params:
            hpo_algotithm (HPOAlgorithm=HPOAlgorithm.GridSearch) : GridSearch, RandomizedSearch or BayesianSearch
                If GridSearch, then return param list;
                If RandomizedSearch, then return param distributions;
                If BayesianSearch, then return param range.
        Returns:
            param_config (List[Dict[str]]): A list of dictionary with parameters names as keys and lists
                                            of parameter settings to try as values.
        """
        if hpo_algotithm == HPOAlgorithm.GridSearch:
            param_config = [
                {
                    "n_estimators": range(20, 121, 20),
                    "learning_rate": [0.001, 0.01, 0.1, 0.2],
                    "loss": ["deviance", "exponential"],
                    "max_depth": [None, 10, 20, 30, 50],
                    "max_features": [None, "sqrt", "log2"],
                    "min_samples_leaf": [1, 2, 4, 8],
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.RandomSearch:
            param_config = [
                {
                    "n_estimators": sp_randint(10, 1000),
                    "learning_rate": stats.uniform(0.001, 0.2),
                    "loss": ["deviance", "exponential"],
                    "max_depth": sp_randint(5, 50),
                    "max_features": [None, "sqrt", "log2"],
                    "min_samples_leaf": sp_randint(1, 10),
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.BayesianSearch:
            param_config = [
                {
                    "n_estimators": Integer(10, 1000),
                    "learning_rate": Real(0.001, 0.2),
                    "loss": Categorical(["deviance"]),
                    "max_depth": Integer(5, 50),
                    "max_features": Categorical([None, "sqrt", "log2"]),
                    "min_samples_leaf": Integer(1, 10),
                }
            ]

        return param_config

    @staticmethod
    def MLP_param_config(hpo_algotithm: HPOAlgorithm) -> List[Dict]:
        """
        Define the hyper parameter configuration space of multi-layer perceptron Classifier.

        Hyper parameters:
        hidden_layer_sizes: tuple, length = n_layers - 2, default=(100,)
            The ith element represents the number of neurons in the ith hidden layer.
        activation: {"logistic", "tanh", "relu"}, default="relu"
            Activation function for the hidden layer.
        alpha: float, default=0.0001
            Strength of the L2 regularization term.
        batch_size: int,
            Size of minibatches for stochastic optimizers.
            When set to "auto", batch_size=min(200, n_samples).
        learning_rate_init: float, default=0.001
            The initial learning rate used. It controls the step-size in updating the weights.
        max_iter: int, default = 200
            Number of epochs.

        Params:
            hpo_algotithm (HPOAlgorithm) :  GridSearch or RandomizedSearch
                If GridSearchCV, then return param list;
                If RandomizedSearchCV, then return param distributions.
        Return:
            param_config (List[Dict[str]]): A list of dictionary with parameters names as keys and lists
                                            of parameter settings to try as values.
        """
        if hpo_algotithm == HPOAlgorithm.GridSearch:
            param_config = [
                {
                    "hidden_layer_sizes": [
                        (50,),
                        (100,),
                        (200,),
                        (100, 50),
                        (200, 100),
                        (200, 100, 50),
                    ],
                    "activation": ["tanh", "relu", "logistic"],
                    "alpha": [0.00001, 0.0001, 0.001],
                    "batch_size": [16, 32, 64],
                    "learning_rate_init": [0.0001, 0.0005, 0.001, 0.005, 0.01],
                    "max_iter": [100, 200],
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.RandomSearch:
            param_config = [
                {
                    "hidden_layer_sizes": [
                        (50,),
                        (100,),
                        (200,),
                        (100, 50),
                        (200, 100),
                        (200, 100, 50),
                    ],
                    "activation": ["tanh", "relu", "logistic"],
                    "alpha": stats.uniform(0.00001, 0.001),
                    "batch_size": [16, 32, 64],
                    "learning_rate_init": stats.uniform(0.0001, 0.01),
                    "max_iter": [100, 200],
                }
            ]
        elif hpo_algotithm == HPOAlgorithm.BayesianSearch:
            param_config = [
                {
                    "hidden_layer_sizes": list(range(50, 201, 10)),
                    "activation": Categorical(["tanh", "relu", "logistic"]),
                    "alpha": Real(0.00001, 0.001),
                    "batch_size": Categorical([16, 32, 64]),
                    "learning_rate_init": Real(0.0001, 0.01),
                    "max_iter": Categorical([100, 200]),
                }
            ]

        return param_config


def check_X_y(X, y) -> None:
    """
    This function is to check if input `X` is of type pd.DataFrame and `y` is of type np.ndarray.

    Params:
        X: Input training vectors.
        y: Class labels.
    """
    if not isinstance(X, pd.DataFrame):
        raise TypeError("X should be of type pd.DataFrame!")
    if not isinstance(y, np.ndarray):
        raise TypeError("y should be of type np.ndarray!")


def get_scoring_str(metrics, average) -> str:
    """
    This function is to check if input `metrics` is of type Metrics and `average` is type of Average.
    Then return the scoring parameters for cross validation according to the `metrics` and `average`.

    Params:
        metrics: Metrics to measure model's performance.
        average: When the metrics in [Metrics.Fscore, Metrics.Precision, Metrics.Recall],
                average determines the type of averaging performed on the data.
    Return:
        scoring (str): Define the strategy to evaluate the performance of the cross-validation.
    """
    if not isinstance(metrics, Metrics):
        raise TypeError("metrics should be of type Metrics!")
    if not isinstance(average, Average):
        raise TypeError("average should be of type Average!")

    if metrics == Metrics.Accuracy or average == Average.Binary:
        scoring = metrics.value
    else:
        scoring = metrics.value + "_" + average.value

    return scoring


def display_classifier_attr(classifier: BaseEstimator, best_params: Dict):
    """
    Display classifier's parameters and attributes.
    SVM: Support vectors.
    LogisticRegression: Coefficient and bias.
    RandomForestClassifier: The impurity-based feature importances.
    KNN: Distance metric.
    DecisionTreeClassifier: The impurity-based feature importances.
    GBDT: Feature importances and train score.

    Params:
        classifier (estimator object): Scikit-learn classifier.
        best_params (Dict[str, str]): Parameters with the best performance via CV on training set.
    """
    if not isinstance(classifier, BaseEstimator):
        raise TypeError("classifier should be of type sklearn.BaseEstimator!")
    if not isinstance(best_params, Dict):
        raise TypeError("best_params should be of type Dict!")
    for param in best_params.keys():
        if param not in classifier.get_params():
            raise ValueError(f"classifier doesn't contrain parameter {param}!")
    check_is_fitted(classifier)

    print("Best parameters selected:", best_params)

    if isinstance(classifier, SVC):
        print("Support vectors:", classifier.__dict__["support_vectors_"])
    elif isinstance(classifier, LogisticRegression):
        print("Coefficient:", classifier.__dict__["coef_"])
        print("Bias:", classifier.__dict__["intercept_"])
    elif isinstance(classifier, RandomForestClassifier):
        print("Feature importances:", classifier.feature_importances_.tolist())
    elif isinstance(classifier, KNeighborsClassifier):
        print("Distance metric:", classifier.__dict__["effective_metric_"])
    elif isinstance(classifier, DecisionTreeClassifier):
        print("Feature importances:", classifier.feature_importances_.tolist())
    elif isinstance(classifier, GradientBoostingClassifier):
        print("Feature importances:", classifier.feature_importances_.tolist())
        print("Train score:", classifier.__dict__["train_score_"])
    elif isinstance(classifier, MLPClassifier):
        fig, ax = plt.subplots()
        ax.plot(classifier.__dict__["loss_curve_"])
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Loss")
        ax.set_title("Loss curve")


def plot_performance(
    classifier: BaseEstimator,
    test_X: np.ndarray,
    test_y: np.ndarray,
    roc: bool = True,
    pos_label: int = 1,
    ax: Axes = None,
):
    """
    Plot Receiver operating characteristic (ROC) curve and Precision Recall Curve.
    For binary classification:
        a. ROC curve.
        b. P-R curve.
    For multi class classification:
        a. ROC curve with given label of the positive class.
        b. P-R curve with given label of the positive class.

    Params:
        classifier (BaseEstimator): Trained classifier.
        test_X (np.ndarray of shape (n_samlpes, n_features)): Features of test set.
        test_y (np.ndarray of shape (n_samlpes, 1)): Labels of test set.
        roc (bool = True):
            True: Plot Receiver operating characteristic (ROC) curve.
            False: Plot Precision Recall Curve.
        pos_label (int = 1): The class considered as the positive class.
        average (Average = Average.Binary): average determines the type of averaging performed on the data.
                If Average.Binary, then only report results for the class == pos_label.
                If Average.Micro, then calculate metrics globally by counting the total true positives,
                false negatives and false positives.
                If Average.Macro, then calculate metrics for each label, and find their mean.
        ax (Axes = None): Axes object to plot on. If `None`, a new figure and axes is created.
    """
    # Check if input are valid.
    if not isinstance(classifier, BaseEstimator):
        raise TypeError("Classifier must be of type BaseEstimator!")
    if not isinstance(test_X, np.ndarray):
        raise TypeError("test_X should be of type np.ndarray!")
    if not isinstance(test_y, np.ndarray):
        raise TypeError("test_y should be of type np.ndarray!")
    if pos_label not in test_y:
        raise ValueError("Positive label doesn't exist in test_y!")
    if ax and not isinstance(ax, Axes):
        raise TypeError("Input ax must be of type Axes!")

    predict_y = classifier.predict_proba(test_X)
    pos_class_idx = classifier.classes_.tolist().index(pos_label)
    if roc:
        fpr, tpr, thresholds = roc_curve(
            y_true=test_y, y_score=predict_y[:, pos_class_idx], pos_label=pos_label
        )
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(fpr, tpr)
        ax.set_xlabel("FPR")
        ax.set_ylabel("TPR")
        ax.set_title("ROC curve")

    else:
        precision, recall, thresholds = precision_recall_curve(
            y_true=test_y, probas_pred=predict_y[:, pos_class_idx], pos_label=pos_label
        )

        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(recall, precision)
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_title("P-R curve")
