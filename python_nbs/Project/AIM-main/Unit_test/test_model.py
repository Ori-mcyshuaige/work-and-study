#!/usr/bin/env python3
###############################################################
## Copyright: GeneGenieDx Corp 2022
## Author: whgu
# Date of creation: 04/04/2022
# Date of revision: 06/24/2022
#
## AIM
## Description: Unit test for functions of training and testing models.
#
###############################################################
import os
os.environ["PYTHONWARNINGS"] = "ignore"
import unittest
import sys

cur_folder = os.path.dirname(os.path.realpath(__file__))
last_folder = os.path.dirname(cur_folder)

sys.path.append(last_folder)
from sklearn import datasets

from model import *
from utils import *
from generate_feature_matrix import *


class TestModel(unittest.TestCase):
    """Test model.py"""

    def setUp(self) -> None:
        """
        Load the iris dataset (classification).
        Classes: 3,
        Samples per class: 50,
        Samples total: 150,
        Dimensionality: 4,
        Features: real, positive.
        """
        data = datasets.load_iris()
        self.X = data.data
        self.y = data.target

    def test_calculate_prediction_performance_input(self) -> None:
        """
        Test if the input of calculate_prediction_performance() are valid.
        """
        calculate_prediction_performance(y_true=[0, 1, 2, 3], y_pred=[0, 2, 1, 3])
        calculate_prediction_performance(
            y_true=np.array([0, 1, 2, 3]), y_pred=np.array([0, 2, 1, 3])
        )
        calculate_prediction_performance(
            y_true=np.array([0, 1, 2, 3]), y_pred=[0, 2, 1, 3]
        )
        calculate_prediction_performance(
            y_true=[0, 1, 2, 3], y_pred=[0, 2, 1, 3], metrics=Metrics.Fscore
        )
        calculate_prediction_performance(
            y_true=[0, 1, 2, 3],
            y_pred=[0, 2, 1, 3],
            metrics=Metrics.Fscore,
            average=Average.Micro,
        )
        with self.assertRaises(TypeError):
            calculate_prediction_performance(y_true=1, y_pred=0)
        with self.assertRaises(ValueError):
            calculate_prediction_performance(y_true=[0, 1, 2, 3], y_pred=[0, 2, 1])
        with self.assertRaises(ValueError):
            calculate_prediction_performance(y_true=[0, 1, 2], y_pred=[0, 2, 1, 3])
        with self.assertRaises(TypeError):
            calculate_prediction_performance(
                [0, 1, 2, 3], y_pred=[0, 2, 1, 3], metrics="accuracy"
            )
        with self.assertRaises(TypeError):
            calculate_prediction_performance(
                [0, 1, 2, 3], y_pred=[0, 2, 1, 3], average="micro"
            )

    def test_calculate_accuracy(self) -> None:
        """
        Test function of calculating accuracy in calculate_prediction_performance().
        """
        y_true = [0, 1, 2, 3]
        y_pred = [0, 2, 1, 3]

        accuracy = calculate_prediction_performance(y_true, y_pred)
        self.assertEqual(accuracy, 0.5)

    def test_calculate_fscore(self) -> None:
        """
        Test function of calculating fscore in calculate_prediction_performance().
        """
        y_true = [0, 1, 2, 0, 1, 2]
        y_pred = [0, 2, 1, 0, 0, 1]

        # Compute micro f1 score.
        micro_f1 = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Fscore, average=Average.Micro
        )
        self.assertEqual(micro_f1, 1 / 3)

        # Compute macro f1 score.
        macro_f1 = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Fscore, average=Average.Macro
        )
        self.assertEqual(macro_f1, 4 / 15)

        # Compute weighted f1 score.
        weighted_f1 = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Fscore, average=Average.Weighted
        )
        self.assertEqual(weighted_f1, 4 / 15)

        # Compute binary f1 score.
        y_true = [0, 1, 0, 0, 1, 1]
        y_pred = [0, 1, 1, 0, 0, 1]

        binary_f1 = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Fscore, average=Average.Binary
        )
        self.assertEqual(binary_f1, 2 / 3)

    def test_calculate_precision(self) -> None:
        """
        Test function of calculating precision in calculate_prediction_performance().
        """
        y_true = [0, 1, 2, 0, 1, 2]
        y_pred = [0, 2, 1, 0, 0, 1]

        # Compute micro precision.
        micro_precision = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Precision, average=Average.Micro
        )
        self.assertEqual(micro_precision, 1 / 3)

        # Compute macro precision.
        macro_precision = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Precision, average=Average.Macro
        )
        self.assertEqual(macro_precision, 2 / 9)

        # Compute weighted precision.
        weighted_precision = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Precision, average=Average.Weighted
        )
        self.assertEqual(weighted_precision, 2 / 9)

        # Compute binary precision.
        y_true = [0, 1, 0, 0, 1, 1]
        y_pred = [0, 1, 1, 0, 0, 1]

        binary_precision = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Precision, average=Average.Binary
        )
        self.assertEqual(binary_precision, 2 / 3)

    def test_calculate_recall(self) -> None:
        """
        Test function of calculating recall in calculate_prediction_performance().
        """
        y_true = [0, 1, 2, 0, 1, 2]
        y_pred = [0, 2, 1, 0, 0, 1]

        # Compute micro recall.
        micro_recall = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Recall, average=Average.Micro
        )
        self.assertEqual(micro_recall, 1 / 3)

        # Compute macro recall.
        macro_recall = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Recall, average=Average.Macro
        )
        self.assertEqual(macro_recall, 1 / 3)

        # Compute weighted recall.
        weighted_recall = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Recall, average=Average.Weighted
        )
        self.assertEqual(weighted_recall, 1 / 3)

        # Compute binary recall.
        y_true = [0, 1, 0, 0, 1, 1]
        y_pred = [0, 1, 1, 0, 0, 1]

        binary_recall = calculate_prediction_performance(
            y_true, y_pred, metrics=Metrics.Recall, average=Average.Binary
        )
        self.assertEqual(binary_recall, 2 / 3)

    def test_HPO_input(self) -> None:
        """
        Test if the input of HPO() are valid.
        """
        HPO(
            self.X, self.y, SVC(), ParamConfig.get_param_config(ModelType.SVM),
        )
        HPO(
            self.X,
            self.y,
            SVC(),
            ParamConfig.get_param_config(
                ModelType.SVM, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            hpo_search_iter=10,
        )
        HPO(
            self.X,
            self.y,
            SVC(),
            ParamConfig.get_param_config(
                ModelType.SVM, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            hpo_search_iter=10,
        )
        HPO(
            self.X,
            self.y,
            SVC(),
            ParamConfig.get_param_config(ModelType.SVM),
            metrics=Metrics.Fscore,
            average=Average.Macro,
        )
        HPO(self.X, self.y, SVC(), ParamConfig.get_param_config(ModelType.SVM), cv=5)
        with self.assertRaises(TypeError):
            HPO(
                self.X.tolist(),
                self.y.tolist(),
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X, self.y[:-1], SVC(), ParamConfig.get_param_config(ModelType.SVM),
            )
        with self.assertRaises(TypeError):
            HPO(
                self.X, self.y, "svm", ParamConfig.get_param_config(ModelType.SVM),
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                [
                    {
                        "penalty": ["l2"],
                        "kernel": ["linear", "poly"],
                        "gamma": ["scale"],
                    },
                ],
            )
        with self.assertRaises(TypeError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                hpo_algorithm="grid_search",
            )
        with self.assertRaises(TypeError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                metrics="accuracy",
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                cv=2.2,
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X, self.y, SVC(), ParamConfig.get_param_config(ModelType.SVM), cv=1
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                cv=1000,
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                n_jobs=2.0,
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                n_jobs=0,
            )
        with self.assertRaises(AssertionError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                hpo_search_iter=-1,
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(
                    ModelType.SVM, hpo_algotithm=HPOAlgorithm.RandomSearch
                ),
                hpo_algorithm=HPOAlgorithm.GridSearch,
            )
        with self.assertRaises(TypeError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(
                    ModelType.SVM, hpo_algotithm=HPOAlgorithm.BayesianSearch
                ),
                hpo_algorithm=HPOAlgorithm.RandomSearch,
            )

    def test_HPO_SVM(self) -> None:
        """
        Test function of hyperparameter optimization for SVM.
        We compare the cross-validation scores on the iris data before and after optimization.
        The expected result is that the score after optimization is greater than the score before optimization.
        """
        # Mean cross validataion (3 fold) accuracy with default hyper-parameters.
        score_before_HPO = cross_val_score(
            SVC(), self.X, self.y, scoring="accuracy", cv=3
        ).mean()

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though GridSearchCV.
        best_params, best_score, cv_results = HPO(
            self.X,
            self.y,
            SVC(),
            ParamConfig.get_param_config(ModelType.SVM),
            cv=3,
            n_jobs=20,
        )

        # Mean cross validataion (3 fold) accuracy after hyper-parameter optimization.
        score_after_HPO = cross_val_score(
            SVC(
                kernel=best_params["kernel"],
                C=best_params["C"],
                gamma=best_params["gamma"],
            ),
            self.X,
            self.y,
            cv=3,
        ).mean()

        self.assertGreaterEqual(score_after_HPO, score_before_HPO)
        self.assertEqual(len(cv_results.keys()), 14)

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though RandomSearchCV.
        HPO(
            self.X,
            self.y,
            SVC(),
            ParamConfig.get_param_config(
                ModelType.SVM, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though BayesianSearchCV.
        HPO(
            self.X,
            self.y,
            SVC(),
            ParamConfig.get_param_config(
                ModelType.SVM, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )

    def test_HPO_LR(self) -> None:
        """
        Test function of hyperparameter optimization for Logistic Regression.
        We compare the cross-validation scores on the iris data before and after optimization.
        The expected result is that the score after optimization is greater than the score before optimization.
        """
        # Mean cross validataion (3 fold) accuracy with default hyper-parameters.
        score_before_HPO = cross_val_score(
            LogisticRegression(), self.X, self.y, scoring="accuracy", cv=3, n_jobs=20
        ).mean()

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data through GridSearchCV.
        best_params, best_score, cv_results = HPO(
            self.X,
            self.y,
            LogisticRegression(),
            ParamConfig.get_param_config(ModelType.LR),
            cv=3,
        )

        # Mean cross validataion (3 fold) accuracy after hyper-parameter optimization.
        score_after_HPO = cross_val_score(
            LogisticRegression(
                penalty=best_params["penalty"],
                C=best_params["C"],
                multi_class=best_params["multi_class"],
                solver=best_params["solver"],
                max_iter=best_params["max_iter"],
                l1_ratio=best_params["l1_ratio"],
            ),
            self.X,
            self.y,
            cv=3,
        ).mean()

        self.assertGreaterEqual(score_after_HPO, score_before_HPO)
        self.assertEqual(len(cv_results.keys()), 17)

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though RandomSearchCV.
        HPO(
            self.X,
            self.y,
            LogisticRegression(),
            ParamConfig.get_param_config(
                ModelType.LR, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though BayesianSearchCV.
        HPO(
            self.X,
            self.y,
            LogisticRegression(),
            ParamConfig.get_param_config(
                ModelType.LR, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )


    def test_HPO_RF(self) -> None:
        """
        Test function of hyperparameter optimization for Random Forest.
        We compare the cross-validation scores on the iris data before and after optimization.
        The expected result is that the score after optimization is greater than the score before optimization.
        """
        # Mean cross validataion (3 fold) accuracy with default hyper-parameters.
        score_before_HPO = cross_val_score(
            RandomForestClassifier(random_state=123),
            self.X,
            self.y,
            scoring="accuracy",
            cv=3,
        ).mean()

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though GridSearchCV.
        best_params, best_score, cv_results = HPO(
            self.X,
            self.y,
            RandomForestClassifier(random_state=123),
            ParamConfig.get_param_config(ModelType.RF),
            cv=3,
            n_jobs=20,
        )

        # Mean cross validataion (3 fold) accuracy after hyper-parameter optimization.
        score_after_HPO = cross_val_score(
            RandomForestClassifier(
                n_estimators=best_params["n_estimators"],
                bootstrap=best_params["bootstrap"],
                criterion=best_params["criterion"],
                max_depth=best_params["max_depth"],
                max_features=best_params["max_features"],
                min_samples_leaf=best_params["min_samples_leaf"],
                random_state=123,
            ),
            self.X,
            self.y,
            cv=3,
        ).mean()

        self.assertGreaterEqual(score_after_HPO, score_before_HPO)
        self.assertEqual(len(cv_results.keys()), 17)

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though RandomSearchCV.
        HPO(
            self.X,
            self.y,
            RandomForestClassifier(),
            ParamConfig.get_param_config(
                ModelType.RF, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though BayesianSearchCV.
        HPO(
            self.X,
            self.y,
            RandomForestClassifier(),
            ParamConfig.get_param_config(
                ModelType.RF, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )

    def test_HPO_KNN(self) -> None:
        """
        Test function of hyperparameter optimization for KNN.
        We compare the cross-validation scores on the iris data before and after optimization.
        The expected result is that the score after optimization is greater than the score before optimization.
        """
        # Mean cross validataion (3 fold) accuracy with default hyper-parameters.
        score_before_HPO = cross_val_score(
            KNeighborsClassifier(), self.X, self.y, scoring="accuracy", cv=3, n_jobs=20
        ).mean()

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though GridSearchCV.
        best_params, best_score, cv_results = HPO(
            self.X,
            self.y,
            KNeighborsClassifier(),
            ParamConfig.get_param_config(ModelType.KNN),
            cv=3,
        )

        # Mean cross validataion (3 fold) accuracy after hyper-parameter optimization.
        score_after_HPO = cross_val_score(
            KNeighborsClassifier(
                n_neighbors=best_params["n_neighbors"], weights=best_params["weights"]
            ),
            self.X,
            self.y,
            cv=3,
        ).mean()

        self.assertGreaterEqual(score_after_HPO, score_before_HPO)
        self.assertEqual(len(cv_results.keys()), 13)

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though RandomSearchCV.
        HPO(
            self.X,
            self.y,
            KNeighborsClassifier(),
            ParamConfig.get_param_config(
                ModelType.KNN, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though BayesianSearchCV.
        HPO(
            self.X,
            self.y,
            KNeighborsClassifier(),
            ParamConfig.get_param_config(
                ModelType.KNN, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )

    def test_HPO_DT(self) -> None:
        """
        Test function of hyperparameter optimization for Decision Tree.
        We compare the cross-validation scores on the iris data before and after optimization.
        The expected result is that the score after optimization is greater than the score before optimization.
        """
        # Mean cross validataion (3 fold) accuracy with default hyper-parameters.
        score_before_HPO = cross_val_score(
            DecisionTreeClassifier(random_state=123),
            self.X,
            self.y,
            scoring="accuracy",
            cv=3,
        ).mean()

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though GridSearchCV.
        best_params, best_score, cv_results = HPO(
            self.X,
            self.y,
            DecisionTreeClassifier(random_state=123),
            ParamConfig.get_param_config(ModelType.DT),
            cv=3,
            n_jobs=20,
        )
        # Mean cross validataion (3 fold) accuracy after hyper-parameter optimization.
        score_after_HPO = cross_val_score(
            DecisionTreeClassifier(
                criterion=best_params["criterion"],
                splitter=best_params["splitter"],
                max_depth=best_params["max_depth"],
                max_features=best_params["max_features"],
                min_samples_leaf=best_params["min_samples_leaf"],
                random_state=123,
            ),
            self.X,
            self.y,
            cv=3,
        ).mean()
        self.assertGreaterEqual(score_after_HPO, score_before_HPO)
        self.assertEqual(len(cv_results.keys()), 16)

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though RandomSearchCV.
        HPO(
            self.X,
            self.y,
            DecisionTreeClassifier(),
            ParamConfig.get_param_config(
                ModelType.DT, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though BayesianSearchCV.
        HPO(
            self.X,
            self.y,
            DecisionTreeClassifier(),
            ParamConfig.get_param_config(
                ModelType.DT, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )

    def test_HPO_GBDT(self) -> None:
        """
        Test function of hyperparameter optimization for Gradient Boosting.
        We compare the cross-validation scores on the iris data before and after optimization.
        The expected result is that the score after optimization is greater than the score before optimization.
        """
        # Mean cross validataion (3 fold) accuracy with default hyper-parameters.
        score_before_HPO = cross_val_score(
            GradientBoostingClassifier(random_state=123),
            self.X,
            self.y,
            scoring="accuracy",
            cv=3,
        ).mean()

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though GridSearchCV.
        best_params, best_score, cv_results = HPO(
            self.X,
            self.y,
            GradientBoostingClassifier(random_state=123),
            ParamConfig.get_param_config(ModelType.GBDT),
            cv=3,
            n_jobs=20,
        )

        # Mean cross validataion (3 fold) accuracy after hyper-parameter optimization.
        score_after_HPO = cross_val_score(
            GradientBoostingClassifier(
                n_estimators=best_params["n_estimators"],
                learning_rate=best_params["learning_rate"],
                loss=best_params["loss"],
                max_depth=best_params["max_depth"],
                max_features=best_params["max_features"],
                min_samples_leaf=best_params["min_samples_leaf"],
                random_state=123,
            ),
            self.X,
            self.y,
            cv=3,
        ).mean()

        self.assertGreaterEqual(score_after_HPO, score_before_HPO)
        self.assertEqual(len(cv_results.keys()), 17)

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though RandomSearchCV.
        HPO(
            self.X,
            self.y,
            GradientBoostingClassifier(),
            ParamConfig.get_param_config(
                ModelType.GBDT, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though BayesianSearchCV.
        HPO(
            self.X,
            self.y,
            GradientBoostingClassifier(),
            ParamConfig.get_param_config(
                ModelType.GBDT, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )

    def test_HPO_MLP(self) -> None:
        """
        Test function of hyperparameter optimization for MLP Classifier.
        We compare the cross-validation scores on the iris data before and after optimization.
        The expected result is that the score after optimization is greater than the score before optimization.
        """
        # Load digits dataset.
        X, y = datasets.load_digits(return_X_y=True)

        # Mean cross validataion (3 fold) accuracy with default hyper-parameters.
        score_before_HPO = cross_val_score(
            MLPClassifier(random_state=123), X, y, scoring="accuracy", cv=3,
        ).mean()

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though GridSearchCV.
        best_params, best_performance, cv_results = HPO(
            X,
            y,
            MLPClassifier(random_state=123),
            ParamConfig.get_param_config(ModelType.MLP),
            cv=3,
            n_jobs=20,
        )

        # Mean cross validataion (3 fold) accuracy after hyper-parameter optimization.
        score_after_HPO = cross_val_score(
            MLPClassifier(
                hidden_layer_sizes=best_params["hidden_layer_sizes"],
                activation=best_params["activation"],
                alpha=best_params["alpha"],
                batch_size=best_params["batch_size"],
                learning_rate_init=best_params["learning_rate_init"],
                max_iter=best_params["max_iter"],
                random_state=123,
            ),
            X,
            y,
            cv=3,
        ).mean()
        self.assertGreaterEqual(score_after_HPO, score_before_HPO)

        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though RandomSearchCV.
        HPO(
            self.X,
            self.y,
            MLPClassifier(),
            ParamConfig.get_param_config(
                ModelType.MLP, hpo_algotithm=HPOAlgorithm.RandomSearch
            ),
            hpo_algorithm=HPOAlgorithm.RandomSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )
        # Search the parameter setting that gave the best mean cross-validated score on the hold out data though BayesianSearchCV.
        HPO(
            self.X,
            self.y,
            MLPClassifier(),
            ParamConfig.get_param_config(
                ModelType.MLP, hpo_algotithm=HPOAlgorithm.BayesianSearch
            ),
            hpo_algorithm=HPOAlgorithm.BayesianSearch,
            cv=3,
            n_jobs=20,
            hpo_search_iter=10,
        )

    def test_train_input(self) -> None:
        """
        Test if the input of train() are valid.
        """
        data, feature_maps = generate_mock_feature_matrix(200, 3, [10, 20, 30])
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            data, ["f0", "f10", "f20", "f30"]
        )
        train(train_X, train_y, test_X, test_y, ModelType.SVM)
        train(
            train_X,
            train_y,
            test_X,
            test_y,
            ModelType.SVM,
            hpo_algorithm=HPOAlgorithm.RandomSearch,
        )
        train(train_X, train_y, test_X, test_y, ModelType.SVM, metrics=Metrics.Fscore)
        train(
            train_X,
            train_y,
            test_X,
            test_y,
            ModelType.SVM,
            metrics=Metrics.Fscore,
            average=Average.Macro,
            n_jobs=20,
        )
        train(train_X, train_y, test_X, test_y, ModelType.SVM, cv=10)
        train(train_X, train_y, test_X, test_y, ModelType.SVM, display_attribute=False)
        train(train_X, train_y, test_X, test_y, ModelType.SVM, display_cv_results=True)
        with self.assertRaises(TypeError):
            train(train_X.values, train_y, test_X, test_y, ModelType.SVM)
        with self.assertRaises(TypeError):
            train(train_X.values, pd.DataFrame(train_y), test_X, test_y, ModelType.SVM)
        with self.assertRaises(TypeError):
            train(train_X, train_y, test_X.values, test_y, ModelType.SVM)
        with self.assertRaises(TypeError):
            train(train_X.values, train_y, test_X, pd.DataFrame(test_y), ModelType.SVM)
        with self.assertRaises(ValueError):
            train(train_X, train_y[1:], test_X, test_y, ModelType.SVM)
        with self.assertRaises(ValueError):
            train(train_X, train_y, test_X, test_y[:-1], ModelType.SVM)
        with self.assertRaises(TypeError):
            train(train_X, train_y, test_X, test_y, "svm")
        with self.assertRaises(TypeError):
            train(train_X, train_y, test_X, test_y, ModelType.SVM, metrics="accuracy")
        with self.assertRaises(TypeError):
            train(
                train_X,
                train_y,
                test_X,
                test_y,
                ModelType.SVM,
                metrics=Metrics.Fscore,
                average="macro",
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                cv=2.2,
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X, self.y, SVC(), ParamConfig.get_param_config(ModelType.SVM), cv=1
            )
        with self.assertRaises(ValueError):
            HPO(
                self.X,
                self.y,
                SVC(),
                ParamConfig.get_param_config(ModelType.SVM),
                cv=1000,
            )

    def test_train_output(self) -> None:
        """
        Test whether the output of train() are correct.
        """
        data, feature_maps = generate_mock_feature_matrix(200, 3, [10, 20, 30])
        train_X, train_y, test_X, test_y = prepare_train_test_data(
            data, ["f0", "f10", "f20", "f30"]
        )

        for model_type in ModelType:
            classifier, accuracy, feature_importance = train(
                train_X,
                train_y,
                test_X,
                test_y,
                model_type=model_type,
                metrics=Metrics.Accuracy,
                n_jobs=20,
            )

            check_is_fitted(classifier)
            self.assertLessEqual(accuracy, 1)
            self.assertGreaterEqual(accuracy, 0)
            self.assertEqual(len(feature_importance), 4)

    def test_get_scoring_str(self) -> None:
        """
        Test get_scoring_str().
        """
        scoring = get_scoring_str(metrics=Metrics.Accuracy, average=Average.Micro)
        self.assertEqual(scoring, "accuracy")
        scoring = get_scoring_str(metrics=Metrics.Accuracy, average=Average.Macro)
        self.assertEqual(scoring, "accuracy")
        scoring = get_scoring_str(metrics=Metrics.Accuracy, average=Average.Binary)
        self.assertEqual(scoring, "accuracy")

        scoring = get_scoring_str(metrics=Metrics.Fscore, average=Average.Binary)
        self.assertEqual(scoring, "f1")
        scoring = get_scoring_str(metrics=Metrics.Fscore, average=Average.Micro)
        self.assertEqual(scoring, "f1_micro")
        scoring = get_scoring_str(metrics=Metrics.Fscore, average=Average.Macro)
        self.assertEqual(scoring, "f1_macro")

        scoring = get_scoring_str(metrics=Metrics.Recall, average=Average.Binary)
        self.assertEqual(scoring, "recall")
        scoring = get_scoring_str(metrics=Metrics.Recall, average=Average.Micro)
        self.assertEqual(scoring, "recall_micro")
        scoring = get_scoring_str(metrics=Metrics.Recall, average=Average.Macro)
        self.assertEqual(scoring, "recall_macro")

        scoring = get_scoring_str(metrics=Metrics.Precision, average=Average.Binary)
        self.assertEqual(scoring, "precision")
        scoring = get_scoring_str(metrics=Metrics.Precision, average=Average.Micro)
        self.assertEqual(scoring, "precision_micro")
        scoring = get_scoring_str(metrics=Metrics.Precision, average=Average.Macro)
        self.assertEqual(scoring, "precision_macro")

        with self.assertRaises(TypeError):
            get_scoring_str("accuracy", average=Average.Micro)
        with self.assertRaises(TypeError):
            get_scoring_str(metrics=Metrics.Accuracy, average="micro")

    def test_k_fold_cross_validation_input(self) -> None:
        """
        Test if the input of k_fold_cross_validation() are valid.
        """
        k_fold_cross_validation(
            self.X, self.y, SVC(),
        )
        k_fold_cross_validation(
            self.X, self.y, KNeighborsClassifier(),
        )
        k_fold_cross_validation(
            self.X, self.y, SVC(), metrics=Metrics.Fscore, average=Average.Macro,
        )
        k_fold_cross_validation(self.X, self.y, SVC(), cv=2)
        k_fold_cross_validation(self.X, self.y, SVC(), cv=10)
        k_fold_cross_validation(self.X, self.y, SVC(), n_jobs=1)
        k_fold_cross_validation(self.X, self.y, SVC(), n_jobs=10)

        with self.assertRaises(TypeError):
            k_fold_cross_validation(
                self.X.tolist(), self.y.tolist(), SVC(),
            )
        with self.assertRaises(ValueError):
            k_fold_cross_validation(self.X, self.y[:-1], SVC())
        with self.assertRaises(TypeError):
            k_fold_cross_validation(self.X, self.y, "svm")
        with self.assertRaises(TypeError):
            k_fold_cross_validation(
                self.X, self.y, SVC(), metrics="accuracy",
            )
        with self.assertRaises(ValueError):
            k_fold_cross_validation(self.X, self.y, SVC(), cv=2.2)
        with self.assertRaises(ValueError):
            k_fold_cross_validation(self.X, self.y, SVC(), cv=1)
        with self.assertRaises(ValueError):
            k_fold_cross_validation(self.X, self.y, SVC(), cv=1000)
        with self.assertRaises(ValueError):
            k_fold_cross_validation(
                self.X, self.y, SVC(), n_jobs=2.0,
            )
        with self.assertRaises(ValueError):
            k_fold_cross_validation(self.X, self.y, SVC(), n_jobs=0)

    def test_k_fold_cross_validation_output(self) -> None:
        """
        Test if the output of k_fold_cross_validation() are correct.
        """
        # Test if the output is correct when cv = 2
        k_fold_performance = k_fold_cross_validation(self.X, self.y, SVC(), cv=2)
        self.assertEqual(k_fold_performance.shape, (2,))
        for performance in k_fold_performance:
            self.assertGreaterEqual(performance, 0)
            self.assertLessEqual(performance, 1)

        # Test if the output is correct when cv = 10
        k_fold_performance = k_fold_cross_validation(self.X, self.y, SVC(), cv=10)
        self.assertEqual(k_fold_performance.shape, (10,))
        for performance in k_fold_performance:
            self.assertGreaterEqual(performance, 0)
            self.assertLessEqual(performance, 1)


def main() -> None:
    unittest.main()


if __name__ == "__main__":
    main()
