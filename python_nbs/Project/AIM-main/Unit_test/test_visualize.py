#!/usr/bin/env python3
###############################################################
## Copyright: GeneGenieDx Corp 2022
## Author: whgu
# Date of creation: 05/30/2022
#
## AIM
## Description: Unit test for displaying classifiers' parameters and attributes,
# plotting classifiers' performance.
#
###############################################################
import os
os.environ["PYTHONWARNINGS"] = "ignore"
import unittest
from unittest.mock import patch
import sys

cur_folder = os.path.dirname(os.path.realpath(__file__))
last_folder = os.path.dirname(cur_folder)

sys.path.append(last_folder)
from sklearn import datasets
from sklearn.exceptions import NotFittedError

from model import *


class TestDisplay(unittest.TestCase):
    """Test utils.display_classifier_attr() and utils.plot_performance()"""

    def setUp(self) -> None:
        """
        Load the iris dataset(multi class) and breast cancer dataset (binary).
        iris data:
        Classes: 3,
        Samples per class: 50,
        Samples total: 150,
        Dimensionality: 4,
        Features: real, positive.

        breast cancer dataset:
        Classes: 2,
        Samples per class: 212(M),357(B)
        Samples total: 569
        Dimensionality: 30,
        Features: real, positive.
        """
        # Load iris dataset.
        iris_data = datasets.load_iris()
        self.X_iris = iris_data.data
        self.y_iris = iris_data.target

        # Load breast cancer dataset.
        breast_cancer_data = datasets.load_breast_cancer()
        self.X_breast_cancer = breast_cancer_data.data
        self.y_breast_cancer = breast_cancer_data.target

    def test_display_input(self) -> None:
        """
        Test whether the input of display_classifier_attr() is valid.
        """
        best_params, best_score, cv_results = HPO(
            self.X_iris, self.y_iris, SVC(), ParamConfig.get_param_config(ModelType.SVM)
        )
        classifier = SVC(
            kernel=best_params["kernel"],
            C=best_params["C"],
            gamma=best_params["gamma"],
        )
        classifier.fit(self.X_iris, self.y_iris)
        display_classifier_attr(classifier, best_params)

        with self.assertRaises(TypeError):
            display_classifier_attr("SVM", best_params)
        with self.assertRaises(TypeError):
            display_classifier_attr(classifier, list(best_params))
        with self.assertRaises(ValueError):
            display_classifier_attr(classifier, {"penalty": "l2"})

        # Test classifier without fitting data.
        with self.assertRaises(NotFittedError):
            display_classifier_attr(SVC(), {})

    def test_display_print(self) -> None:
        """
        Test whether the print of display_classifier_attr() is correct.
        """
        # Test for SVM.
        best_params, best_score, cv_results = HPO(
            self.X_iris, self.y_iris, SVC(), ParamConfig.get_param_config(ModelType.SVM)
        )
        classifier = SVC(
            kernel=best_params["kernel"],
            C=best_params["C"],
            gamma=best_params["gamma"],
        )
        classifier.fit(self.X_iris, self.y_iris)

        with patch("utils.print") as mocked_print:
            display_classifier_attr(classifier, best_params)
            mocked_print.assert_any_call("Best parameters selected:", best_params)
            mocked_print.assert_called_with(
                "Support vectors:", classifier.__dict__["support_vectors_"],
            )

        # Test for LogisticRegression.
        best_params, best_score, cv_results = HPO(
            self.X_iris,
            self.y_iris,
            LogisticRegression(),
            ParamConfig.get_param_config(ModelType.LR),
        )
        classifier = LogisticRegression(
            penalty=best_params["penalty"],
            C=best_params["C"],
            multi_class=best_params["multi_class"],
            solver=best_params["solver"],
            max_iter=best_params["max_iter"],
            l1_ratio=best_params["l1_ratio"],
        )
        classifier.fit(self.X_iris, self.y_iris)

        with patch("utils.print") as mocked_print:
            display_classifier_attr(classifier, best_params)
            mocked_print.assert_any_call("Best parameters selected:", best_params)
            mocked_print.assert_any_call(
                "Coefficient:", classifier.__dict__["coef_"],
            )
            mocked_print.assert_called_with(
                "Bias:", classifier.__dict__["intercept_"],
            )

        # Test for RandomForest.
        best_params, best_score, cv_results = HPO(
            self.X_iris,
            self.y_iris,
            RandomForestClassifier(),
            ParamConfig.get_param_config(ModelType.RF),
            n_jobs=20,
        )
        classifier = RandomForestClassifier(
            n_estimators=best_params["n_estimators"],
            bootstrap=best_params["bootstrap"],
            criterion=best_params["criterion"],
            max_depth=best_params["max_depth"],
            max_features=best_params["max_features"],
            min_samples_leaf=best_params["min_samples_leaf"],
        )
        classifier.fit(self.X_iris, self.y_iris)

        with patch("utils.print") as mocked_print:
            display_classifier_attr(classifier, best_params)
            mocked_print.assert_any_call("Best parameters selected:", best_params)
            mocked_print.assert_called_with(
                "Feature importances:", classifier.feature_importances_.tolist(),
            )

        # Test for KNN.
        best_params, best_score, cv_results = HPO(
            self.X_iris,
            self.y_iris,
            KNeighborsClassifier(),
            ParamConfig.get_param_config(ModelType.KNN),
            n_jobs=20,
        )
        classifier = KNeighborsClassifier(n_neighbors=best_params["n_neighbors"])
        classifier.fit(self.X_iris, self.y_iris)

        with patch("utils.print") as mocked_print:
            display_classifier_attr(classifier, best_params)
            mocked_print.assert_any_call("Best parameters selected:", best_params)
            mocked_print.assert_called_with(
                "Distance metric:", classifier.__dict__["effective_metric_"]
            )

        # Test for DecisionTree.
        best_params, best_score, cv_results = HPO(
            self.X_iris,
            self.y_iris,
            DecisionTreeClassifier(),
            ParamConfig.get_param_config(ModelType.DT),
            n_jobs=20,
        )
        classifier = DecisionTreeClassifier(
            criterion=best_params["criterion"],
            splitter=best_params["splitter"],
            max_depth=best_params["max_depth"],
            max_features=best_params["max_features"],
            min_samples_leaf=best_params["min_samples_leaf"],
        )
        classifier.fit(self.X_iris, self.y_iris)

        with patch("utils.print") as mocked_print:
            display_classifier_attr(classifier, best_params)
            mocked_print.assert_any_call("Best parameters selected:", best_params)
            mocked_print.assert_called_with(
                "Feature importances:", classifier.feature_importances_.tolist(),
            )

        # Test for GBDT.
        best_params, best_score, cv_results = HPO(
            self.X_iris,
            self.y_iris,
            GradientBoostingClassifier(),
            ParamConfig.get_param_config(ModelType.GBDT),
            n_jobs=20,
        )
        classifier = GradientBoostingClassifier(
            n_estimators=best_params["n_estimators"],
            learning_rate=best_params["learning_rate"],
            loss=best_params["loss"],
            max_depth=best_params["max_depth"],
            max_features=best_params["max_features"],
            min_samples_leaf=best_params["min_samples_leaf"],
        )
        classifier.fit(self.X_iris, self.y_iris)

        with patch("utils.print") as mocked_print:
            display_classifier_attr(classifier, best_params)
            mocked_print.assert_any_call("Best parameters selected:", best_params)
            mocked_print.assert_any_call(
                "Feature importances:", classifier.feature_importances_.tolist(),
            )
            mocked_print.assert_called_with(
                "Train score:", classifier.__dict__["train_score_"]
            )

        # Test for MLP
        best_params, best_performance, cv_results = HPO(
            self.X_iris,
            self.y_iris,
            MLPClassifier(),
            ParamConfig.get_param_config(ModelType.MLP),
            n_jobs=20,
        )
        classifier = MLPClassifier(
            hidden_layer_sizes=best_params["hidden_layer_sizes"],
            activation=best_params["activation"],
            alpha=best_params["alpha"],
            batch_size=best_params["batch_size"],
            learning_rate_init=best_params["learning_rate_init"],
            max_iter=best_params["max_iter"],
        )

        classifier.fit(
            self.X_iris, self.y_iris,
        )
        with patch("utils.Axes.plot") as mocked_plot:
            display_classifier_attr(classifier, best_params)
            mocked_plot.assert_any_call(classifier.__dict__["loss_curve_"])

    def test_plot_input(self) -> None:
        """
        Test whether the input of plot_performance() is valid.
        """
        SVM_iris = SVC(probability=True)  # Classifier for iris dataset.
        SVM_breast_cancer = SVC(
            probability=True
        )  # Classifier for breast cancer dataset.

        # Test whether the input classifier is valid.
        with self.assertRaises(NotFittedError):
            plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=1)
        with self.assertRaises(TypeError):
            plot_performance("SVM", self.X_iris, self.y_iris, pos_label=1)
        SVM_iris.fit(self.X_iris, self.y_iris)
        SVM_breast_cancer.fit(self.X_breast_cancer, self.y_breast_cancer)
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=1)
        plot_performance(
            SVM_breast_cancer, self.X_breast_cancer, self.y_breast_cancer, pos_label=1,
        )
        KNN_iris = KNeighborsClassifier()
        KNN_iris.fit(self.X_iris, self.y_iris)
        plot_performance(KNN_iris, self.X_iris, self.y_iris, pos_label=1, roc=True)
        LR_iris = LogisticRegression()
        LR_iris.fit(self.X_iris, self.y_iris)
        plot_performance(LR_iris, self.X_iris, self.y_iris, pos_label=1, roc=True)
        DT_iris = DecisionTreeClassifier()
        DT_iris.fit(self.X_iris, self.y_iris)
        plot_performance(DT_iris, self.X_iris, self.y_iris, pos_label=1, roc=True)

        with self.assertRaises(AttributeError):
            classifier = SVC()
            classifier.fit(self.X_iris, self.y_iris)
            plot_performance(classifier, self.X_iris, self.y_iris, pos_label=1)

        # Test whether the input test_X and test_y are valid.
        with self.assertRaises(TypeError):
            plot_performance(SVM_iris, self.X_iris.tolist(), self.y_iris, pos_label=1)
        with self.assertRaises(TypeError):
            plot_performance(SVM_iris, self.X_iris, self.y_iris.tolist(), pos_label=1)
        with self.assertRaises(ValueError):
            plot_performance(SVM_iris, self.X_iris, self.y_iris[:-1], pos_label=1)
        with self.assertRaises(ValueError):
            plot_performance(SVM_breast_cancer, self.X_iris, self.y_iris)

        # Test whether the input pos_label is valid.
        plot_performance(SVM_breast_cancer, self.X_breast_cancer, self.y_breast_cancer)
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=0)
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=1)
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=2)
        with self.assertRaises(ValueError):
            plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=3)
        with self.assertRaises(ValueError):
            plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label="2")

        # Test whether the input roc is valid.
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=1, roc=True)
        plot_performance(
            SVM_breast_cancer,
            self.X_breast_cancer,
            self.y_breast_cancer,
            pos_label=1,
            roc=True,
        )
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=1, roc=False)
        plot_performance(
            SVM_breast_cancer,
            self.X_breast_cancer,
            self.y_breast_cancer,
            pos_label=1,
            roc=True,
        )

        # Test whether the input axis is valid.
        fig, ax = plt.subplots()
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=0, ax=ax)
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=1, ax=ax)
        plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=2, ax=ax)
        with self.assertRaises(TypeError):
            plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=0, ax="ax0")

    def test_plot_output(self) -> None:
        """
        Test whether the plotting of plot_performance() is valid.
        """
        SVM_iris = SVC(probability=True)  # SVM for iris dataset.
        SVM_breast_cancer = SVC(probability=True)  # SVM for breast cancer dataset.
        SVM_iris.fit(self.X_iris, self.y_iris)
        SVM_breast_cancer.fit(self.X_breast_cancer, self.y_breast_cancer)

        # Plot ROC curve for SVM on iris datset(multi class).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = SVM_iris.predict_proba(self.X_iris)
            # Specify class 0 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_iris, y_score=predict_y[:, 0], pos_label=0
            )
            plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=0)
            assert np.array_equal(mocked_plt.call_args_list[0][0], (fpr, tpr))
            # Specify class 1 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_iris, y_score=predict_y[:, 1], pos_label=1
            )
            plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=1)
            assert np.array_equal(mocked_plt.call_args_list[1][0], (fpr, tpr))
            # Specify class 2 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_iris, y_score=predict_y[:, 2], pos_label=2
            )
            plot_performance(SVM_iris, self.X_iris, self.y_iris, pos_label=2)
            assert np.array_equal(mocked_plt.call_args_list[2][0], (fpr, tpr))

        # Plot ROC curve for SVM on breast cancer datset(binary).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = SVM_breast_cancer.predict_proba(self.X_breast_cancer)
            # Specify class 0 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_breast_cancer, y_score=predict_y[:, 0], pos_label=0
            )
            plot_performance(
                SVM_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                pos_label=0,
            )
            assert np.array_equal(mocked_plt.call_args_list[0][0], (fpr, tpr))
            # Specify class 1 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_breast_cancer, y_score=predict_y[:, 1], pos_label=1
            )
            plot_performance(
                SVM_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                pos_label=1,
            )
            assert np.array_equal(mocked_plt.call_args_list[1][0], (fpr, tpr))

        # Plot P-R curve for SVM on iris datset(multi class).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = SVM_iris.predict_proba(self.X_iris)
            # Specify class 0 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_iris, probas_pred=predict_y[:, 0], pos_label=0
            )
            plot_performance(SVM_iris, self.X_iris, self.y_iris, roc=False, pos_label=0)
            assert np.array_equal(mocked_plt.call_args_list[0][0], (recall, precision))
            # Specify class 1 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_iris, probas_pred=predict_y[:, 1], pos_label=1
            )
            plot_performance(SVM_iris, self.X_iris, self.y_iris, roc=False, pos_label=1)
            assert np.array_equal(mocked_plt.call_args_list[1][0], (recall, precision))
            # Specify class 2 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_iris, probas_pred=predict_y[:, 2], pos_label=2
            )
            plot_performance(SVM_iris, self.X_iris, self.y_iris, roc=False, pos_label=2)
            assert np.array_equal(mocked_plt.call_args_list[2][0], (recall, precision))

        # Plot P-R curve for SVM on breast cancer datset(binary).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = SVM_breast_cancer.predict_proba(self.X_breast_cancer)
            # Specify class 0 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_breast_cancer, probas_pred=predict_y[:, 0], pos_label=0
            )
            plot_performance(
                SVM_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                roc=False,
                pos_label=0,
            )
            assert np.array_equal(mocked_plt.call_args_list[0][0], (recall, precision))
            # Specify class 1 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_breast_cancer, probas_pred=predict_y[:, 1], pos_label=1
            )
            plot_performance(
                SVM_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                roc=False,
                pos_label=1,
            )
            assert np.array_equal(mocked_plt.call_args_list[1][0], (recall, precision))

        RF_iris = RandomForestClassifier()  # random forest for iris dataset.
        RF_breast_cancer = (
            RandomForestClassifier()
        )  # random forest for breast cancer dataset.
        RF_iris.fit(self.X_iris, self.y_iris)
        RF_breast_cancer.fit(self.X_breast_cancer, self.y_breast_cancer)

        # Plot ROC curve for RF on iris datset(multi class).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = RF_iris.predict_proba(self.X_iris)
            # Specify class 0 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_iris, y_score=predict_y[:, 0], pos_label=0
            )
            plot_performance(RF_iris, self.X_iris, self.y_iris, pos_label=0)
            assert np.array_equal(mocked_plt.call_args_list[0][0], (fpr, tpr))
            # Specify class 1 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_iris, y_score=predict_y[:, 1], pos_label=1
            )
            plot_performance(RF_iris, self.X_iris, self.y_iris, pos_label=1)
            assert np.array_equal(mocked_plt.call_args_list[1][0], (fpr, tpr))
            # Specify class 2 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_iris, y_score=predict_y[:, 2], pos_label=2
            )
            plot_performance(RF_iris, self.X_iris, self.y_iris, pos_label=2)
            assert np.array_equal(mocked_plt.call_args_list[2][0], (fpr, tpr))

        # Plot ROC curve for RF on breast cancer datset(binary).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = RF_breast_cancer.predict_proba(self.X_breast_cancer)
            # Specify class 0 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_breast_cancer, y_score=predict_y[:, 0], pos_label=0
            )
            plot_performance(
                RF_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                pos_label=0,
            )
            assert np.array_equal(mocked_plt.call_args_list[0][0], (fpr, tpr))
            # Specify class 1 as pos label.
            fpr, tpr, thresholds = roc_curve(
                y_true=self.y_breast_cancer, y_score=predict_y[:, 1], pos_label=1
            )
            plot_performance(
                RF_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                pos_label=1,
            )
            assert np.array_equal(mocked_plt.call_args_list[1][0], (fpr, tpr))

        # Plot P-R curve for RF on iris datset(multi class).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = RF_iris.predict_proba(self.X_iris)
            # Specify class 0 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_iris, probas_pred=predict_y[:, 0], pos_label=0
            )
            plot_performance(RF_iris, self.X_iris, self.y_iris, roc=False, pos_label=0)
            assert np.array_equal(mocked_plt.call_args_list[0][0], (recall, precision))
            # Specify class 1 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_iris, probas_pred=predict_y[:, 1], pos_label=1
            )
            plot_performance(RF_iris, self.X_iris, self.y_iris, roc=False, pos_label=1)
            assert np.array_equal(mocked_plt.call_args_list[1][0], (recall, precision))
            # Specify class 2 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_iris, probas_pred=predict_y[:, 2], pos_label=2
            )
            plot_performance(RF_iris, self.X_iris, self.y_iris, roc=False, pos_label=2)
            assert np.array_equal(mocked_plt.call_args_list[2][0], (recall, precision))

        # Plot P-R curve for RF on breast cancer datset(binary).
        with patch("utils.Axes.plot") as mocked_plt:
            predict_y = RF_breast_cancer.predict_proba(self.X_breast_cancer)
            # Specify class 0 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_breast_cancer, probas_pred=predict_y[:, 0], pos_label=0
            )
            plot_performance(
                RF_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                roc=False,
                pos_label=0,
            )
            assert np.array_equal(mocked_plt.call_args_list[0][0], (recall, precision))
            # Specify class 1 as pos label.
            precision, recall, thresholds = precision_recall_curve(
                y_true=self.y_breast_cancer, probas_pred=predict_y[:, 1], pos_label=1
            )
            plot_performance(
                RF_breast_cancer,
                self.X_breast_cancer,
                self.y_breast_cancer,
                roc=False,
                pos_label=1,
            )
            assert np.array_equal(mocked_plt.call_args_list[1][0], (recall, precision))


def main() -> None:
    unittest.main()


if __name__ == "__main__":
    main()
