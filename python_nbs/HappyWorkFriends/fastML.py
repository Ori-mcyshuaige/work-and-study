#from .InviteFriends import *
import os, sys, json, re, time, datetime, pickle, random, requests, subprocess, tempfile
import matplotlib as mpl
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import missingno
import seaborn as sns
import scipy.stats as scs

from io import StringIO
from copy import copy, deepcopy
from glob import glob
from collections import defaultdict, OrderedDict
from functools import reduce
from mpl_toolkits.mplot3d import Axes3D

# Feature Selection and Encoding
from sklearn.feature_selection import RFE, RFECV
from sklearn.svm import SVR
from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder, LabelEncoder, label_binarize

# Machine learning 
import sklearn.ensemble as ske
from sklearn import datasets, model_selection, tree, preprocessing, metrics, linear_model
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge, Lasso, SGDClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import train_test_split, ShuffleSplit, cross_val_score
import tensorflow as tf

# Grid and Random Search
from scipy.stats import randint as sp_randint
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV

# Metrics
from sklearn.metrics import precision_recall_fscore_support, roc_curve, auc

# Managing Warnings 
import warnings
warnings.filterwarnings('ignore')

def standardize_data(data, mean, std):
    return (data - mean)/std

def Recursive_Feature_Elimination(dataset_con_enc, predclass):
    # Calculating RFE for non-discretised dataset, and graphing the Importance for each feature, per dataset
    selector1 = RFECV(LogisticRegression(), step=1, cv=5, n_jobs=-1)
    selector1 = selector1.fit(dataset_con_enc.drop(predclass, axis=1).values, dataset_con_enc[predclass].values)
    print("Feature Ranking For Non-Discretised: %s" % selector1.ranking_)
    print("Optimal number of features : %d" % selector1.n_features_)
    # Plot number of features VS. cross-validation scores
    # plt.style.use('seaborn-whitegrid')
    # plt.figure(figsize=(20,5)) 
    # plt.xlabel("Number of features selected - Non-Discretised")
    # plt.ylabel("Cross validation score (nb of correct classifications)")
    # plt.plot(range(1, len(selector1.grid_scores_) + 1), selector1.grid_scores_);

    # Feature space could be subsetted like so:
    dataset_con_enc = dataset_con_enc[dataset_con_enc.columns[np.concatenate([selector1.support_, [True]])]]
    return(dataset_con_enc)

def fit_ml_algo(algo, X_train, X_test, y_train, y_test, cv=10):
    # Function that runs the requested algorithm and returns the accuracy metrics
    model = algo.fit(X_train, y_train)
    test_pred = model.predict(X_test)
    if (isinstance(algo, (LogisticRegression, 
                          KNeighborsClassifier, 
                          GaussianNB, 
                          DecisionTreeClassifier, 
                          RandomForestClassifier,
                          LinearSVC,
                          GradientBoostingClassifier,
                          LinearDiscriminantAnalysis))):
        probs = model.predict_proba(X_test)[:,1]
    else:
        probs = "Not Available"
    acc = round(model.score(X_test, y_test) * 100, 2) 
    # CV 
    train_pred = model_selection.cross_val_predict(algo, 
                                                  X_train, 
                                                  y_train, 
                                                  cv=cv, 
                                                  n_jobs = -1)
    acc_cv = round(metrics.accuracy_score(y_train, train_pred) * 100, 2)
    return train_pred, test_pred, acc, acc_cv, probs

# Create function to fit and score models
def fit_and_score(X_train, X_test, y_train, y_test, models=None):
    """
    Fits and evaluates given machine learning models.
    models : a dict of different Scikit-Learn machine learning models
    X_train : training data
    X_test : testing data
    y_train : labels assosciated with training data
    y_test : labels assosciated with test data
    """
    # Random seed for reproducible results
    np.random.seed(42)
    # Put models in a dictionary
    if models is None:
        models = {"Logistic Regression": LogisticRegression(), 
                  "KNN": KNeighborsClassifier(),
                  #"Linear SVC": LinearSVC(), 
                  "Random Forest": RandomForestClassifier(),
                  "Gradient Boosting Trees": GradientBoostingClassifier(),
                  "LDA": LinearDiscriminantAnalysis()}
    # Make a list to keep model scores
    model_scores, model_probs = {}, {}
    # Loop through models
    for name, model in models.items():
        # Fit the model to the data
        train_pred, test_pred, acc, acc_cv, probs = fit_ml_algo(model, X_train, X_test, y_train, y_test, 10)
        # Evaluate the model and append its score to model_scores
        model_scores[name] = acc
        model_probs[name]  = probs
    return model_scores, model_probs