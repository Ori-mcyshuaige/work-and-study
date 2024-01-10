#from .InviteFriends import *
import os, sys, json, re, time, datetime, pickle, random, requests, subprocess, tempfile
import matplotlib as mpl
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import missingno
import seaborn as sns
import scipy.stats as stats
from scipy.stats import randint as sp_randint
from skopt.space import Real, Categorical, Integer
from skopt import BayesSearchCV
from itertools import combinations

from io import StringIO
from copy import copy, deepcopy
from glob import glob
from collections import defaultdict, OrderedDict
from functools import reduce
from mpl_toolkits.mplot3d import Axes3D
from enum import Enum
from typing import List, Dict,Tuple

# Feature Selection and Encoding
from sklearn.feature_selection import RFE, RFECV,SelectFromModel, SelectKBest, f_classif,VarianceThreshold
from sklearn.svm import SVR
from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder, LabelEncoder, label_binarize, StandardScaler,MinMaxScaler,MaxAbsScaler,RobustScaler

# Machine learning 
import sklearn.ensemble as ske
from sklearn.base import BaseEstimator
from sklearn.neural_network import MLPClassifier
from sklearn.utils.validation import check_is_fitted

from sklearn.pipeline import Pipeline
from sklearn import datasets, model_selection, tree, preprocessing, metrics, linear_model
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.svm import LinearSVC,SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge, Lasso, SGDClassifier,LassoCV,RidgeCV
from sklearn.tree import DecisionTreeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import train_test_split, ShuffleSplit, cross_val_score
from sklearn.covariance import OAS
from sklearn.inspection import permutation_importance
from sklearn.datasets import make_blobs
import tensorflow as tf

# Grid and Random Search
from sklearn.model_selection import GridSearchCV,cross_val_score,RandomizedSearchCV

# Metrics
from sklearn.metrics import precision_recall_fscore_support, roc_curve, auc,precision_recall_curve,f1_score, accuracy_score, precision_score, recall_score,get_scorer
#
from utils.utils_method import *
from utils.preprocess import *
from utils.model import *
from utils.select_feature import *

# Managing Warnings 
import warnings
warnings.filterwarnings('ignore')

def standardize_data(data, mean, std):
    return (data - mean)/std

def Models(if_p_far_gt_N=0):
    # return models use to fit the data
    models = {ModelType.LR: LogisticRegression(solver='sag',max_iter=3000), 
              ModelType.KNN: KNeighborsClassifier(),
              ModelType.RF: RandomForestClassifier(min_samples_leaf=5),
              ModelType.GBDT: GradientBoostingClassifier(),
              ModelType.LDA: LinearDiscriminantAnalysis()}
    if if_p_far_gt_N:
        oa = OAS(store_precision=False, assume_centered=False)
        models['LDA'] = LinearDiscriminantAnalysis(solver="lsqr", covariance_estimator=oa)
    return(models)

def var_filter(thr,data):
    datatran = pd.DataFrame(data).apply(lambda x: x/x.mean())
    clf = VarianceThreshold(thr)
    clf.fit(datatran)
    feature = clf.get_support()
    data = pd.DataFrame(data).loc[:,feature]
    return data

def regulatization(X,y,penalty):
    if penalty == 'l1':
        clf = Lasso(alpha = .01)
    elif penalty =='l2':
        clf = Ridge(alpha = 1)
    slct = SelectFromModel(clf)
    slct.fit(X,y)
    feature = slct.get_feature_names_out()
    return feature

def Recursive_Feature_Elimination(dataset_con_enc, method_use, predclass):
    # Calculating RFE for non-discretised dataset, and graphing the Importance for each feature, per dataset
    selector1 = RFECV(method_use, step=1, cv=5, n_jobs=-1)
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

def get_best_model_and_accuracy(model,params,X,y):
    grid = GridSearchCV(model,params,error_score=0.)
    grid.fit(X,y)
    return grid.best_score_, grid.best_params_

def sel_f_model(x_train,y_train,x_test,y_test):
    selectors = {
        "LogisticRegression": LogisticRegression(solver='sag',max_iter=3000), 
        "LinearSVC": LinearSVC(), 
        "RandomForest": RandomForestClassifier(),
        "LDA": LinearDiscriminantAnalysis(),
        "LassoCV": LassoCV(fit_intercept=True, max_iter=1000),
        "RidgeCV": RidgeCV(fit_intercept=True),
        "Lasso1": Lasso(alpha = 1,fit_intercept=True, max_iter=1000),
        # "Ridge1": Ridge(alpha = 1,fit_intercept=True),
        "Lasso.1": Lasso(alpha = .1,fit_intercept=True, max_iter=1000),
        # "Ridge50": Ridge(alpha = 50,fit_intercept=True)
            }

    models = {
        ModelType.LR: LogisticRegression(),
        ModelType.KNN: KNeighborsClassifier(),
        # "Linear SVC": LinearSVC(),
        ModelType.RF: RandomForestClassifier(),
        # "Decision Tree": DecisionTreeClassifier(),
        ModelType.GBDT: GradientBoostingClassifier(),
        ModelType.LDA: LinearDiscriminantAnalysis()
        }
    features = []
    md = []
    baseauc = 0
    result = pd.DataFrame(columns = list(selectors.keys()),index = [i.value for i in list(models.keys())]+['features'])
    best_params = pd.DataFrame(columns = list(selectors.keys()),index = [i.value for i in list(models.keys())])
    x_train_var = var_filter(0.05,x_train)
    x_test_var = x_test[x_train_var.columns]
    for slctname,selector in selectors.items():
        clf = selector
        if slctname in ["LassoCV","Lasso1","Lasso.1"]:
            X = x_train
        else:
            X = x_train_var
        clf = clf.fit(X, y_train)
        model = SelectFromModel(clf, prefit=True,max_features=None)
        model.feature_names_in_ = X.columns
        features_slct = model.get_feature_names_out()
        
        # clf = RFECV(selector, step=1, cv=5, n_jobs=-1)
        # clf = clf.fit(x_train, y_train)
        # features_slct = x_train.columns[clf.support_]
        if len(features_slct) == 0:
            continue
        X_train = x_train[features_slct]
        X_test = x_test[features_slct]
        
        print('检查循环和训练集....')
        print(slctname)
        print(X_train.shape)
        model_scores, model_probs,best_param = fit_and_score(X_train, X_test, y_train, y_test,models)
        result.at['features',slctname] = ';'.join(features_slct)
        best_params[slctname] = best_param
        for modelname,prob in model_probs.items():
            fpr, tpr, threshold = metrics.roc_curve(y_test, prob)
            roc_auc = auc(fpr,tpr)
            result.at[modelname.value,slctname] = roc_auc
            if roc_auc > baseauc:
                baseauc = roc_auc
                md = models[modelname]
                mdn = modelname
                features = list(features_slct)
                model_probs_save = model_probs
            elif roc_auc == baseauc:
                if len(list(features_slct))<len(features):
                    md = models[modelname]
                    mdn = modelname
                    features = list(features_slct)
                    model_probs_save = model_probs
    return features, md, result, mdn, best_params,model_probs_save
    # clf = selectors["Logistic Regression"]
    
    # clf = clf.fit(x_train, y_train)
    # model = SelectFromModel(clf, prefit=True)
    # model.feature_names_in_ = x_train.columns
    # return x_train.columns

    

def fit_ml_algo(model_type, algo, X_train, X_test, y_train, y_test, cv=10):
    # Function that runs the requested algorithm and returns the accuracy metrics
    best_params,best_scores,cv_result = HPO(
            X_train.values,
            y_train.values,
            algo,
            ParamConfig.get_param_config(model_type, HPOAlgorithm.GridSearch),
            hpo_algorithm=HPOAlgorithm.GridSearch,
            metrics=Metrics.Accuracy,
            average=Average.Micro,
            cv=5,
            n_jobs=8,
            hpo_search_iter=100,
        )
    model = algo.set_params(**best_params)
    # best_params = []
    model = algo.fit(X_train, y_train)
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
    test_pred = model.predict(X_test)
    # CV 
    train_pred = model_selection.cross_val_predict(algo, 
                                                  X_train, 
                                                  y_train, 
                                                  cv=cv, 
                                                  n_jobs = -1)
    acc_cv = round(metrics.accuracy_score(y_train, train_pred) * 100, 2)
    return train_pred, test_pred, acc, acc_cv, probs,best_params

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
        models = {"LogisticRegression": LogisticRegression(), 
                  "KNN": KNeighborsClassifier(n_neighbors=2),
                  #"Linear SVC": LinearSVC(), 
                  "RandomForest": RandomForestClassifier(),
                  "DecisionTree": DecisionTreeClassifier(),
                  "GradientBoosting": GradientBoostingClassifier()
                  }
    # Make a list to keep model scores
    model_scores, model_probs = {}, {}
    # Loop through models
    best_params = []
    for name, model in models.items():
        # Fit the model to the data
        print(f'begin gridsearchcv-{name}')
        train_pred, test_pred, acc, acc_cv, probs,best_param = fit_ml_algo(name, model, X_train, X_test, y_train, y_test, 10)
        print(f'finish gridsearchcv-{name}')
        # Evaluate the model and append its score to model_scores
        model_scores[name] = acc
        model_probs[name]  = probs
        best_params.append(best_param)
    return model_scores, model_probs,best_params




class Biomarker_Evaluate:
    
    def __init__(self,y_true,x):
        self.y_true = y_true
        self.x = x
    
    def roc_(self):
        fpr,tpr,threshold_roc = metrics.roc_curve(self.y_true,self.x)
        auc_ = metrics.auc(fpr,tpr)
        return fpr,tpr,threshold_roc,auc_
    def prroc(self):
        precision,recall,threshold_pr = metrics.precision_recall_curve(self.y_true,self.x)
        ap = metrics.average_precision_score(self.y_true, self.x)
        return precision,recall,threshold_pr,ap
    def evaluate(self):
        fpr,tpr,threshold_roc,auc_ = self.roc_()
        precision,recall,threshold_pr,ap = self.prroc()
        
        bestpr_idx = np.argmax((2*precision*recall)/(precision+recall))
        best_threshold_pr = threshold_pr[bestpr_idx]
        best_precision = precision[bestpr_idx]
        best_recall = recall[bestpr_idx]
        pred_pr = self.x >= best_threshold_pr
        
        bestroc_idx = np.argmax(tpr-fpr)
        best_threshold_roc = threshold_roc[bestroc_idx]
        best_fpr = fpr[bestroc_idx]
        best_tpr = tpr[bestroc_idx]
        pred_roc = self.x >= best_threshold_roc
        
        evaluate = {'auc':auc_,'fpr':best_fpr,'tpr':best_tpr,'threshold_roc':best_threshold_roc,\
                'ap':ap,'precision':best_precision,'recall':best_recall,'threshold_pr':best_threshold_pr}
        return pred_roc,pred_pr,evaluate
    def confusion_matrix_(self,y_pred):
        y_true = self.y_true
        tp = sum(y_pred[[i for i,j in list(enumerate(y_true)) if j == 1]]==1)
        tn = sum(y_pred[[i for i,j in list(enumerate(y_true)) if j == 0]]==0)
        fp = sum(y_pred[[i for i,j in list(enumerate(y_true)) if j == 0]]==1)
        fn = sum(y_pred[[i for i,j in list(enumerate(y_true)) if j == 1]]==0)
        confusion_matrix=pd.DataFrame(columns = pd.MultiIndex.from_product([['Predict'],['0','1']]),
                 index = pd.MultiIndex.from_product([['True Value'],['0','1']]))
        confusion_matrix.loc['True Value','1']['Predict','1'] = tp
        confusion_matrix.loc['True Value','0']['Predict','0'] = tn
        confusion_matrix.loc['True Value','0']['Predict','1'] = fp
        confusion_matrix.loc['True Value','1']['Predict','0'] = fn
        return confusion_matrix
    def plot_roc(self,path):
        fpr,tpr,threshold_roc,auc_ = self.roc_()
        plt.figure(figsize=(4,4),dpi=600,facecolor="w")
        plt.grid(False)
        bwith = 2 #边框宽度设置为2
        ax = plt.gca()#获取边框
        ax.spines[['bottom','left','top','right']].set_color('black')  # 设置上‘脊梁’为红色
        ax.spines[['bottom','left','top','right']].set_linewidth(bwith)
        plt.title('Receiver Operating Characteristic')
        plt.plot(1-fpr, tpr, 'b', label = 'auc = %.2f\nspecificity = %.2f\nsensitivity = %.2f' % (auc_,tpr[np.argmax(tpr-fpr)],1-fpr[np.argmax(tpr-fpr)]))
        plt.legend(loc = 'lower right')
        # plt.plot([1, 0], [0, 1],'r--')
        plt.xlim([1.01, -0.01])
        plt.ylim([-0.01, 1.01])
        plt.ylabel('Specificity')
        plt.xlabel('Sensitivity')
        
        plt.savefig(path,bbox_inches='tight',dpi=600)

    def plot_prroc(self,path):
        precision,recall,threshold_pr,ap = self.prroc()
        plt.figure(figsize=(4,4),dpi=600,facecolor="w")
        plt.grid(False)
        bwith = 2 #边框宽度设置为2
        ax = plt.gca()#获取边框
        ax.spines[['bottom','left','top','right']].set_color('black')  # 设置上‘脊梁’为红色
        ax.spines[['bottom','left','top','right']].set_linewidth(bwith)
        plt.title('Precision-Recall')
        plt.plot(recall, precision, 'b', 
                 label = 'ap = %.2f\nprecision = %.2f\nrecall = %.2f'%(ap,\
                                                            precision[np.argmax((2*precision*recall)/(precision+recall))],\
                                                            recall[np.argmax((2*precision*recall)/(precision+recall))]))
        plt.legend(loc = 'lower right')
        plt.plot([0, 1], [1, 0],'r--')
        plt.xlim([-0.01, 1.01])
        plt.ylim([-0.01, 1.01])
        plt.ylabel('Recall')
        plt.xlabel('Precision')
        
        plt.savefig(path,bbox_inches='tight',dpi=600)

    def plot_confusion_matrix(self,y_pred, classes, path, title='Confusion Matrix'):
        cm = self.confusion_matrix_(y_pred)
        cm = cm.values.astype('int')
        plt.figure(figsize=(6,5),dpi=600,facecolor="w")
        np.set_printoptions(precision=2)

        # 在混淆矩阵中每格的概率值
        ind_array = np.arange(len(classes))
        x, y = np.meshgrid(ind_array, ind_array)
        for x_val, y_val in zip(x.flatten(), y.flatten()):
            c = cm[y_val][x_val]
            if c > 0.001:
                plt.text(x_val, y_val, "%d" % (c,), color='black', fontsize=15, va='center', ha='center')

        plt.imshow(cm, interpolation='nearest', cmap=plt.cm.viridis)
        plt.title(title)
        plt.colorbar()
        xlocations = np.array(range(len(classes)))
        plt.xticks(xlocations, classes, rotation=90)
        plt.yticks(xlocations, classes)
        plt.ylabel('Actual label')
        plt.xlabel('Predict label')

        # offset the tick
        tick_marks = np.array(range(len(classes))) + 0.5
        plt.gca().set_xticks(tick_marks, minor=True)
        plt.gca().set_yticks(tick_marks, minor=True)
        plt.gca().xaxis.set_ticks_position('none')
        plt.gca().yaxis.set_ticks_position('none')
        # plt.grid(True, which='minor', linestyle='-')
        plt.grid(False)
        plt.gcf().subplots_adjust(bottom=0.15)
        
        plt.savefig(path,bbox_inches='tight',dpi=600)
        


