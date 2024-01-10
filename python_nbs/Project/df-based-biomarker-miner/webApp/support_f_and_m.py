import os,uuid,sys,base64,math
import pandas as pd
import numpy as np
import GEOparse
import zipfile
import re
import itertools as it


from flask import Flask,render_template,flash,redirect,url_for,\
request,session,jsonify,send_from_directory,send_file,make_response
from flask_ckeditor import CKEditor,upload_success,upload_fail
from flask_dropzone import Dropzone
from flask_wtf.csrf import validate_csrf
from wtforms import ValidationError
from markupsafe import escape
from werkzeug.utils import secure_filename
from werkzeug.datastructures import FileStorage
from matplotlib.figure import Figure
from Bio import Entrez
from scipy import stats

from utils.PLOTS import *
from utils.fastML import *
from utils.utils_method import *
from utils.preprocess import *
from utils.model import *
from utils.select_feature import *


from forms import *
from io import BytesIO

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r,pandas2ri
from rpy2.robjects.conversion import localconverter




def gene_data_preprocess(df_raw, predclass, test, con, imputeNA=1, filt_outlier=1):
    ## 检查特征名称是否合规，包括空值、同义名、重名
    # df_raw=df_raw.T
    # df_raw.iloc[:,:-1][df_raw.iloc[:,:-1].astype('float')<0] = 0.000000001
    # df_raw=df_raw.T
    df_raw=df_raw[~pd.isna(df_raw.index)]
    df_raw=df_raw[df_raw.index.str.contains('/')==False]
    ro.globalenv['rdf_raw']=pandas2ri.py2rpy(df_raw.drop(predclass).astype('float'))
    ro.globalenv['genes']=pandas2ri.py2rpy(df_raw.drop(predclass).index)
    rdf_raw=ro.r('''df_raw = as.data.frame(avereps(rdf_raw,ID=genes))''')
    df_raw=pd.concat([rdf_raw,pd.DataFrame(df_raw.loc[predclass,:]).T])
    # df_raw=pd.concat([df_raw.drop(predclass).groupby(level=0).mean(),pd.DataFrame(df_raw.loc[predclass,:]).T])
    df_raw=df_raw.T        
    ## 按照每个观察值缺失或为0的特征数，按百分位法去离群值
    if filt_outlier == 1:
        # within 2 standerd deviation
        # df_dat[np.abs(df_dat.A - df_dat.A.mean()) <= (2 * df_dat.A.std())]
        countsNa=~((df_raw.drop([predclass],axis=1)==0) ==\
               (pd.isnull(df_raw.drop([predclass],axis=1)))).sum(axis=1)
        iqr=np.percentile(countsNa,75,axis=0)-np.percentile(countsNa,25,axis=0)
        lowerlimt,upperlimt = np.percentile(countsNa,25,axis=0)-1.5*iqr,np.percentile(countsNa,75,axis=0)+1.5*iqr
        df_raw=df_raw[[x<=upperlimt for x in countsNa]]
    ## KNN法补值
    if imputeNA == 1:
        imp_knn = KNNImputer(n_neighbors=5)
        df_raw = pd.concat([pd.DataFrame(imp_knn.fit_transform(df_raw.drop([predclass],axis=1)),\
                  columns = df_raw.columns.drop(predclass),index=df_raw.index),df_raw[predclass]],axis=1)
    ## 检查数据是否已log2转换
    qx = np.percentile(df_raw.iloc[:,:-1],[0,25,50,75,99,100])
    LogC=((qx[4]>100)|(((qx[5]-qx[0])>50) & (qx[1]>0)))
    if LogC:
        df1 = df_raw.drop(predclass,axis=1).astype('float')
        df1[df1<0] = 0.000000001
        df1[predclass] = df_raw[predclass]
        df_raw = df1
        df_raw[df_raw.columns.drop(predclass)]=np.log2(df_raw.drop(predclass,axis=1).astype('float')+1).values 
    df_raw.iloc[:,:-1]=pd.DataFrame(ro.r.normalizeBetweenArrays\
                            (df_raw.drop([predclass],axis=1).T)).T.values
    df_raw[predclass] = [i.replace(' ','_') for i in df_raw[predclass]]
    return df_raw

def gene_analysis(app,matrix_files, predclass, con, test,clinical,quantitative,split, testP=0.25,\
                    imputeNA=1, filt_outlier=1):
    director = os.path.join(app.config['GENE_OUTPUT_PATH'],matrix_files)
    ro.r('''library(limma)
            ''')
    
    ##加载并合并多项研究  
    pandas2ri.activate()
    matrix_files_ = os.listdir(f'{app.config["GENE_UPLOAD_PATH"]}/{matrix_files}')
    
    if clinical[0] != 'None':
        for i in clinical:
            print(clinical)
            matrix_files_.remove(i)
        for i in range(len(clinical)):
            clinical_file = f'{app.config["GENE_UPLOAD_PATH"]}/{matrix_files}/{clinical[i]}'
            if clinical_file.endswith('xlsx')|clinical_file.endswith('xls'):
                df_cli = pd.read_excel(clinical_file,index_col = 0)
            elif clinical_file.endswith('txt'):
                df_cli = pd.read_table(clinical_file,index_col = 0)
            else:
                df_cli = pd.read_csv(clinical_file,index_col = 0)
            if i==0:
                alldf=df_cli
            else:
                alldf=pd.concat([alldf,df_cli],join ='inner')
        df_cli = alldf
        
    # 划分训练集测试集
    if test != 'None':
        matrix_files_.remove(test)
        for i in range(len(matrix_files_)):
            matrix_file = f'{app.config["GENE_UPLOAD_PATH"]}/{matrix_files}/{matrix_files_[i]}'
            if matrix_file.endswith('xlsx')|matrix_file.endswith('xls'):
                df_raw = pd.read_excel(matrix_file,index_col = 0)
            elif matrix_file.endswith('txt'):
                df_raw = pd.read_table(matrix_file,index_col = 0)
            else:
                df_raw = pd.read_csv(matrix_file,index_col = 0)
            df_raw = gene_data_preprocess(df_raw, predclass, con, test, imputeNA=1, filt_outlier=1)
            df_raw['batchType'] = i+1
            ## 合并多项研究
            if i==0:
                alldf=df_raw
            else:
                alldf=pd.concat([alldf,df_raw],join ='outer')
        df_raw = alldf
        if len(matrix_files_)>1:
            ##过滤空值过多的gene，剩余以0填充
            if filt_outlier == 1:
                dfNa=pd.isnull(df_raw)
                cutoff = 0.1
                df_raw=df_raw.loc[:,dfNa.mean()<=cutoff]
            if imputeNA == 1:
                df_raw.fillna(0,inplace=True)
             ## R包sva中的combat函数矫正批次效应
            df_raw[df_raw.columns.drop([predclass,'batchType'])]=pd.DataFrame\
            (ro.r.removeBatchEffect(df_raw.drop([predclass,'batchType'],axis=1)\
                                    .T.values,df_raw['batchType'])).T.values

        df_raw = df_raw.drop(['batchType'],axis=1)
        
        df_train = df_raw
        test = f'{app.config["GENE_UPLOAD_PATH"]}/{matrix_files}/{test}'
        if test.endswith('xlsx')|test.endswith('xls'):
            df_test = pd.read_excel(test,index_col = 0)
        elif test.endswith('txt'):
            df_test = pd.read_table(test,index_col = 0)
        else:
            df_test = pd.read_csv(test,index_col = 0)
        df_test = gene_data_preprocess(df_test, predclass, con, test, imputeNA=1, filt_outlier=1)
    else:
        for i in range(len(matrix_files_)):
            matrix_file = f'{app.config["GENE_UPLOAD_PATH"]}/{matrix_files}/{matrix_files_[i]}'
            if matrix_file.endswith('xlsx')|matrix_file.endswith('xls'):
                df_raw = pd.read_excel(matrix_file,index_col = 0).T
            elif matrix_file.endswith('txt'):
                df_raw = pd.read_table(matrix_file,index_col = 0).T
            else:
                df_raw = pd.read_csv(matrix_file,index_col = 0).T
            df_train = pd.concat([df_raw[df_raw[predclass]==list(set(df_raw[predclass]))[0]].sample(frac=.75),\
                              df_raw[df_raw[predclass]==list(set(df_raw[predclass]))[1]].sample(frac=.75)])
            df_test = pd.concat([df_raw,df_train]).drop_duplicates(keep=False).T
            df_train = df_train.T
            df_train = gene_data_preprocess(df_train, predclass, con, test, imputeNA=1, filt_outlier=1)
            df_test = gene_data_preprocess(df_test, predclass, con, test, imputeNA=1, filt_outlier=1)
            df_train['batchType'] = i+1
            df_test['batchType'] = i+1
            ## 合并多项研究
            if i==0:
                alltrain=df_train
                alltest = df_test
            else:
                alltrain=pd.concat([alltrain,df_train],join ='outer')
                alltest = pd.concat([alltest,df_test],join ='outer')
        df_train = alltrain
        df_test = alltest
        if len(matrix_files_)>1:
            ##过滤空值过多的gene，剩余以0填充
            if filt_outlier == 1:
                dfNa=pd.isnull(df_train)
                cutoff = 0.1
                df_train=df_train.loc[:,dfNa.mean()<=cutoff]
                dfNa=pd.isnull(df_test)
                cutoff = 0.1
                df_test=df_test.loc[:,dfNa.mean()<=cutoff]
            if imputeNA == 1:
                df_train.fillna(0,inplace=True)
                df_test.fillna(0,inplace=True)
             ## R包sva中的combat函数矫正批次效应
            df_train[df_train.columns.drop([predclass,'batchType'])]=pd.DataFrame\
            (ro.r.removeBatchEffect(df_train.drop([predclass,'batchType'],axis=1).T.values,\
                                    df_train['batchType'])).T.values
            df_test[df_test.columns.drop([predclass,'batchType'])]=pd.DataFrame\
            (ro.r.removeBatchEffect(df_test.drop([predclass,'batchType'],axis=1).T.values,\
                                    df_test['batchType'])).T.values

        df_train = df_train.drop(['batchType'],axis=1)
        df_test = df_test.drop(['batchType'],axis=1)
    df_train.T.to_csv(path_or_buf=f'{director}/df_train.csv')
    df_test.T.to_csv(path_or_buf=f'{director}/df_test.csv')
    #卡方检验
    if clinical[0] != 'None':
        df_cli = df_cli.loc[list(df_train.index)+list(df_test.index),]
        df_cli['dataset'] = ['train']*df_train.shape[0]+['test']*df_test.shape[0]
        df_cli[predclass] = list(df_train[predclass])+list(df_test[predclass])
        for i,j in zip(quantitative,split):
            j = j.strip(' ').strip('-')
            j = j.split('-')
            j = [int(x.strip()) for x in j if x != '']
            j.sort()
            j = [0]+j+[10000]
            df_cli[i] = [str(x) for x in pd.cut(df_cli[i],j)]
        df_cli.dropna(inplace=True)
        cli_train = df_cli.loc[df_train.index,:]
        cli_test = df_cli.loc[df_test.index,:]
        df_cli.to_csv(path_or_buf=f'{director}/clinical.csv')
        chi2contingency = pd.DataFrame(index = ['cli_class','cli_dataset','train_class','test_class'],\
                                       columns=df_cli.columns[:-2])
        for i in df_cli.columns[:-2]:
            chi2data = pd.crosstab(df_cli[predclass],df_cli[i])
            chi2contingency.loc['cli_class',i] = stats.chi2_contingency(chi2data)[1]
            chi2data = pd.crosstab(df_cli['dataset'],df_cli[i])
            chi2contingency.loc['cli_dataset',i] = stats.chi2_contingency(chi2data)[1]
            chi2data = pd.crosstab(cli_train[predclass],cli_train[i])
            chi2contingency.loc['train_class',i] = stats.chi2_contingency(chi2data)[1]
            chi2data = pd.crosstab(cli_test[predclass],cli_test[i])
            chi2contingency.loc['test_class',i] = stats.chi2_contingency(chi2data)[1]
        chi2contingency.to_csv(path_or_buf=f'{director}/chi2_contingency.csv')

    disease = list(set(df_train[predclass]))
    disease.remove(con)
    disease = disease[0]
    ## limma差异分析
    ro.globalenv['rpred'] = predclass
    ro.globalenv['rcon'] = con.replace(' ', '_')
    ro.globalenv['rdisease'] = disease.replace(' ', '_')
    ro.globalenv['rdf_train'] = pandas2ri.py2rpy(df_train)
    ro.r('''
    design <- as.data.frame(model.matrix(~0+factor(rdf_train[,rpred])))
    rownames(design) <- rownames(rdf_train)
    colnames(design) <- unique(rdf_train[,rpred])
    data <- as.data.frame(t(rdf_train[,-match(rpred,colnames(rdf_train))]))
    colnames(data) <- rownames(rdf_train[,-match(rpred,colnames(rdf_train))])
    rownames(data) <- colnames(rdf_train[,-match(rpred,colnames(rdf_train))])
    fit <- lmFit(data,design)
    cont.matrix <- makeContrasts(paste0(rdisease,'-',rcon),levels=design)
    fit2 <- contrasts.fit(fit,cont.matrix)
    fit2 <- eBayes(fit2)
    trainstat <- topTable(fit2,adjust = 'fdr',number=200000)
    ''')
    trainstat = ro.r('trainstat')
    trainstat.to_csv(path_or_buf=f'{director}/trainstat.csv')
    
    ro.globalenv['rdf_test'] = pandas2ri.py2rpy(df_test)
    ro.r('''
    design <- as.data.frame(model.matrix(~0+factor(rdf_test[,rpred])))
    rownames(design) <- rownames(rdf_test)
    colnames(design) <- unique(rdf_test[,rpred])
    data <- as.data.frame(t(rdf_test[,-match(rpred,colnames(rdf_test))]))
    colnames(data) <- rownames(rdf_test[,-match(rpred,colnames(rdf_test))])
    rownames(data) <- colnames(rdf_test[,-match(rpred,colnames(rdf_test))])
    fit <- lmFit(data,design)
    cont.matrix <- makeContrasts(paste0(colnames(design)[1],'-',colnames(design)[2]),levels=design)
    fit2 <- contrasts.fit(fit,cont.matrix)
    fit2 <- eBayes(fit2)
    teststat <- topTable(fit2,adjust = 'fdr',number=200000)
    ''')
    teststat = ro.r('teststat')
    teststat.to_csv(path_or_buf=f'{director}/teststat.csv')
    
    diffgene = list(trainstat.index[(trainstat['adj.P.Val']<0.05) & (abs(trainstat['logFC'])>math.log2(1.2))]\
                    .intersection(teststat.index[(teststat['adj.P.Val']<1) &(abs(teststat['logFC'])>=math.log2(1))]))
    if len(diffgene)<10:
        diffgene = list(trainstat.index[(trainstat['P.Value']<0.05) & (abs(trainstat['logFC'])>math.log2(1.2))]\
                    .intersection(teststat.index[(teststat['P.Value']<1) &(abs(teststat['logFC'])>=math.log2(1))]))
    # diffgene = list(trainstat.index[(trainstat['adj.P.Val']<0.001) & (abs(trainstat['logFC'])>math.log2(2))])
    df_train[diffgene+[predclass]].T.to_csv(path_or_buf=f'{director}/df_diff.csv')
    
    ###特征选择
    
    df_train[predclass] = [0 if x == con else 1 for x in df_train[predclass]]
    df_test[predclass] = [0 if x == con else 1 for x in df_test[predclass]]
    
    x_train,x_test = df_train[diffgene],df_test[diffgene]
    y_train,y_test = df_train[predclass],df_test[predclass]
    # x_train[:] = preprocess(x_train.values,scale_by_minmax=True)
    # x_test[:] = preprocess(x_test.values,scale_by_minmax=True)
    # x_train = var_filter(0.1,x_train)
    # x_test = x_test[x_train.columns]
    # x_train[:] = preprocess(x_train.values,scale_by_minmax=False)
    # x_test[:] = preprocess(x_test.values,scale_by_minmax=False)
    ##wrapper
    # df_dat = Recursive_Feature_Elimination(df_train[diffgene], predclass)
    # features_slct = list(df_dat.columns.drop(predclass))
    # df_dat[predclass] = [0 if x == con else 1 for x in df_dat[predclass]]
    ##SelectFromModel
    features_slct,model,result,modelname,best_params,model_probs = sel_f_model(x_train,y_train,x_test,y_test)
    result.to_csv(path_or_buf=f'{director}/feature_slct.csv')
    best_params.to_csv(path_or_buf=f'{director}/best_params.csv')
    df_dat = df_train[features_slct+[predclass]]

    x_train, x_test = df_train[features_slct],df_test[features_slct]
    y_train,y_test = df_train[predclass],df_test[predclass]
    
    trainauc = {}
    testauc = {}
    for gene in diffgene:
        train_evaluate = Biomarker_Evaluate(df_train[predclass],df_train[gene].astype('float'))
        test_evaluate = Biomarker_Evaluate(df_test[predclass],df_test[gene].astype('float'))
        
        trainpred_roc,trainpred_pr,trainauc[gene] = train_evaluate.evaluate()
        testpred_roc,testpred_pr,testauc[gene] = test_evaluate.evaluate()
    
    trainroc = pd.DataFrame(trainauc).T
    testroc = pd.DataFrame(testauc).T
    df_roc = pd.concat([trainroc,testroc],axis =1,keys=['train','test'])
    df_features = pd.concat([trainstat.loc[features_slct,['logFC','adj.P.Val']],\
                             teststat.loc[features_slct,['logFC','adj.P.Val']]],axis=1,keys = ['train','test'])
    df_features = pd.concat([df_features,df_roc.loc[features_slct,:]],axis =1)
    df_roc.to_csv(path_or_buf=f'{director}/df_roc.csv')
    df_features.to_csv(path_or_buf=f'{director}/df_features.csv')
    
    return df_dat, x_train, x_test, y_train, y_test,model,modelname,disease,model_probs

def data_preprocess(df_raw, predclass, con, imputeNA=1, filt_outlier=1):
      
    ## 按照每个观察值缺失或为0的特征数，按百分位法去离群值
    if filt_outlier == 1:
        # within 2 standerd deviation
        # df_dat[np.abs(df_dat.A - df_dat.A.mean()) <= (2 * df_dat.A.std())]
        countsNa=~((df_raw.drop([predclass],axis=1)==0) ==\
               (pd.isnull(df_raw.drop([predclass],axis=1)))).sum(axis=1)
        iqr=np.percentile(countsNa,75,axis=0)-np.percentile(countsNa,25,axis=0)
        lowerlimt,upperlimt = np.percentile(countsNa,25,axis=0)-1.5*iqr,np.percentile(countsNa,75,axis=0)+1.5*iqr
        df_raw=df_raw[[x<=upperlimt for x in countsNa]]
    ## KNN法补值
    if imputeNA == 1:
        imp_knn = KNNImputer(n_neighbors=5)
        df_raw = pd.concat([pd.DataFrame(imp_knn.fit_transform(df_raw.drop([predclass],axis=1)),\
                  columns = df_raw.columns.drop(predclass),index=df_raw.index),df_raw[predclass]],axis=1)
    ## 标准化
    # df_raw[:] = StandardScaler().fit_transform(df_raw)
    # df_raw[:] = RobustScaler().fit_transform(df_raw)
    # df_raw[:] = MaxAbsScaler().fit_transform(df_raw)
    # df_raw[:] = MinMaxScaler().fit_transform(df_raw)

    return df_raw

def matrix_analysis(app,matrix_files, predclass, con, testP=0.25, imputeNA=1, filt_outlier=1):
    director = os.path.join(app.config['MATRIX_OUTPUT_PATH'],matrix_files)
    
    ##加载并合并多项研究  
    matrix_files_ = os.listdir(f'{app.config["MATRIX_UPLOAD_PATH"]}/{matrix_files}')
    
    # 划分训练集测试集
    for i in range(len(matrix_files_)):
        matrix_file = f'{app.config["MATRIX_UPLOAD_PATH"]}/{matrix_files}/{matrix_files_[i]}'
        if matrix_file.endswith('xlsx')|matrix_file.endswith('xls'):
            df_raw = pd.read_excel(matrix_file)
        elif matrix_file.endswith('txt'):
            df_raw = pd.read_table(matrix_file)
        else:
            df_raw = pd.read_csv(matrix_file)
        if isinstance(df_raw.iloc[0,0],str):
            df_raw.index = df_raw.iloc[:,0]
            df_raw = df_raw.iloc[:,1:]
        df_train = pd.concat([df_raw[df_raw[predclass]==list(set(df_raw[predclass]))[0]].sample(frac=.75),\
                          df_raw[df_raw[predclass]==list(set(df_raw[predclass]))[1]].sample(frac=.75)])
        df_test = pd.concat([df_raw,df_train]).drop_duplicates(keep=False)
    df_train.to_csv(path_or_buf=f'{director}/df_train.csv')
    df_test.to_csv(path_or_buf=f'{director}/df_test.csv')
    df_train = data_preprocess(df_train, predclass, con, imputeNA=1, filt_outlier=1)
    df_test = data_preprocess(df_test, predclass, con, imputeNA=1, filt_outlier=1)
    ###特征选择
    while True:
        try:
            con = int(con)
            disease = list(set(df_train[predclass]))
            disease.remove(con)
            disease = disease[0]
            break
        except ValueError:
            disease = list(set(df_train[predclass]))
            disease.remove(con)
            disease = disease[0]
            df_train[predclass] = [0 if x == con else 1 for x in df_train[predclass]]
            df_test[predclass] = [0 if x == con else 1 for x in df_test[predclass]]
            break

    x_train,x_test = df_train.drop(predclass,axis=1),df_test.drop(predclass,axis=1)
    y_train,y_test = df_train[predclass],df_test[predclass]
    
    # x_train[:] = preprocess(x_train.values,scale_by_minmax=True)
    # x_test[:] = preprocess(x_test.values,scale_by_minmax=True)
    # x_train = var_filter(0.1,x_train)
    # x_test = x_test[x_train.columns]
    # x_train[:] = preprocess(x_train.values,scale_by_minmax=False)
    # x_test[:] = preprocess(x_test.values,scale_by_minmax=False)
    
    ##wrapper
    # df_dat = Recursive_Feature_Elimination(df_train, predclass)
    # features_slct = list(df_dat.columns.drop(predclass))
    # df_dat[predclass] = [0 if x == con else 1 for x in df_dat[predclass]]
    ##SelectFromModel
    features_slct,model,result,modelname,best_params,model_probs = sel_f_model(x_train,y_train,x_test,y_test)
    result.to_csv(path_or_buf=f'{director}/feature_slct.csv')
    best_params.to_csv(path_or_buf=f'{director}/best_params.csv')
    df_dat = df_train[features_slct+[predclass]]

    x_train, x_test = df_train[features_slct],df_test[features_slct]
    y_train,y_test = df_train[predclass],df_test[predclass]
    
    trainauc = {}
    testauc = {}
    for x in df_train.columns.drop(predclass):
        train_evaluate = Biomarker_Evaluate(df_train[predclass],df_train[x].astype('float'))
        test_evaluate = Biomarker_Evaluate(df_test[predclass],df_test[x].astype('float'))
        
        trainpred_roc,trainpred_pr,trainauc[x] = train_evaluate.evaluate()
        testpred_roc,testpred_pr,testauc[x] = test_evaluate.evaluate()
    
    trainroc = pd.DataFrame(trainauc).T
    testroc = pd.DataFrame(testauc).T
    df_roc = pd.concat([trainroc,testroc],axis =1,keys=['train','test'])
    df_roc.to_csv(path_or_buf=f'{director}/df_roc.csv')
    
    return df_dat, x_train, x_test, y_train, y_test,model,modelname,disease,model_probs