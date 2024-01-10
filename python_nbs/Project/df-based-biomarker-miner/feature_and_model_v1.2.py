#!/data/users/machiyu/anaconda3/envs/Flask/bin/python

from support_f_and_m import *

app = Flask(__name__)
app.secret_key = os.getenv('SECRET_KEY', 'secret string')
app.jinja_env.trim_blocks =True
app.jinja_env.lstrip_blocks = True

app.config['MATRIX_UPLOAD_PATH'] = os.path.join(app.root_path, 'upload_file/matrix')
app.config['GENE_UPLOAD_PATH'] = os.path.join(app.root_path, 'upload_file/gene')
app.config['MATRIX_OUTPUT_PATH'] = os.path.join(app.root_path, 'output_result/matrix')
app.config['GENE_OUTPUT_PATH'] = os.path.join(app.root_path, 'output_result/gene')
app.config['MATRIX_ALLOWED_EXTENSIONS'] = ['xls','xlsx','csv','tsv']

os.makedirs(app.config['MATRIX_UPLOAD_PATH'],exist_ok = True)
os.makedirs(app.config['GENE_UPLOAD_PATH'],exist_ok = True)
os.makedirs(app.config['MATRIX_OUTPUT_PATH'],exist_ok = True)
os.makedirs(app.config['GENE_OUTPUT_PATH'],exist_ok = True)

@app.route('/', methods = ['GET', 'POST'])
def index():
    return render_template('index.html')

@app.route('/gene', methods = ['GET', 'POST'])
def gene():
    form = geneUpload()
    if form.validate_on_submit():
        f = form.file.data
        filename = secure_filename(f.filename)
        f.save(os.path.join(app.config['GENE_UPLOAD_PATH'],filename))
        zip_file = zipfile.ZipFile(f)
        for names in zip_file.namelist():  # 解压
            zip_file.extract(names, path=os.path.join(app.config['GENE_UPLOAD_PATH'],\
                                                      '.'.join(filename.split('.')[:-1]))) # 解压后名称
        zip_file.close()
        if os.path.exists(os.path.join(app.config['GENE_UPLOAD_PATH'],filename)):
            os.remove(os.path.join(app.config['GENE_UPLOAD_PATH'],filename))
        if not os.path.exists(os.path.join(app.config['GENE_OUTPUT_PATH'],'.'.join(filename.split('.')[:-1]))):
            os.makedirs(os.path.join(app.config['GENE_OUTPUT_PATH'],\
                                     '.'.join(filename.split('.')[:-1])),exist_ok = True)
        flash('Upload success')
        session['filenames'] = [filename]
        predclass = form.predclass.data
        test = form.test.data
        con = form.con.data
        clinical = form.clinical.data.strip(' ').strip(';')
        quantitative = form.quantitative.data.strip(' ').strip(';')
        split = form.split.data.strip(' ').strip(';')
        return redirect(url_for('gene_feature_select',dat_files=filename,predclass=predclass,con=con,test = test,\
                               clinical =clinical,quantitative=quantitative,split=split))
    return render_template('gene_upload.html', form=form)

@app.route('/gene_feature_select/<dat_files>,<predclass>,<con>,<test>,<clinical>,<quantitative>,<split>',methods = ['GET','POST'])
# @app.route('/matrix_feature_select',defaults={'dat_files': None,'predclass':'pred','con':'control','test':None})
def gene_feature_select(dat_files, predclass, con, test,clinical,quantitative,split, testP=0.25, imputeNA=1, filt_outlier=1):
    form = download()
    if form.validate_on_submit():
        return redirect(url_for('download_file',filename = dat_files))
    dat_file='.'.join(dat_files.split('.')[:-1])
    director = os.path.join(app.config['GENE_OUTPUT_PATH'],f'{dat_file}/result')
    if not os.path.exists(director):
        os.makedirs(director,exist_ok = True)

    clinical = clinical.split(';')
    clinical = [i.strip() for i in clinical if i != '']
    quantitative = quantitative.split(';')
    quantitative = [i.strip() for i in quantitative if i != '']
    split = split.split(';')
    split = [i.strip() for i in split if i != '']

    df_dat, X_train, X_test, y_train, y_test,model,modelname,disease,model_probs = gene_analysis\
    (app,dat_file, predclass, con, test,clinical,quantitative,split, testP=0.25, imputeNA=1, filt_outlier=1)

    fig = Figure(figsize=(16,32),dpi=600)
    axes = fig.subplots(nrows=1, ncols=2)
    #missingno.matrix(df_dat, figsize = (30,5), ax=axes[0])

    top50_features = Plot_feature_importance(df_dat, predclass, ax=axes[0])

    # selection model
    model_evaluate = {}

    md_test = Biomarker_Evaluate(y_test,model_probs[modelname])
    pred_roc,pred_pr,model_evaluate[modelname.value] = md_test.evaluate()
    roc_confusion_matrix = md_test.confusion_matrix_(pred_roc)
    pr_confusion_matrix = md_test.confusion_matrix_(pred_pr)
    roc_confusion_matrix.to_csv(path_or_buf = f'{director}/{modelname.value}_roc_confusion_matrix.csv')
    pr_confusion_matrix.to_csv(path_or_buf = f'{director}/{modelname.value}_pr_confusion_matrix.csv')
    model_evaluate = pd.DataFrame(model_evaluate).T
    model_evaluate.to_csv(path_or_buf = f'{director}/model_evaluate.csv')
    md_test.plot_roc(f'{director}/{modelname.value}_roc_plot.png')
    md_test.plot_prroc(f'{director}/{modelname.value}_pr_plot.png')
    md_test.plot_confusion_matrix(pred_roc,[con,disease],f'{director}/{modelname.value}_roc_confusion_matrix.png')
    md_test.plot_confusion_matrix(pred_pr,[con,disease],f'{director}/{modelname.value}_pr_confusion_matrix.png')

    Plot_ROC_Curves_Web([i.value for i in list(model_probs.keys())], list(model_probs.values()), y_test, ax=axes[1])

    directory = f"{app.config['GENE_OUTPUT_PATH']}/{dat_file}"
    zip_file = zipfile.ZipFile(f"{directory}.zip",'w')
    for path, dirnames, filenames in os.walk(directory):
        zippath = path.replace(directory,'')
        if filenames:
            for filename in filenames:
                zip_file.write(os.path.join(path, filename), os.path.join(zippath, filename),\
                               compress_type=zipfile.ZIP_DEFLATED)
    zip_file.close()

    ##返回一个demo
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # buf.seek(0)
    # Embed the result in the html output.
    img = base64.b64encode(buf.getbuffer()).decode("ascii")
    # data = base64.b64encode(buf.getvalue()).decode()
    return render_template('file_select.html',importance_data = top50_features.T.to_html(),img = img,form = form)

@app.route('/patent_assistant',methods = ['GET','POST'])
def patent_assistant():
    form = patentAssistant()
    if form.validate_on_submit():
        pname = form.project.data
        gene = form.geneID.data.strip(' ').strip(';')
        analysis = form.analysis.data
        begin = int(form.combinationbegin.data)
        end = int(form.combinationend.data)
        predclass = form.predclass.data
        con = form.con.data
        return redirect(url_for('plot_combination_roc_box',pname=pname,analysis=analysis,\
                                gene=gene,predclass=predclass,con=con,begin=begin,end=end))
    return render_template('patent_assistant.html', form=form)

@app.route('/plot_combination_roc_box/<pname>,<analysis>,<gene>,<predclass>,<con>,<begin>,<end>',methods = ['GET','POST'])
def plot_combination_roc_box(pname,analysis,gene,predclass,con,begin,end):
    path = f"{app.config['GENE_OUTPUT_PATH']}/{pname}"
    path = path.replace('gene',analysis)
    savepath = os.path.join(path,gene)
    savepath = savepath.replace(';','&')
    if not os.path.exists(savepath):
        os.makedirs(savepath,exist_ok = True)
    gene = gene.split(';')
    gene = [i.strip() for i in gene if i != '']
    df_train = pd.read_csv(f'{path}/df_train.csv',index_col=0).T
    df_test = pd.read_csv(f'{path}/df_test.csv',index_col=0).T
    disease = list(set(df_train[predclass]))
    disease.remove(con)
    disease = disease[0]
    y_train,y_test = df_train[predclass],df_test[predclass]
    df_train,df_test = df_train[gene+[predclass]],df_test[gene+[predclass]]
    for i,j in zip(df_train.columns.drop(predclass),df_test.columns.drop(predclass)):
        databox = df_train[[i,predclass]]
        databox['genename'] = i
        databox.columns=['expression',predclass,'genename']
        databox['expression']=databox['expression'].astype('float64')
        databox[['genename',predclass]]=databox[['genename',predclass]].astype('category')
        fig = plot_box(databox,'genename','expression',predclass,f'{savepath}/{i}_train.png')
        databox = df_train[[j,predclass]]
        databox['genename'] = j
        databox.columns=['expression',predclass,'genename']
        databox['expression']=databox['expression'].astype('float64')
        databox[['genename',predclass]]=databox[['genename',predclass]].astype('category')
        plot_box(databox,'genename','expression',predclass,f'{savepath}/{j}_test.png')
    y_train = [0 if x == con else 1 for x in y_train]
    y_test = [0 if x == con else 1 for x in y_test]
    model_evaluate = {}
    for i in range(int(begin),int(end)+1):
        for features in it.combinations(gene,i):
            features = list(features)
            fname = ";".join(features)
            x_train,x_test = df_train[features],df_test[features]
            if len(features)>1:
                models = {"Logistic Regression": LogisticRegression()}
                model_scores, model_probs, best_params = fit_and_score(x_train, x_test, y_train, y_test, models)
                biomarker_evaluate = Biomarker_Evaluate(y_test,model_probs["Logistic Regression"])
                pred_roc,pred_pr,model_evaluate[fname] = biomarker_evaluate.evaluate()
                biomarker_evaluate.plot_roc(f'{savepath}/{fname}_roc_plot.png')
                biomarker_evaluate.plot_prroc(f'{savepath}/{fname}_pr_plot.png')
                biomarker_evaluate.plot_confusion_matrix(pred_roc,[con,disease],\
                                                     f'{savepath}/{fname}_roc_confusion_matrix.png')
                biomarker_evaluate.plot_confusion_matrix(pred_pr,[con,disease],\
                                                         f'{savepath}/{fname}_pr_confusion_matrix.png')
            else:
                biomarker_evaluate = Biomarker_Evaluate(y_test,x_test.iloc[:,0])
                pred_roc,pred_pr,model_evaluate[fname] = biomarker_evaluate.evaluate()
                biomarker_evaluate.plot_roc(f'{savepath}/{fname}_roc_plot.png')
                biomarker_evaluate.plot_prroc(f'{savepath}/{fname}_pr_plot.png')
                biomarker_evaluate.plot_confusion_matrix(pred_roc,[con,disease],\
                                                         f'{savepath}/{fname}_roc_confusion_matrix.png')
                biomarker_evaluate.plot_confusion_matrix(pred_pr,[con,disease],\
                                                         f'{savepath}/{fname}_pr_confusion_matrix.png')
    model_evaluate = pd.DataFrame(model_evaluate).T
    model_evaluate.to_csv(path_or_buf = f'{savepath}/all_result.csv')
    directory = savepath
    zip_file = zipfile.ZipFile(f"{directory}.zip".replace('/','@').replace('@','/',directory.count('/')-1),'w')
    for path, dirnames, filenames in os.walk(directory):
        zippath = path.replace(directory,'')
        if filenames:
            for filename in filenames:
                zip_file.write(os.path.join(path, filename), os.path.join(zippath, filename),\
                               compress_type=zipfile.ZIP_DEFLATED)
    zip_file.close()
    filename = f"{directory}.zip".replace('/','@').replace('@','/',directory.count('/')-1).replace(f"{app.config['GENE_OUTPUT_PATH']}/",'')
    return redirect(url_for('download_file',filename = filename))

@app.route('/matrix', methods = ['GET', 'POST'])
def matrix():
    form = matrixUpload()
    if form.validate_on_submit():
        f = form.file.data
        filename = secure_filename(f.filename)
        if not os.path.exists(os.path.join(app.config['MATRIX_OUTPUT_PATH'],'.'.join(filename.split('.')[:-1]))):
            os.makedirs(os.path.join(app.config['MATRIX_OUTPUT_PATH'],\
                                     '.'.join(filename.split('.')[:-1])),exist_ok = True)
        if not os.path.exists(os.path.join(app.config['MATRIX_UPLOAD_PATH'],'.'.join(filename.split('.')[:-1]))):
            os.makedirs(os.path.join(app.config['MATRIX_UPLOAD_PATH'],\
                                     '.'.join(filename.split('.')[:-1])),exist_ok = True)
        f.save(f"{app.config['MATRIX_UPLOAD_PATH']}/{'.'.join(filename.split('.')[:-1])}/{filename}")
        flash('Upload success')
        session['filenames'] = [filename]
        predclass = form.predclass.data
        con = form.con.data
        return redirect(url_for('matrix_feature_select',dat_files=filename,predclass=predclass,con=con))
    return render_template('matrix_upload.html', form=form)


# def matrix():
#     form = matrixUpload()
#     if form.validate_on_submit():
#         studyid=form.studyid.data
#         studyid=studyid.strip(';').strip()
#         study_result=os.path.join(app.config['MATRIX_OUTPUT_PATH'],studyid)
#         if os.path.exists(study_result):
#             return render_template('matrix_upload.html', form=form)
#         else:
#             os.makedirs(study_result,exist_ok = True)
#             for geoid in studyid.split(';'):
#                 geoid=geoid.strip()
#                 if not os.path.exists(os.path.join(app.config['MATRIX_UPLOAD_PATH'],geoid)):
#                     os.makedirs(os.path.join(app.config['MATRIX_UPLOAD_PATH'],geoid))
#                     while True:
#                         try:
#                             gse=GEOparse.get_GEO(geo=geoid, destdir=os.path.join(app.config['MATRIX_UPLOAD_PATH'],geoid),
#                                                    silent=False,include_data = True,annotate_gpl=True)
#                             break
#                         except OSError:
#                             continue

#     return render_template('matrix_upload.html', form=form)


@app.route('/matrix_feature_select/<dat_files>,<predclass>,<con>',methods = ['GET','POST'])
def matrix_feature_select(dat_files, predclass, con, testP=0.25, imputeNA=1, filt_outlier=1):
    form = download()
    if form.validate_on_submit():
        return redirect(url_for('download_file',filename = dat_files))
    dat_file='.'.join(dat_files.split('.')[:-1])
    director = os.path.join(app.config['MATRIX_OUTPUT_PATH'],f'{dat_file}/result')
    if not os.path.exists(director):
        os.makedirs(director,exist_ok = True)

    df_dat, X_train, X_test, y_train, y_test,model,modelname,disease,model_probs = matrix_analysis\
    (app,dat_file, predclass, con, testP=0.25,
     imputeNA=1, filt_outlier=1)

    fig = Figure(figsize=(16,32),dpi=600)
    axes = fig.subplots(nrows=1, ncols=2)
    #missingno.matrix(df_dat, figsize = (30,5), ax=axes[0])

    features_importance = Plot_feature_importance(df_dat, predclass, ax=axes[0])
    features_importance.to_csv(path_or_buf=f'{director}/features_importance.csv')

    # selection model
    model_evaluate = {}

    md_test = Biomarker_Evaluate(y_test,model_probs[modelname])
    pred_roc,pred_pr,model_evaluate[modelname.value] = md_test.evaluate()
    roc_confusion_matrix = md_test.confusion_matrix_(pred_roc)
    pr_confusion_matrix = md_test.confusion_matrix_(pred_pr)
    roc_confusion_matrix.to_csv(path_or_buf = f'{director}/{modelname.value}_roc_confusion_matrix.csv')
    pr_confusion_matrix.to_csv(path_or_buf = f'{director}/{modelname.value}_pr_confusion_matrix.csv')
    model_evaluate = pd.DataFrame(model_evaluate).T
    model_evaluate.to_csv(path_or_buf = f'{director}/model_evaluate.csv')
    md_test.plot_roc(f'{director}/{modelname.value}_roc_plot.png')
    md_test.plot_prroc(f'{director}/{modelname.value}_pr_plot.png')
    md_test.plot_confusion_matrix(pred_roc,[con,disease],f'{director}/{modelname.value}_roc_confusion_matrix.png')
    md_test.plot_confusion_matrix(pred_pr,[con,disease],f'{director}/{modelname.value}_pr_confusion_matrix.png')

    Plot_ROC_Curves_Web([i.value for i in list(model_probs.keys())], list(model_probs.values()), y_test, ax=axes[1])

    directory = f"{app.config['MATRIX_OUTPUT_PATH']}/{dat_file}"
    zip_file = zipfile.ZipFile(f"{directory}.zip",'w')
    for path, dirnames, filenames in os.walk(directory):
        zippath = path.replace(directory,'')
        if filenames:
            for filename in filenames:
                zip_file.write(os.path.join(path, filename), os.path.join(zippath, filename),\
                               compress_type=zipfile.ZIP_DEFLATED)
    zip_file.close()

    ##返回一个demo
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # buf.seek(0)
    # Embed the result in the html output.
    img = base64.b64encode(buf.getbuffer()).decode("ascii")
    # data = base64.b64encode(buf.getvalue()).decode()
    return render_template('file_select.html',importance_data = features_importance.T.to_html(),img = img,form = form)

@app.route('/deving', methods = ['GET', 'POST'])
def deving():
    return '<h1>我什么都没有做哈哈哈哈哈哈哈哈哈哈哈哈哈<h1>'

@app.route('/download_file/<filename>',methods = ['GET','POST'])
def download_file(filename):
    directory = app.config['GENE_OUTPUT_PATH']
    response = make_response(send_from_directory(directory, filename, as_attachment=True))
    response.headers["Content-Disposition"] = f"attachment; filename={filename.encode().decode('latin-1')}"
    return response

if __name__ == '__main__':
    app.run('10.10.10.199', 5001)