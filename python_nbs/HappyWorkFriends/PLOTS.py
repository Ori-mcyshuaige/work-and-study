from .InviteFriends import *
import sklearn.metrics as metrics

def Plot_feature_importance(dataset_con_enc, predclass, top=50, ax=None):
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier()
    clf.fit(dataset_con_enc.drop(predclass, axis=1), dataset_con_enc[predclass])

    plt.style.use('seaborn-whitegrid')
    importance = clf.feature_importances_
    importance = pd.DataFrame(importance, index=dataset_con_enc.drop(predclass, axis=1).columns, columns=["Importance"])
    importance.sort_values(by='Importance', ascending=True).head(top).plot(kind='barh', figsize=(20,len(importance)/2), ax=ax)
    
    return importance.sort_values(by='Importance', ascending=False)

def Heatmap_Correlation(dataset_con_enc):
    plt.style.use('seaborn-whitegrid')
    fig = plt.figure(figsize=(12,10))
    mask = np.zeros_like(dataset_con_enc.corr(), dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(dataset_con_enc.corr(), 
                vmin=-1, vmax=1, 
                square=True, 
                cmap=sns.color_palette("RdBu_r", 100), 
                mask=mask, 
                linewidths=.5)

def Plot_Linear_Relationship(df, x, y, xlim=None, ylim=None, title='', do_R2=True, **kwargs):
    with sns.plotting_context('paper', font_scale=1.8):
        g = sns.jointplot(x=x, y=y, data=df, kind='reg', xlim=xlim, ylim=ylim, stat_func=scs.pearsonr, annot_kws=dict(frameon=False), **kwargs)
        r2 = r2_score(df[y], df[x])
        #r2 = scs.pearsonr(df[x], df[y])[0]**2
        if title: g.ax_marg_x.set_title(title)
        if do_R2: g.ax_marg_x.set_title(title + ' $R^2$: %.3f' %(r2))
        plt.show()
        
def plot_roc_curve(y_test, preds):
    fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
    roc_auc = metrics.auc(fpr, tpr)
    plt.title('Receiver Operating Characteristic')
    plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.show()
        
def Plot_ROC_Curves(models, probs, y_test, ax=None):
    plt.style.use('seaborn-whitegrid')
    fig = plt.figure(figsize=(8,8))
    colors = [
    'blue',
    'green',
    'red',
    'cyan',
    'magenta',
    'yellow',
    ]

    plt.title('Receiver Operating Characteristic')
    plt.plot([0, 1], [0, 1], color='gray', linestyle='dashed')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')

    def plot_roc_curves(y_test, prob, model, ax=None):
        fpr, tpr, threshold = metrics.roc_curve(y_test, prob)
        roc_auc = metrics.auc(fpr, tpr)
        plt.plot(fpr, tpr, 'b', label = model + ' AUC = %0.2f' % roc_auc, color=colors[i], ax=ax)
        plt.legend(loc = 'lower right')

    for i, model in list(enumerate(models)):
        plot_roc_curves(y_test, probs[i], model, ax)

    plt.show()
    
def Plot_ROC_Curves_Web(models, probs, y_test, ax=None):
    plt.style.use('seaborn-whitegrid')
    fig = plt.figure(figsize=(10,10))
    colors = [
    'blue',
    'green',
    'red',
    'cyan',
    'magenta',
    'yellow',
    ]

    plt.title('Receiver Operating Characteristic')
    plt.plot([0, 1], [0, 1], color='black', linestyle='dashed')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')

    def plot_roc_curves(y_test, prob, model, ax=None):
        fpr, tpr, threshold = metrics.roc_curve(y_test, prob)
        roc_auc = metrics.auc(fpr, tpr)
        ax.plot(fpr, tpr, 'b', label = model + ' AUC = %0.2f' % roc_auc, color=colors[i])
        ax.legend(loc = 'lower right')

    for i, model in list(enumerate(models)):
        plot_roc_curves(y_test, probs[i], model, ax)