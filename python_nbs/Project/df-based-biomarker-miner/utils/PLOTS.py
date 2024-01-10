import matplotlib as mpl
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import seaborn as sns
import sklearn.metrics as metrics
from statannotations.Annotator import Annotator

def Plot_feature_importance(dataset_con_enc, predclass, top=50, ax=None):
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier()
    clf.fit(dataset_con_enc.drop(predclass, axis=1), dataset_con_enc[predclass])

    plt.style.use('seaborn-whitegrid')
    importance = clf.feature_importances_
    importance = pd.DataFrame(importance, index=dataset_con_enc.drop(predclass, axis=1).columns, columns=["Importance"])
    importance.sort_values(by='Importance', ascending=True).head(top).plot(kind='barh', figsize=(12,12), ax=ax)
    
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
        
def Plot_ROC_Curves(models, probs, y_test):
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

    def plot_roc_curves(y_test, prob, model):
        fpr, tpr, threshold = metrics.roc_curve(y_test, prob)
        roc_auc = metrics.auc(fpr, tpr)
        plt.plot(fpr, tpr, 'b', label = model + ' AUC = %0.3f' % roc_auc, color=colors[i])
        plt.legend(loc = 'lower right')

    for i, model in list(enumerate(models)):
        plot_roc_curves(y_test, probs[i], model)

    plt.show()
    
def Plot_ROC_Curves_Web(models, probs, y_test, ax=None):
    plt.style.use('seaborn-whitegrid')
    fig = plt.figure(figsize=(12,12))
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
        ax.plot(fpr, tpr, 'b', label = model + ' AUC = %.5f' % roc_auc, color=colors[i])
        ax.legend(loc = 'lower right')

    for i, model in list(enumerate(models)):
        plot_roc_curves(y_test, probs[i], model, ax)

def plot_box(df,x,y,group,path):
    fig,ax = plt.subplots(figsize=(5,4),dpi=600,facecolor="w")
    ax = sns.boxplot(data=df, x=x, y=y,ax=ax,hue=group,fliersize=0,saturation=0.75,whis=1.5,width=0.8)
    ax = sns.stripplot(x=x, y=y, hue=group, data=df, dodge=True, alpha=0.6, ax=ax)#palette=sns.color_palette("plasma_r",n_colors=4)

    pairs=[((i,list(set(df[group]))[0]),(i,list(set(df[group]))[1])) for i in set(df[x])]


    annotator = Annotator(ax, pairs, data=df, x=x, y=y, hue=group)
    annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=1,loc='outside',verbose=2)
    annotator.apply_and_annotate()

    ax.tick_params(which='major',direction='in',length=3,width=1.,labelsize=15,bottom=False)
    for spine in ["top","left","right"]:
        ax.spines[spine].set_visible(True)
    ax.spines[['bottom',"top","left","right"]].set_linewidth(1.5)
    # ax.grid(axis='y',ls='--',c='gray')
    ax.legend(loc = [1,.7])
    ax.grid(False)
    ax.set_axisbelow(False)
    plt.ylabel(None)
    plt.xlabel(None)
    
    plt.savefig(path,bbox_inches='tight',dpi=600)


def Plot_tSNE(df_plot, prefix, opath='/data/users/danxu/Data/imgs', class_label='pred'):
    from sklearn.manifold import TSNE
    df_plot[class_label] = df_plot[class_label].astype('category')
    X_std = df_plot.drop(class_label, axis=1).values
    Target = df_plot[class_label]
    tsne = TSNE(n_components=2)
    tsne_results = tsne.fit_transform(X_std)
    plt.figure()
    sns.scatterplot(x=tsne_results[:,0], y=tsne_results[:,1],
                    palette=sns.color_palette("hls", len(Target.unique())),
                    hue=Target, alpha=0.8)
    fig_file = f'{opath}/{prefix}.t-SNE.png'
    plt.savefig(fig_file, dpi=500 ,bbox_inches='tight')
    return(fig_file)

