#!/usr/bin/env python
# coding: utf-8

# In[1]:


# get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

import sys

sys.path.append('code')

from sklearn import metrics
from sklearn import model_selection 
from sklearn import utils
from sklearn.metrics import precision_recall_fscore_support

from sklearn.dummy import DummyClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV

from imblearn import over_sampling

import warnings
warnings.filterwarnings("ignore") 

sns.set_context('paper')
sns.set_palette('colorblind')
sns.set_style('whitegrid')


# # Read in Formatted Data

# In[2]:


SR_encoder = {'S':0, 'I': 1, 'R': 1, 'RA': 1}

ecoli_ast = pd.read_csv('/Users/karatsang/Desktop/ML_analysis/ecoli/ast.tsv', sep='\t').set_index('Sample').replace(SR_encoder)
ecoli_rgi = pd.read_pickle('/Users/karatsang/Desktop/ML_analysis/ecoli/rgi_encoded.pkl')
pseudo_ast = pd.read_csv('/Users/karatsang/Desktop/ML_analysis/pseudomonas/ast.tsv', sep='\t').set_index('Sample').replace(SR_encoder).drop('ertapenem', axis=1)
pseudo_rgi = pd.read_pickle('/Users/karatsang/Desktop/ML_analysis/pseudomonas/rgi_encoded.pkl')


# # Choose Model Type using E. coli

# In[3]:


def fit_interpretable_clf_per_abx(rgi, ast, classifiers, params):
    """
    Fit intepretable models to each and get best overall performance
    """
    classifier_perf = {}

    for clf_name, clf in classifiers.items():

        for abx in ast.columns:
            # check balance if less than 10% are minority class then just train all
            props = ast[abx].value_counts()
            if props.min() / props.max() < 0.1:
                x_train_res = rgi
                x_test = rgi
                y_train_res = ast[abx]
                y_test = ast[abx]
                abx = abx + "*"
            # otherwise take it and use SMOTE to rebalance the training set
            else:
                # 20% totally held out test-set
                x_train, x_test, y_train, y_test = model_selection.train_test_split(rgi, ast[abx],
                                                      test_size = 0.2, stratify=ast[abx],
                                                      random_state=42)
                # resampled to even class balance via SMOTE
                sm = over_sampling.SMOTE(random_state=42, ratio = 1.0)
                x_train_res, y_train_res = sm.fit_sample(x_train, y_train)
        
            # use SSS to fit best interpreter per model
            sss = model_selection.StratifiedShuffleSplit(n_splits=3, random_state=42)
    
            param_search = model_selection.GridSearchCV(clf, cv = sss,
                                                        param_grid=params[clf_name],
                                                        scoring='neg_log_loss')
            try:
                param_search.fit(x_train_res, y_train_res)
                  
                if clf_name not in classifier_perf:
                    classifier_perf.update({clf_name: {abx: (param_search.best_params_, 
                                                                param_search.best_score_)}})
                else:
                    classifier_perf[clf_name].update({abx: (param_search.best_params_,
                                                            param_search.best_score_)})  
            except ValueError:
                print(abx, 'failed')
                clf.fit(x_train_res, y_train_res)
                y_pred = clf.predict(x_train_res)
                
                score = metrics.log_loss(y_train_res, y_pred)
                
                if clf_name not in classifier_perf:
                    classifier_perf.update({clf_name: {abx: ('failed', 
                                                            score)}})
                else:
                    classifier_perf[clf_name].update({abx: ('failed',
                                                            score)})  
    
          
    plot_classifier_cv_perf(classifier_perf)


def expand_fields(df):
    dfs = []
    for row in df:
        temp_df = df[row].apply(pd.Series)
        temp_df = temp_df.rename(columns={0: row + '_params',
                                          1: row + '_score'})
        dfs.append(temp_df)
    return pd.concat(dfs, axis=1)
    

def plot_classifier_cv_perf(classifier_perf):
  
    df = pd.DataFrame(classifier_perf)
    df = expand_fields(df)

    score_df = df[['DT_score', 'LR_score', 'NB_score', 'RF_score']]
    score_df = pd.melt(score_df.reset_index(), id_vars='index')
    sns.violinplot(data = score_df, y='value', x='variable')
    plt.title('Per Antibiotic Binary Accuracy Distribution (neg log loss)')


# In[4]:


classifiers = {'LR': LogisticRegression(class_weight='balanced'),
                   'DT': DecisionTreeClassifier(class_weight='balanced'),
                   'RF': RandomForestClassifier(class_weight='balanced'),
                   'NB': MultinomialNB()}
    
params = {'LR': {'penalty': ['l1', 'l2'],
                     'C': [0.2, 0.5, 1, 1.5]},
              'DT': {"max_depth": [5, 3, None],
                     "max_features": [1, 5, 11],
                     "min_samples_split": [2, 5, 11],
                     "min_samples_leaf": [1, 5, 11],
                     "criterion": ["gini", "entropy"]},
              'RF': {"max_depth": [5, 3, None],
                     "max_features": [1, 5, 11],
                 "min_samples_split": [2, 5, 11],
                     "min_samples_leaf": [1, 5, 11],
                     "bootstrap": [True, False],
                     "criterion": ["gini", "entropy"]},
              'NB': {"alpha": [0.3, 0.7, 1],
                    'fit_prior': [True, False]}}

fit_interpretable_clf_per_abx(ecoli_rgi, ecoli_ast, classifiers, params)


# As logistic regression and random forest were largely equivalent in performance

# # Fit and Evaluate E. coli Models

# In[7]:


def fit_LR_model_per_abx(rgi, ast):
    """
    Fit logistic regression
    """
    models = {}
    models_av_prec = {}
    models_auc = {}
    
    for abx in ast.columns:        
        #filter rgi to just genomes with AST due to NAs
        abx_rgi = rgi[~ast[abx].isna()]
        abx_ast = ast.loc[~ast[abx].isna(), abx]
        
        props = abx_ast.value_counts()
        label_prop = props.min() / props.max() 
        
        # if all one class (and in our data this is only all resistant) then just return 
        # resistant for everything but add two asterisks
        if props.max() == abx_ast.shape[0]:
            abx = abx + "**"
            models.update({abx: DummyClassifier(strategy='constant', constant=1).fit(abx_rgi, abx_ast)})
            models_av_prec.update({abx: 1})
            models_auc.update({abx: 1})
            continue
        
        # check balance if less than 10% are minority class then just train all
        if label_prop < 0.1:
            x_train_res = abx_rgi
            x_test = abx_rgi
            y_train_res = abx_ast
            y_test = abx_ast
            abx = abx + "*"
        else:
            # if balanced in any shape then have 20% totally held out test-set
            x_train, x_test, y_train, y_test = model_selection.train_test_split(abx_rgi, abx_ast,
                                                      test_size = 0.2, stratify=abx_ast,
                                                      random_state=42)
            # resampled the training set to even class balance via SMOTE
            sm = over_sampling.SMOTE(random_state=42, ratio = 1.0)
            x_train_res, y_train_res = sm.fit_sample(x_train, y_train)
        
        # train and tune an LR model using 3-fold CV
        try:
            clf = LogisticRegressionCV(cv=3, solver='liblinear', class_weight='balanced')

            clf.fit(x_train_res, y_train_res)
            
        # if this fails as there are too few labels of one class to use CV then just 
        # (over)fit the sklearn default LR
        except ValueError as e:
            if str(e).startswith("Class label "):
                clf = LogisticRegression()
                clf.fit(x_train_res, y_train_res)
            else:
                raise
        
        # generate prediction metrics
        y_pred = clf.predict_proba(x_test)[:,1]
    
        prec, recall, thresholds = metrics.precision_recall_curve(y_test, y_pred)
        av_prec = metrics.average_precision_score(y_test, y_pred)
    
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        auc = metrics.roc_auc_score(y_test, y_pred)
    
        f, axarr = plt.subplots(1, 2)
        plot_prec_recall(prec, recall, av_prec, abx, axarr[0])
        plot_roc(fpr, tpr, auc, abx, axarr[1])
        plt.tight_layout()
        
        models.update({abx: clf})
        models_av_prec.update({abx: av_prec})
        models_auc.update({abx: auc})
    return models, models_av_prec, models_auc
    
def plot_prec_recall(precision, recall, av_prec, abx, ax):
    ax.step(recall, precision, color='b', alpha=0.2,
             where='post')
    ax.fill_between(recall, precision, step='post', alpha=0.2,
                 color='b')

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_ylim([0.0, 1.05])
    ax.set_xlim([0.0, 1.0])
    ax.set_title('{0} PR AP={1:0.2f}'.format(abx, 
          av_prec))
    
def plot_roc(fpr, tpr, auc, abx, ax):    
    ax.plot(fpr, tpr, color='darkorange', lw=2)
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('{0} ROC AUC={1:0.2f})'.format(abx, auc))


# In[8]:


#Get metrics for tuned LR models for ecoli
ecoli_models, ecoli_av_prec, ecoli_auc = fit_LR_model_per_abx(ecoli_rgi, ecoli_ast)
plt.show()

print(pd.Series(ecoli_auc).describe())
pd.Series(ecoli_auc).sort_values().plot(kind='bar')
plt.title('E. coli Area Under ROC')
plt.show()

pd.Series(ecoli_av_prec).sort_values().plot(kind='bar')
print(pd.Series(ecoli_av_prec).describe())
plt.title('E. coli Average Precision')
plt.show()


# In[9]:


# feature importances
def plot_feature_importances(classifiers, rgi):
    feature_importance = pd.DataFrame()
    for abx in classifiers:
        try:
            feature_importance[abx] =  pd.Series(classifiers[abx].coef_[0], name=abx)
        except AttributeError:
            continue
        
    feature_importance.index = rgi.columns
    for ix, abx in enumerate(feature_importance):
        top20 = feature_importance[abx].nlargest(20)
        plt.bar( top20.index, top20.values)
        plt.title(abx)
        plt.xticks(rotation=90)
        plt.ylabel('Weight')
        # plt.show()


# In[10]:


plot_feature_importances(ecoli_models, ecoli_rgi)


# # Fit and Evaluate Pseudomonas Models

# In[11]:


#Get metrics for tuned LR models for pseudo
pseudo_models, pseudo_av_prec, psuedo_auc = fit_LR_model_per_abx(pseudo_rgi, pseudo_ast)
plt.show()

print(pd.Series(psuedo_auc).describe())
pd.Series(psuedo_auc).sort_values().plot(kind='bar')
plt.title('P. aeruginosa Area Under ROC')
plt.show()

pd.Series(pseudo_av_prec).sort_values().plot(kind='bar')
print(pd.Series(pseudo_av_prec).describe())
plt.title('P. aeruginosa Average Precision')
plt.show()


# In[12]:


plot_feature_importances(pseudo_models, pseudo_rgi)


# # Generate Main Metrics

# In[13]:


def over_under_pred_tally(ast, rgi, classifier_dict):
    
    data = {'antibiotic': [],
            'classifier': [],
            'Over Prediction': [],
            'Under Prediction': [],
            'Complete Prediction': []}
    
    # reorder to match paper plot
    try:
        ast = ast[['ampicillin',
              'amoxicillin-clavulanic_acid',
              'piperacillin-tazobactam',
              'cefazolin',
              'cefalotin',
              'ceftriaxone',
              'ceftazidime',
              'cefixime',
              'cefoxitin',
              'ertapenem',
              'meropenem',
              'nitrofurantoin',
              'tetracycline',
              'trimethoprim-sulfamethoxazole',
              'ciprofloxacin',
              'gentamicin_C',
              'amikacin',
              'tobramycin']]
        add_ert = False
        
    # deal with ertapenem not being done for pseudomonas
    except KeyError:
        add_ert = True
        ast = ast[['ampicillin',
              'amoxicillin-clavulanic_acid',
              'piperacillin-tazobactam',
              'cefazolin',
              'cefalotin',
              'ceftriaxone',
              'ceftazidime',
              'cefixime',
              'cefoxitin',
              'meropenem',
              'nitrofurantoin',
              'tetracycline',
              'trimethoprim-sulfamethoxazole',
              'ciprofloxacin',
              'gentamicin_C',
              'amikacin',
              'tobramycin']]

    
    # rename AST sequences to noted names
    for clf_abx in classifier_dict:
        for abx_ast in ast:
            if abx_ast.startswith(clf_abx[:-3]):
                ast = ast.rename(columns={abx_ast:clf_abx})


    for abx in ast:
        abx_rgi = rgi[~ast[abx].isna()]
        abx_ast = ast.loc[~ast[abx].isna(), abx]

        ast_pred = classifier_dict[abx].predict(abx_rgi)
        
        complete, under, over = 0, 0, 0
        for ix, ast_pred in enumerate(ast_pred):
            if ast_pred == 1 and abx_ast[ix] == 1:
                complete += 1
            elif ast_pred == 0 and abx_ast[ix] == 1:
                under += 1
            elif ast_pred == 1 and abx_ast[ix] == 0:
                over += 1
        
        data['antibiotic'].append(abx.capitalize())
        data['classifier'].append(classifier_dict[abx])
        data['Over Prediction'].append(over)
        data['Under Prediction'].append(under)
        data['Complete Prediction'].append(complete)
        
        if add_ert and abx.startswith('cefoxitin'):
            data['antibiotic'].append('Ertapenem')
            data['classifier'].append([])
            data['Over Prediction'].append(0)
            data['Under Prediction'].append(0)
            data['Complete Prediction'].append(0)      
    
    return pd.DataFrame(data)


# In[14]:


ecoli_perf = over_under_pred_tally(ecoli_ast, ecoli_rgi, ecoli_models)
complete='#822340'
over='#F5B233'
grey='#797979'
ecoli_perf.plot(kind='barh', stacked=True, colors=[complete, over, grey])
_ = plt.yticks(range(18), ecoli_perf['antibiotic'])
plt.legend(loc=0)
plt.xlabel('Number of E. coli isolates (n=115)')
plt.tight_layout()
plt.savefig('ecoli.svg')
plt.savefig('ecoli.png')
plt.savefig('ecoli.pdf')


# In[15]:


pseudo_perf = over_under_pred_tally(pseudo_ast, pseudo_rgi, pseudo_models)
complete='#822340'
over='#F5B233'
grey='#797979'
pseudo_perf.plot(kind='barh', stacked=True, colors=[complete, over, grey])
_ = plt.yticks(range(18), pseudo_perf['antibiotic'])
plt.xlabel('Number of P. aeruginosa isolates (n=102)')
plt.tight_layout()
plt.savefig('pseudo.svg')
plt.savefig('pseudo.png')
plt.savefig('pseudo.pdf')

