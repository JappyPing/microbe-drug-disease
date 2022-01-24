'''
Author: Pengyao PING
Date: 2022-01-24 23:57:17
LastEditors: Pengyao PING
LastEditTime: 2022-01-25 02:59:31
Email: Pengyao.Ping@student.uts.edu.au
Description: 
'''
import scipy.io
from sklearn.metrics import precision_recall_curve, roc_curve,auc,f1_score
import numpy as np
import pandas as pd

def auc_aupr(i,type):
    s="./ten_fold_predict_result/"+str(type)+"/Predict_result"+str(i)
    p = scipy.io.loadmat(s)
    if str(type)=='BNNR':
        p = p['F1']
    if str(type)=='DR2DI':
        p = p['F2']
    if str(type) == 'MBiRW':
        p = p['F3']
    if str(type) == 'MSBMF':
        p = p['F4']
    dis_drug = pd.read_csv("./Data/label.csv")
    dis_drug = dis_drug.values
    label = []
    result = []
    for i in range(np.size(dis_drug, axis=0)):
        for j in range(np.size(dis_drug, axis=1)):
            label.append(dis_drug[i][j])
            result.append(p[i][j])
    fpr, tpr, thr = roc_curve(label, result)
    auc_val = auc(fpr, tpr)
    precision,recall, thresholds=precision_recall_curve(label, result)
    aupr_val = auc(recall, precision)

    #f1
    paddr=np.add(precision,recall)              #p+r
    pmultiplyr=np.multiply(precision,recall)    #p*r
    f1score=np.divide(2*pmultiplyr,paddr)
    return auc_val,aupr_val,f1score,fpr,tpr,precision,recall,thresholds

def auc_aupr_ten(type):
    auc_val = 0
    aupr_val = 0
    f1_val = []
    fpr_val = []
    tpr_val = []
    precision_val = []
    recall_val = []
    thresholds_val = []
    for j in range(1,11):
        auc1, aupr1, f1, fpr, tpr, precision, recall, thresholds = auc_aupr(j,type)
        auc_val = auc_val + auc1
        aupr_val = aupr_val + aupr1
        f1_val.append(f1)
        fpr_val.append(fpr)
        tpr_val.append(tpr)
        precision_val.append(precision)
        recall_val.append(recall)
        thresholds_val.append(thresholds)
    #scipy.io.savemat('./Evaluate_data/'+str(type)+'/f1score_'+str(type)+'.mat', mdict={'f1_val': f1_val})
    # scipy.io.savemat('./Evaluate_data/'+str(type)+'/fpr_'+str(type)+'.mat', mdict={'fpr_val': fpr_val})
    # scipy.io.savemat('./Evaluate_data/'+str(type)+'/tpr_'+str(type)+'.mat', mdict={'tpr_val': tpr_val})
    #scipy.io.savemat('./Evaluate_data/'+str(type)+'/precision_'+str(type)+'.mat', mdict={'precision_val': precision_val})
    #scipy.io.savemat('./Evaluate_data/'+str(type)+'/recall_'+str(type)+'.mat', mdict={'recall_val': recall_val})
    #scipy.io.savemat('./Evaluate_data/'+str(type)+'/f1_thresholds_'+str(type)+'.mat', mdict={'thresholds_val': thresholds_val})
    print(str(type)+"_avg_auc:",auc_val/10)
    print(str(type) + "_avg_aupr:", aupr_val/10)