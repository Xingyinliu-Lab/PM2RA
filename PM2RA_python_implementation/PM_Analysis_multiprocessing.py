import pandas as pd
from scipy.integrate import quad
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import numpy as np
from scipy.stats import ks_2samp
import datetime
import fdrcorrection
import load_save_project
import sys
import copy
np.seterr(divide='ignore',invalid='ignore')
fileplace=sys.argv[1]
logs_place=sys.argv[2]
datafilename=sys.argv[3]
minimum_coexist_taxa_number=int(sys.argv[4])
minimum_taxa_number=int(sys.argv[5])
minimum_taxa_prevlance=float(sys.argv[6])
control=sys.argv[7]
treat=sys.argv[8]
condition=sys.argv[9]
processers=int(sys.argv[10])
fdr=sys.argv[11]
ilr=sys.argv[12]
confidence_level=float(sys.argv[13])
wetherFindBestBandwidth=sys.argv[14]
kernel_str=sys.argv[15]
#
# fileplace='D:\\project\\exe\\pm\\Project\\12\\'
# logs_place='D:\\project\\exe\\pm\\Project\\12\\logs\\'
# datafilename='D:\\project\\exe\\pm\\demodata_small.csv'
# minimum_coexist_taxa_number=10
# minimum_taxa_number=10
# minimum_taxa_prevlance=0.000
# control= 'H2029'
# treat = 'crc'
# condition= 'condition'
# processers=6
# fdr='ON'
# ilr='ON'
# confidence_level=0.95
# wetherFindBestBandwidth='ON'
# kernel_str='gaussian'

#
# fileplace='D:\\project\\exe\\pm\\Project\\12\\'
# logs_place='D:\\project\\exe\\pm\\Project\\12\\logs\\'
# datafilename='D:\\project\\exe\\pm\\demoASD2-3.csv'
# minimum_coexist_taxa_number=10
# minimum_taxa_number=10
# minimum_taxa_prevlance=0.0000
# control= 'A2-3'
# treat = 'B2-3'
# condition= 'Age_Group'
# processers=6
# fdr='ON'
# ilr='ON'
# confidence_level=0.95
# wetherFindBestBandwidth='ON'
# kernel_str='gaussian'


demodata = pd.read_csv(datafilename, header=0, index_col=False)
taxaname = demodata.columns
# convert absolute abundance to relative abundance
genuslist = []
for g in taxaname:
    if g != condition:
        genuslist.append(g)
sumabundance = demodata[genuslist].sum(axis=1)
demodata[genuslist] = demodata[genuslist].div(sumabundance, axis='rows')
# ilr

tmpd=demodata[genuslist]
tmpd[tmpd<=0]=np.nan
demodata[genuslist] =tmpd
orireldata=copy.copy(demodata)

if ilr=='ON':
    demodata[genuslist] = np.log(demodata[genuslist])
    gmean = demodata[genuslist].mean(axis=1)
    demodata[genuslist] = demodata[genuslist].subtract(gmean, axis='rows')
# filter low detection taxa
genuslist = []
for g in taxaname:
    if g != condition:
        controlgdata = orireldata.loc[orireldata[condition]
                                     == control, g]
        treatgdata = orireldata.loc[orireldata[condition] == treat, g]
        controlobs = sum(controlgdata > 0)
        treatobs = sum(treatgdata > 0)
        median_control_abundance=np.nanmedian(controlgdata)
        if controlobs>minimum_taxa_number and treatobs>minimum_taxa_number and median_control_abundance>minimum_taxa_prevlance:
            genuslist.append(g)
demodata=demodata[genuslist+list([condition])]
taxaname = demodata.columns


def hotellingt2d(data,data_average,data_cov):
    data_cov = np.matrix(data_cov)
    try:
        data_cov_r=np.linalg.inv(data_cov)
    except:
        data_cov_r = np.linalg.pinv(data_cov)
    cdata=data-data_average.reshape(-1,1)
    t2d=np.diag(np.matmul(np.matmul(np.transpose(cdata), data_cov_r),cdata))
    return t2d

# kernels = ['gaussian', 'tophat', 'epanechnikov']

def getKDE(x):
    kde=None
    if wetherFindBestBandwidth:
        grid = GridSearchCV(KernelDensity(kernel=kernel_str),
                            {'bandwidth': np.linspace(0.05,1, 20)},
                            cv=5)  # 20-fold cross-validation
        grid.fit(x[:, None])
        kde = grid.best_estimator_
    if not wetherFindBestBandwidth:
        kde = KernelDensity(kernel=kernel_str, bandwidth=0.1).fit(x[:, None])
    return kde

def inte_ked(low,high,e_kde):
    def f(x,e_kde):
        return np.exp(e_kde.score_samples(np.array([x]).reshape(-1, 1)))[0]
    v, err = quad(f, low, high, args = (e_kde))
    return v

def inte_2ked(low,high,e_kde1,e_kde2):
    def f(x,e_kde1,e_kde2):
        e1=np.exp(e_kde1.score_samples(np.array([x]).reshape(-1, 1)))[0]
        e2=np.exp(e_kde2.score_samples(np.array([x]).reshape(-1, 1)))[0]
        return np.min([e1,e2])
    v, err = quad(f, low, high, args = (e_kde1,e_kde2),limit=200,)
    return v

def remove_outlier(x):
    Q1 = np.quantile(x, 0.25)
    Q2 = np.quantile(x, 0.50)
    Q3 = np.quantile(x, 0.75)
    upper=Q2+1.5*(Q3-Q1)
    return x[x<upper]


def cal_pm_socre_nD(taxalist,demodata):
    tmpdemodata=demodata[taxalist+[condition]]
    tmpdemodata=tmpdemodata.dropna()
    tmpdemodata=tmpdemodata.reset_index(drop=True)
    control_obs=tmpdemodata[tmpdemodata[condition]==control]
    control_obs=np.transpose(np.asarray(control_obs[taxalist]))
    treat_obs=tmpdemodata[tmpdemodata[condition]==treat]
    treat_obs = np.transpose(np.asarray(treat_obs[taxalist]))
    if (control_obs.shape[1]<minimum_coexist_taxa_number) or (treat_obs.shape[1]<minimum_coexist_taxa_number):
        return None,None,None,None,None,None,None,None,None,None
    control_average=np.nanmean(control_obs,axis=1)
    control_cov=np.cov(control_obs,rowvar=True)
    treat_average=np.average(treat_obs,axis=1)
    treat_cov=np.cov(treat_obs,rowvar=True)
    cc=hotellingt2d(control_obs,control_average,control_cov)
    ct=hotellingt2d(control_obs,treat_average,treat_cov)
    tt=hotellingt2d(treat_obs,treat_average,treat_cov)
    tc=hotellingt2d(treat_obs,control_average,control_cov)
    cc=remove_outlier(cc)
    ct=remove_outlier(ct)
    tt=remove_outlier(tt)
    tc=remove_outlier(tc)
    _,cc_tc_pvalue = ks_2samp(cc,tc)
    _, tt_ct_pvalue = ks_2samp(tt, ct)
    # if (cc_tc_pvalue>1-confidence_level) and (tt_ct_pvalue>1-confidence_level) :
    #     return 0,np.min([cc_tc_pvalue,tt_ct_pvalue]),None,None,None,None,None,None,None,None
    ## cc and tc
    cc_kde, ct_kde, tt_kde, tc_kde=None,None,None,None
    if (np.min(cc)>np.max(tc)) or (np.min(tc)>np.max(cc)):
        cc_tc_pm=1
    else:
        cc_tc_lowerbound=np.max([np.min(cc),np.min(tc)])
        cc_tc_upperbound=np.min([np.max(cc),np.max(tc)])
        cc_kde = getKDE(cc)
        tc_kde = getKDE(tc)
        cc_tc_pm=inte_2ked(cc_tc_lowerbound,cc_tc_upperbound,cc_kde,tc_kde)
    ## ct and tt
    if (np.min(tt)>np.max(ct)) or (np.min(ct)>np.max(tt)):
        tt_ct_pm=1
    else:
        tt_ct_lowerbound=np.max([np.min(tt),np.min(ct)])
        tt_ct_upperbound=np.min([np.max(tt),np.max(ct)])
        ct_kde = getKDE(ct)
        tt_kde = getKDE(tt)
        tt_ct_pm=inte_2ked(tt_ct_lowerbound,tt_ct_upperbound,ct_kde,tt_kde)
    pm_score=np.max([1-tt_ct_pm,1-cc_tc_pm])
    return pm_score,np.min([cc_tc_pvalue,tt_ct_pvalue]),cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde




def cal_1D(comb,i):
    res = pd.DataFrame(columns=['taxa','pm', 'pvalue'])
    run_log=pd.DataFrame(columns=['taxa','numpyindex'])
    cc_list=[]
    ct_list=[]
    tt_list=[]
    tc_list=[]
    cc_kde_list=[]
    ct_kde_list=[]
    tt_kde_list=[]
    tc_kde_list=[]
    count = 0
    for m in comb:
        try:
            m=list(m)
            pm_score, pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde = cal_pm_socre_nD(m, demodata)
            if (pm_score is None) or (cc_kde is None) or (ct_kde is None) or  (tt_kde is None) or (tc_kde is None):
                continue
            res.loc[count, 'taxa'] = m[0]
            res.loc[count, 'pm'] = pm_score
            res.loc[count, 'pvalue'] = pvalue
            run_log.loc[count,'taxa']= m[0]
            run_log.loc[count, 'numpyindex'] = count
            count = count + 1
            cc_list.append(cc)
            ct_list.append(ct)
            tt_list.append(tt)
            tc_list.append(tc)
            cc_kde_list.append(cc_kde)
            ct_kde_list.append(ct_kde)
            tt_kde_list.append(tt_kde)
            tc_kde_list.append(tc_kde)
            if count%3==0:
                finishedtask = {'hasfinished': count,
                                'isfinished':False,
                             }
                load_save_project.save_dict(logs_place + '/'+str(i)+'_task_1D.json', finishedtask)
        except:
            pass
    finishedtask = {'hasfinished': count,
                    'isfinished': True,
                    }
    load_save_project.save_dict(logs_place + '/' + str(i) + '_task_1D.json', finishedtask)
    np.save(logs_place + '/' + str(i) + '_cc_list_1D.npy',cc_list)
    np.save(logs_place + '/' + str(i) + '_ct_list_1D.npy', ct_list)
    np.save(logs_place + '/' + str(i) + '_tt_list_1D.npy', tt_list)
    np.save(logs_place + '/' + str(i) + '_tc_list_1D.npy', tc_list)
    np.save(logs_place + '/' + str(i) + '_cc_kde_list_1D.npy', cc_kde_list)
    np.save(logs_place + '/' + str(i) + '_ct_kde_list_1D.npy', ct_kde_list)
    np.save(logs_place + '/' + str(i) + '_tt_kde_list_1D.npy', tt_kde_list)
    np.save(logs_place + '/' + str(i) + '_tc_kde_list_1D.npy', tc_kde_list)
    run_log.to_csv(logs_place + '/' + str(i) + 'runlog_1D.csv')
    return res

def cal_2D(comb,i):
    res = pd.DataFrame(columns=['taxa1', 'taxa2', 'pm', 'pvalue'])
    run_log=pd.DataFrame(columns=['taxa1', 'taxa2', 'numpyindex'])
    cc_list=[]
    ct_list=[]
    tt_list=[]
    tc_list=[]
    cc_kde_list=[]
    ct_kde_list=[]
    tt_kde_list=[]
    tc_kde_list=[]
    count = 0
    for m in comb:
        try:
            m=list(m)
            pm_score, pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde = cal_pm_socre_nD(m, demodata)
            if (pm_score is None) or (cc_kde is None) or (ct_kde is None) or  (tt_kde is None) or (tc_kde is None):
                continue
            res.loc[count, 'taxa1'] = m[0]
            res.loc[count, 'taxa2'] = m[1]
            res.loc[count, 'pm'] = pm_score
            res.loc[count, 'pvalue'] = pvalue
            run_log.loc[count,'taxa1']= m[0]
            run_log.loc[count, 'taxa2'] = m[1]
            run_log.loc[count, 'numpyindex'] = count
            cc_list.append(cc)
            ct_list.append(ct)
            tt_list.append(tt)
            tc_list.append(tc)
            cc_kde_list.append(cc_kde)
            ct_kde_list.append(ct_kde)
            tt_kde_list.append(tt_kde)
            tc_kde_list.append(tc_kde)
            count = count + 1
            if count%3==0:
                finishedtask = {'hasfinished': count,
                                'isfinished':False,
                             }
                load_save_project.save_dict(logs_place + '/'+str(i)+'_task.json', finishedtask)
        except:
            pass
    finishedtask = {'hasfinished': count,
                    'isfinished': True,
                    }
    load_save_project.save_dict(logs_place + '/' + str(i) + '_task.json', finishedtask)
    np.save(logs_place + '/' + str(i) + '_cc_list.npy',cc_list)
    np.save(logs_place + '/' + str(i) + '_ct_list.npy', ct_list)
    np.save(logs_place + '/' + str(i) + '_tt_list.npy', tt_list)
    np.save(logs_place + '/' + str(i) + '_tc_list.npy', tc_list)
    np.save(logs_place + '/' + str(i) + '_cc_kde_list.npy', cc_kde_list)
    np.save(logs_place + '/' + str(i) + '_ct_kde_list.npy', ct_kde_list)
    np.save(logs_place + '/' + str(i) + '_tt_kde_list.npy', tt_kde_list)
    np.save(logs_place + '/' + str(i) + '_tc_kde_list.npy', tc_kde_list)
    run_log.to_csv(logs_place + '/' + str(i) + 'runlog.csv')
    return res

from multiprocessing import Pool
if __name__ == '__main__':
    success_1D = False
    pm_res=None
    try:
        comb = []
        for a in range(len(genuslist)):
            if genuslist[a] == condition:
                continue
            comb.append([genuslist[a]])
        totaltask = {'totalcomb': len(comb),
                     'isfinished': False
                     }
        load_save_project.save_dict(logs_place + '/totaltask_1D.json', totaltask)
        # print(datetime.datetime.now())
        cl = np.array_split(np.asarray(range(len(comb))), processers, axis=0)
        comb = np.asarray(comb)
        res_list = []
        p = Pool(processers)
        for i in range(processers):
            res_list.append(p.apply_async(cal_1D, args=(comb[cl[i]], i,)))
        p.close()
        p.join()
        pm_res = pd.concat([i.get() for i in res_list])
        tmpPvalue = pm_res['pvalue']
        # fdrcorrection:
        # This covers Benjamini/Hochberg for independent or positively correlated and Benjamini/Yekutieli
        # for general or negatively correlated tests. Both are available in the function multipletests, as method=`fdr_bh`, resp. fdr_by.
        if fdr == 'ON':
            _, qvaluelist = fdrcorrection.fdrcorrection_np(tmpPvalue, alpha=0.05, method='indep', is_sorted=False)
            pm_res['qvalue'] = qvaluelist
        res_1D=pm_res
        pm_res.to_csv(fileplace + '/PM_scores_1D.csv')
        # print(datetime.datetime.now())
        success_1D = True
    except:
        pass
    totaltask = {'success': success_1D
                 }
    load_save_project.save_dict(logs_place + '/task_finished_1D.json', totaltask)

    success=False
    try:
        comb=[]
        for a in range(len(genuslist)):
            for b in range(len(genuslist)):
                if a>=b:
                    continue
                if genuslist[a]==condition:
                    continue
                if genuslist[b]==condition:
                    continue
                tmpdemodata=demodata[[genuslist[a],genuslist[b],condition]]
                tmpdemodata=tmpdemodata.dropna()
                tmpdemodata=tmpdemodata.reset_index(drop=True)
                control_obs=tmpdemodata[tmpdemodata[condition]==control]
                control_obs=np.transpose(np.asarray(control_obs[[genuslist[a],genuslist[b]]]))
                treat_obs=tmpdemodata[tmpdemodata[condition]==treat]
                treat_obs = np.transpose(np.asarray(treat_obs[[genuslist[a],genuslist[b]]]))
                if (control_obs.shape[1]<minimum_coexist_taxa_number) or (treat_obs.shape[1]<minimum_coexist_taxa_number):
                    continue
                comb.append([genuslist[a],genuslist[b]])
        totaltask = {'totalcomb': len(comb),
                     'isfinished': False
                        }
        load_save_project.save_dict(logs_place + '/totaltask.json', totaltask)

        # print(datetime.datetime.now())
        cl = np.array_split(np.asarray(range(len(comb))), processers, axis=0)
        comb=np.asarray(comb)
        res_list = []
        p = Pool(processers)
        for i in range(processers):
            res_list.append(p.apply_async(cal_2D, args=(comb[cl[i]],i,)))
        p.close()
        p.join()
        pm_res = pd.concat([i.get() for i in res_list])
        tmpPvalue = pm_res['pvalue']
        # fdrcorrection:
        # This covers Benjamini/Hochberg for independent or positively correlated and Benjamini/Yekutieli
        # for general or negatively correlated tests. Both are available in the function multipletests, as method=`fdr_bh`, resp. fdr_by.
        if fdr=='ON':
            _, qvaluelist = fdrcorrection.fdrcorrection_np(tmpPvalue, alpha=0.05, method='indep', is_sorted=False)
            pm_res['qvalue']=qvaluelist
        if fdr=='ON':
            tmpres1D1 = res_1D[['taxa', 'pm','pvalue','qvalue']]
            tmpres1D2 = res_1D[['taxa', 'pm','pvalue','qvalue']]
            tmpres1D1 = tmpres1D1.rename(columns={'taxa': 'taxa1', 'pm': 'pm1','pvalue':'pvalue1','qvalue':'qvalue1'})
            tmpres1D2 = tmpres1D2.rename(columns={'taxa': 'taxa2', 'pm': 'pm2','pvalue':'pvalue2','qvalue':'qvalue2'})
        if fdr=='OFF':
            tmpres1D1 = res_1D[['taxa', 'pm','pvalue']]
            tmpres1D2 = res_1D[['taxa', 'pm','pvalue']]
            tmpres1D1 = tmpres1D1.rename(columns={'taxa': 'taxa1', 'pm': 'pm1','pvalue':'pvalue1'})
            tmpres1D2 = tmpres1D2.rename(columns={'taxa': 'taxa2', 'pm': 'pm2','pvalue':'pvalue2'})
        res = pd.merge(left=pm_res, right=tmpres1D1, how='left', on='taxa1')
        res = pd.merge(left=res, right=tmpres1D2, how='left', on='taxa2')
        res=res.rename(columns={'pm':'raw_pm_2d'})
        def co_PM_2D_fun(raw_pm_2d, pm1, pm2):
            return np.max([0, np.min([raw_pm_2d - pm1, raw_pm_2d - pm2])])
        res['co_PM_2D'] = res.apply(lambda x: co_PM_2D_fun(x.raw_pm_2d, x.pm1, x.pm2), axis=1)
        res.to_csv(fileplace+'/PM_scores.csv')
        # print(datetime.datetime.now())
        success=True
    except:
        pass
    totaltask = {'success': success
                    }
    load_save_project.save_dict(logs_place + '/task_finished.json', totaltask)
