import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import copy
from matplotlib.patches import Ellipse


# resfilename='D:\\project\\exe\\pm\\Project\\ASD2_3\\PM_scores.csv'
# pm_res=pd.read_csv(resfilename,header=0,index_col=0)
# pm_res=pm_res[pm_res['co_PM_2D']>0]
# pm_res=pm_res.sort_values(by='co_PM_2D',ascending=False)
# pm_res=pm_res.reset_index(drop=True)
# pm_res=pm_res[0:10]
# datafilename='D:\\project\\exe\\pm\\demoASD2-3.csv'
# control_label='A2-3'
# treat_label='B2-3'
# condition='Age_Group'
# threads_count=6
# logs_place='D:\\project\\exe\\pm\\Project\\ASD2_3\\logs'
# fileplace='D:\\project\\exe\\pm\\Project\\ASD2_3'
# image_pdf=fileplace+'/taxa_relationship.pdf'





def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Groups')



def plot_scatter_exe(project_dict,topnum,networktype,filter_ss):
    fileplace = project_dict.get('fileplace')
    fdr = project_dict.get('fdr')
    confidenceI= [0.99, 0.95,0.90][project_dict.get('confidenceI')]
    resfilename = fileplace+'/PM_scores.csv'
    pm_res = pd.read_csv(resfilename, header=0, index_col=0)
    def co_PM_2D_fun(raw_pm_2d,pm1,pm2):
        return np.max([0,np.min([raw_pm_2d-pm1,raw_pm_2d-pm2])])

    if filter_ss == 'ON':
        if fdr == 'ON':
            pm_res.loc[pm_res['qvalue']>1-confidenceI,'raw_pm_2d']=0
            pm_res.loc[pm_res['qvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['qvalue2'] > 1 - confidenceI, 'pm2'] = 0
        if fdr == 'OFF':
            pm_res.loc[pm_res['pvalue'] > 1 - confidenceI, 'raw_pm_2d'] = 0
            pm_res.loc[pm_res['pvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['pvalue2'] > 1 - confidenceI, 'pm2'] = 0

    pm_res['co_PM_2D']=pm_res.apply(lambda x: co_PM_2D_fun(x.raw_pm_2d, x.pm1,x.pm2), axis = 1)
    target = 'raw_pm_2d'
    if networktype == 'Interaction PM scores':
        target = 'co_PM_2D'
    if networktype=='2D PM Scores':
        target = 'raw_pm_2d'
    pm_res=pm_res[pm_res[target]>0]
    pm_res=pm_res.sort_values(by=target,ascending=False)
    pm_res=pm_res.reset_index(drop=True)
    taxnum = np.min([len(pm_res), topnum])
    pm_res = pm_res[0:taxnum]

    datafilename = project_dict.get('datafilename')
    control_label =project_dict.get('controllabel')
    treat_label = project_dict.get('treatlabel')
    condition = project_dict.get('grouplabel')
    multithreading= project_dict.get('multithreading')
    if multithreading=='ON':
        threads_count = int(project_dict.get('threads'))
    if multithreading=='OFF':
        threads_count = 1
    logs_place =project_dict.get('logs_place')
    image_pdf = fileplace + '/taxa_relationship_difference_between_groups.pdf'
    demodata = pd.read_csv(datafilename, header=0, index_col=False)
    taxaname = demodata.columns
    # convert absolute abundance to relative abundance
    genuslist = []
    for g in taxaname:
        if g != condition:
            genuslist.append(g)
    sumabundance = demodata[genuslist].sum(axis=1)
    demodata[genuslist] = demodata[genuslist].div(sumabundance, axis='rows')
    tmpd = demodata[genuslist]
    tmpd[tmpd <= 0] = np.nan
    demodata[genuslist] = tmpd
    ilr_data = copy.copy(demodata)
    ilr_data[genuslist] = np.log(ilr_data[genuslist])
    gmean = ilr_data[genuslist].mean(axis=1)
    ilr_data[genuslist] = ilr_data[genuslist].subtract(gmean, axis='rows')
    def plot_scatter():
        with PdfPages(image_pdf) as pdf:
            for pmi in pm_res.index:
                taxa1=pm_res.loc[pmi,'taxa1']
                taxa2 = pm_res.loc[pmi, 'taxa2']
                tmpdemodata = demodata[[taxa1, taxa2, condition]]
                tmpdemodata = tmpdemodata.dropna()
                tmpdemodata = tmpdemodata.reset_index(drop=True)
                control_a = tmpdemodata.loc[tmpdemodata[condition] == control_label, taxa1]
                control_b = tmpdemodata.loc[tmpdemodata[condition] == control_label, taxa2]
                treat_a = tmpdemodata.loc[tmpdemodata[condition] == treat_label, taxa1]
                treat_b = tmpdemodata.loc[tmpdemodata[condition] == treat_label, taxa2]
                fig= plt.figure(figsize=(20,10.5))
                plt.suptitle(taxa1+' '+taxa2+' PM score:'+str(round(pm_res.loc[pmi,'co_PM_2D'],3)))
                ax1 = plt.subplot(121)
                ax2 = plt.subplot(122)
                ax1.scatter(control_a, control_b, color='red',label=control_label)
                ax1.scatter(treat_a, treat_b, color='green',label=treat_label)
                ax1.set_xlabel(taxa1)
                ax1.set_ylabel(taxa2)
                ax1.legend()
                nstd = 2
                cov = np.cov(control_a, control_b)
                vals, vecs = eigsorted(cov)
                theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                w, h = 2 * nstd * np.sqrt(vals)
                ell = Ellipse(xy=(np.mean(control_a), np.mean(control_b)),
                              width=w, height=h,
                              angle=theta, color='red')
                ell.set_facecolor('red')
                ell.set_alpha(0.2)
                ax1.add_artist(ell)
                nstd = 2
                cov = np.cov(treat_a, treat_b)
                vals, vecs = eigsorted(cov)
                theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                w, h = 2 * nstd * np.sqrt(vals)
                ell_t = Ellipse(xy=(np.mean(treat_a), np.mean(treat_b)),
                              width=w, height=h,
                              angle=theta, color='green')
                ell_t.set_facecolor('green')
                ell_t.set_alpha(0.2)
                ax1.add_artist(ell_t)

                ax1.set_title('Relative abundance scatter plot')
                tmpilrdata = ilr_data[[taxa1, taxa2, condition]]
                tmpilrdata = tmpilrdata.dropna()
                tmpilrdata = tmpilrdata.reset_index(drop=True)
                control_a = tmpilrdata.loc[tmpilrdata[condition] == control_label, taxa1]
                control_b = tmpilrdata.loc[tmpilrdata[condition] == control_label, taxa2]
                treat_a = tmpilrdata.loc[tmpilrdata[condition] == treat_label, taxa1]
                treat_b = tmpilrdata.loc[tmpilrdata[condition] == treat_label, taxa2]
                ax2.scatter(control_a, control_b, color='red', label=control_label)
                ax2.scatter(treat_a, treat_b, color='green', label=treat_label)
                ax2.set_xlabel(taxa1)
                ax2.set_ylabel(taxa2)
                ax2.legend()
                nstd = 2
                cov = np.cov(control_a, control_b)
                vals, vecs = eigsorted(cov)
                theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                w, h = 2 * nstd * np.sqrt(vals)
                ell = Ellipse(xy=(np.mean(control_a), np.mean(control_b)),
                              width=w, height=h,
                              angle=theta, color='red')
                ell.set_facecolor('red')
                ell.set_alpha(0.2)
                ax2.add_artist(ell)
                nstd = 2
                cov = np.cov(treat_a, treat_b)
                vals, vecs = eigsorted(cov)
                theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                w, h = 2 * nstd * np.sqrt(vals)
                ell_t = Ellipse(xy=(np.mean(treat_a), np.mean(treat_b)),
                                width=w, height=h,
                                angle=theta, color='green')
                ell_t.set_facecolor('green')
                ell_t.set_alpha(0.2)
                ax2.add_artist(ell_t)
                ax2.set_title('ilr Relative abundance scatter plot')
                # plt.show()
                pdf.savefig()
                plt.close()

    plot_scatter()
    return None


def plot_t_squared_exe(project_dict, topnum, networktype, filter_ss):
    fileplace = project_dict.get('fileplace')
    fdr = project_dict.get('fdr')
    confidenceI = [0.99, 0.95, 0.90][project_dict.get('confidenceI')]
    resfilename = fileplace + '/PM_scores.csv'
    pm_res = pd.read_csv(resfilename, header=0, index_col=0)

    def co_PM_2D_fun(raw_pm_2d, pm1, pm2):
        return np.max([0, np.min([raw_pm_2d - pm1, raw_pm_2d - pm2])])

    if filter_ss == 'ON':
        if fdr == 'ON':
            pm_res.loc[pm_res['qvalue'] > 1 - confidenceI, 'raw_pm_2d'] = 0
            pm_res.loc[pm_res['qvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['qvalue2'] > 1 - confidenceI, 'pm2'] = 0
        if fdr == 'OFF':
            pm_res.loc[pm_res['pvalue'] > 1 - confidenceI, 'raw_pm_2d'] = 0
            pm_res.loc[pm_res['pvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['pvalue2'] > 1 - confidenceI, 'pm2'] = 0

    pm_res['co_PM_2D'] = pm_res.apply(lambda x: co_PM_2D_fun(x.raw_pm_2d, x.pm1, x.pm2), axis=1)
    target = 'raw_pm_2d'
    if networktype == 'Interaction PM scores':
        target = 'co_PM_2D'
    if networktype == '2D PM Scores':
        target = 'raw_pm_2d'
    pm_res = pm_res[pm_res[target] > 0]
    pm_res = pm_res.sort_values(by=target, ascending=False)
    pm_res = pm_res.reset_index(drop=True)
    taxnum = np.min([len(pm_res), topnum])
    pm_res = pm_res[0:taxnum]

    datafilename = project_dict.get('datafilename')
    control_label = project_dict.get('controllabel')
    treat_label = project_dict.get('treatlabel')
    condition = project_dict.get('grouplabel')
    multithreading = project_dict.get('multithreading')
    if multithreading == 'ON':
        threads_count = int(project_dict.get('threads'))
    if multithreading == 'OFF':
        threads_count = 1
    logs_place = project_dict.get('logs_place')
    image_pdf = fileplace + '/taxa_relationship_difference_between_groups.pdf'
    demodata = pd.read_csv(datafilename, header=0, index_col=False)
    taxaname = demodata.columns
    # convert absolute abundance to relative abundance
    genuslist = []
    for g in taxaname:
        if g != condition:
            genuslist.append(g)
    sumabundance = demodata[genuslist].sum(axis=1)
    demodata[genuslist] = demodata[genuslist].div(sumabundance, axis='rows')
    tmpd = demodata[genuslist]
    tmpd[tmpd <= 0] = np.nan
    demodata[genuslist] = tmpd
    ilr_data = copy.copy(demodata)
    ilr_data[genuslist] = np.log(ilr_data[genuslist])
    gmean = ilr_data[genuslist].mean(axis=1)
    ilr_data[genuslist] = ilr_data[genuslist].subtract(gmean, axis='rows')
    def plot_t_squared():
        with PdfPages(image_pdf) as pdf:
            for pmi in pm_res.index:
                taxa1=pm_res.loc[pmi,'taxa1']
                taxa2 = pm_res.loc[pmi, 'taxa2']
                cc, ct, tt, tc, cc_kde, ct_kde, tt_kde, tc_kde = None, None, None, None, None, None, None, None
                sucess = False
                for i in range(threads_count):
                    tmprunlog=None
                    tmprunlog = pd.read_csv(logs_place + '/' + str(i) + 'runlog.csv',index_col=0,header=0)
                    tmpr=tmprunlog[(tmprunlog['taxa1'] == taxa1)&(tmprunlog['taxa2']==taxa2)]
                    ct_list,ct_list,tt_list,tc_list,cc_kde_list,ct_kde_list,tt_kde_list,tc_kde_list= None, None, None, None, None, None, None, None
                    if len(tmpr)==1:
                        tmpr=tmpr.reset_index(drop=True)
                        numpyindex=tmpr.loc[0,'numpyindex']
                        cc_list=np.load(logs_place + '/' + str(i) + '_cc_list.npy',allow_pickle=True)
                        ct_list = np.load(logs_place + '/' + str(i) + '_ct_list.npy',allow_pickle=True)
                        tt_list = np.load(logs_place + '/' + str(i) + '_tt_list.npy',allow_pickle=True)
                        tc_list = np.load(logs_place + '/' + str(i) + '_tc_list.npy',allow_pickle=True)
                        cc_kde_list = np.load(logs_place + '/' + str(i) + '_cc_kde_list.npy',allow_pickle=True)
                        ct_kde_list = np.load(logs_place + '/' + str(i) + '_ct_kde_list.npy',allow_pickle=True)
                        tt_kde_list = np.load(logs_place + '/' + str(i) + '_tt_kde_list.npy',allow_pickle=True)
                        tc_kde_list = np.load(logs_place + '/' + str(i) + '_tc_kde_list.npy',allow_pickle=True)
                        cc=cc_list[numpyindex]
                        ct=ct_list[numpyindex]
                        tt=tt_list[numpyindex]
                        tc=tc_list[numpyindex]
                        cc_kde=cc_kde_list[numpyindex]
                        ct_kde=ct_kde_list[numpyindex]
                        tt_kde=tt_kde_list[numpyindex]
                        tc_kde=tc_kde_list[numpyindex]
                        # print(len(cc_list),len(ct_list),len(tt_list),len(tc_list),len(cc_kde_list),len(ct_kde_list),len(tt_kde_list),len(tc_kde_list))
                        sucess=True
                        break
                if (cc is not None) and (ct is not None) and (tt is not None) and (tc is not None) and sucess:
                    try:
                        fig= plt.figure(figsize=(20,21))
                        plt.suptitle(taxa1+' '+taxa2+' PM score:'+str(round(pm_res.loc[pmi,target],3)))
                        ax1 = plt.subplot(221)
                        ax2 = plt.subplot(222)

                        ax1.set_title('T-squared statistics base on the '+control_label+' group')
                        ax1.violinplot([cc,tc])
                        ax2.set_title('T-squared statistics base on the '+treat_label+' group')
                        ax2.violinplot([ct,tt])
                        # set style for the axes
                        labels = [control_label,treat_label]
                        for ax in [ax1, ax2]:
                            set_axis_style(ax, labels)
                        if (cc_kde is not None) and (tc_kde is not None):
                            ax5 = plt.subplot(223)
                            minv = np.min([np.min(cc), np.min(tc)])
                            maxv = np.max([np.max(cc), np.max(tc)])
                            X_plot = np.linspace(minv, maxv, 500)
                            cc_log_dens = cc_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                            tc_log_dens = tc_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                            ax5.plot(X_plot, np.exp(cc_log_dens), color='red', label=control_label)
                            ax5.plot(X_plot, np.exp(tc_log_dens), color='green', label=treat_label)
                            ax5.legend()
                            ax5.set_title(control_label + ' based T-squared statistics kernel estimated density')
                        else:
                            ax5 = plt.subplot(223)
                            ax5.hlines(1, 1, 2)
                            ax5.set_title(control_label + ' based T-squared statistics kernel estimated density not plotted')
                        if (ct_kde is not None) and (tt_kde is not None):
                            ax6 = plt.subplot(224)
                            minv = np.min([np.min(ct), np.min(tt)])
                            maxv = np.max([np.max(ct), np.max(tt)])
                            X_plot = np.linspace(minv, maxv, 500)
                            ct_log_dens = ct_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                            tt_log_dens = tt_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                            ax6.plot(X_plot, np.exp(ct_log_dens), color='red', label=control_label)
                            ax6.plot(X_plot, np.exp(tt_log_dens), color='green', label=treat_label)
                            ax6.legend()
                            ax6.set_title(treat_label + ' based T-squared statistics kernel estimated density')
                        else:
                            ax6 = plt.subplot(224)
                            ax6.hlines(1, 1, 2)
                            ax6.set_title(treat_label + ' based T-squared statistics kernel estimated density not plotted')
                        # plt.show()
                        pdf.savefig()
                        plt.close()
                    except Exception as e:
                        print(e)
                        print(taxa1,taxa2)
    plot_t_squared()
    return None


def plot_scatter_and_tsquared_exe(project_dict, topnum, networktype, filter_ss):
    fileplace = project_dict.get('fileplace')
    fdr = project_dict.get('fdr')
    confidenceI = [0.99, 0.95, 0.90][project_dict.get('confidenceI')]
    resfilename = fileplace + '/PM_scores.csv'
    pm_res = pd.read_csv(resfilename, header=0, index_col=0)

    def co_PM_2D_fun(raw_pm_2d, pm1, pm2):
        return np.max([0, np.min([raw_pm_2d - pm1, raw_pm_2d - pm2])])

    if filter_ss == 'ON':
        if fdr == 'ON':
            pm_res.loc[pm_res['qvalue'] > 1 - confidenceI, 'raw_pm_2d'] = 0
            pm_res.loc[pm_res['qvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['qvalue2'] > 1 - confidenceI, 'pm2'] = 0
        if fdr == 'OFF':
            pm_res.loc[pm_res['pvalue'] > 1 - confidenceI, 'raw_pm_2d'] = 0
            pm_res.loc[pm_res['pvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['pvalue2'] > 1 - confidenceI, 'pm2'] = 0

    pm_res['co_PM_2D'] = pm_res.apply(lambda x: co_PM_2D_fun(x.raw_pm_2d, x.pm1, x.pm2), axis=1)
    target = 'raw_pm_2d'
    if networktype == 'Interaction PM scores':
        target = 'co_PM_2D'
    if networktype == '2D PM Scores':
        target = 'raw_pm_2d'
    pm_res = pm_res[pm_res[target] > 0]
    pm_res = pm_res.sort_values(by=target, ascending=False)
    pm_res = pm_res.reset_index(drop=True)
    taxnum = np.min([len(pm_res), topnum])
    pm_res = pm_res[0:taxnum]

    datafilename = project_dict.get('datafilename')
    control_label = project_dict.get('controllabel')
    treat_label = project_dict.get('treatlabel')
    condition = project_dict.get('grouplabel')
    multithreading = project_dict.get('multithreading')
    if multithreading == 'ON':
        threads_count = int(project_dict.get('threads'))
    if multithreading == 'OFF':
        threads_count = 1
    logs_place = project_dict.get('logs_place')
    image_pdf = fileplace + '/taxa_relationship_difference_between_groups.pdf'
    demodata = pd.read_csv(datafilename, header=0, index_col=False)
    taxaname = demodata.columns
    # convert absolute abundance to relative abundance
    genuslist = []
    for g in taxaname:
        if g != condition:
            genuslist.append(g)
    sumabundance = demodata[genuslist].sum(axis=1)
    demodata[genuslist] = demodata[genuslist].div(sumabundance, axis='rows')
    tmpd = demodata[genuslist]
    tmpd[tmpd <= 0] = np.nan
    demodata[genuslist] = tmpd
    ilr_data = copy.copy(demodata)
    ilr_data[genuslist] = np.log(ilr_data[genuslist])
    gmean = ilr_data[genuslist].mean(axis=1)
    ilr_data[genuslist] = ilr_data[genuslist].subtract(gmean, axis='rows')
    def plot_scatter_and_tsquared():
        with PdfPages(image_pdf) as pdf:
            for pmi in pm_res.index:
                taxa1 = pm_res.loc[pmi, 'taxa1']
                taxa2 = pm_res.loc[pmi, 'taxa2']
                cc, ct, tt, tc, cc_kde, ct_kde, tt_kde, tc_kde = None, None, None, None, None, None, None, None
                sucess = False
                for i in range(threads_count):
                    tmprunlog = None
                    tmprunlog = pd.read_csv(logs_place + '/' + str(i) + 'runlog.csv', index_col=0, header=0)
                    tmpr = tmprunlog[(tmprunlog['taxa1'] == taxa1) & (tmprunlog['taxa2'] == taxa2)]
                    ct_list,ct_list,tt_list,tc_list,cc_kde_list,ct_kde_list,tt_kde_list,tc_kde_list= None, None, None, None, None, None, None, None

                    if len(tmpr) == 1:
                        tmpr = tmpr.reset_index(drop=True)
                        numpyindex = tmpr.loc[0, 'numpyindex']
                        cc_list = np.load(logs_place + '/' + str(i) + '_cc_list.npy',allow_pickle=True)
                        ct_list = np.load(logs_place + '/' + str(i) + '_ct_list.npy',allow_pickle=True)
                        tt_list = np.load(logs_place + '/' + str(i) + '_tt_list.npy',allow_pickle=True)
                        tc_list = np.load(logs_place + '/' + str(i) + '_tc_list.npy',allow_pickle=True)
                        cc_kde_list = np.load(logs_place + '/' + str(i) + '_cc_kde_list.npy',allow_pickle=True)
                        ct_kde_list = np.load(logs_place + '/' + str(i) + '_ct_kde_list.npy',allow_pickle=True)
                        tt_kde_list = np.load(logs_place + '/' + str(i) + '_tt_kde_list.npy',allow_pickle=True)
                        tc_kde_list = np.load(logs_place + '/' + str(i) + '_tc_kde_list.npy',allow_pickle=True)
                        cc = cc_list[numpyindex]
                        ct = ct_list[numpyindex]
                        tt = tt_list[numpyindex]
                        tc = tc_list[numpyindex]
                        cc_kde = cc_kde_list[numpyindex]
                        ct_kde = ct_kde_list[numpyindex]
                        tt_kde = tt_kde_list[numpyindex]
                        tc_kde = tc_kde_list[numpyindex]
                        sucess = True
                        break
                if (cc is not None) and (ct is not None) and (tt is not None) and (tc is not None) and sucess:
                    fig = plt.figure(figsize=(20, 30))
                    plt.suptitle(taxa1 + ' ' + taxa2 + ' PM score:' + str(round(pm_res.loc[pmi, target], 3)))
                    ax1 = plt.subplot(321)
                    ax2 = plt.subplot(322)
                    ax3 = plt.subplot(323)
                    ax4 = plt.subplot(324)
                    tmpdemodata = demodata[[taxa1, taxa2, condition]]
                    tmpdemodata = tmpdemodata.dropna()
                    tmpdemodata = tmpdemodata.reset_index(drop=True)
                    control_a = tmpdemodata.loc[tmpdemodata[condition] == control_label, taxa1]
                    control_b = tmpdemodata.loc[tmpdemodata[condition] == control_label, taxa2]
                    treat_a = tmpdemodata.loc[tmpdemodata[condition] == treat_label, taxa1]
                    treat_b = tmpdemodata.loc[tmpdemodata[condition] == treat_label, taxa2]


                    ax1.scatter(control_a, control_b, color='red', label=control_label)
                    ax1.scatter(treat_a, treat_b, color='green', label=treat_label)
                    ax1.set_xlabel(taxa1)
                    ax1.set_ylabel(taxa2)
                    ax1.legend()
                    nstd = 2
                    cov = np.cov(control_a, control_b)
                    vals, vecs = eigsorted(cov)
                    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                    w, h = 2 * nstd * np.sqrt(vals)
                    ell = Ellipse(xy=(np.mean(control_a), np.mean(control_b)),
                                  width=w, height=h,
                                  angle=theta, color='red')
                    ell.set_facecolor('red')
                    ell.set_alpha(0.2)
                    ax1.add_artist(ell)

                    nstd = 2
                    cov = np.cov(treat_a, treat_b)
                    vals, vecs = eigsorted(cov)
                    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                    w, h = 2 * nstd * np.sqrt(vals)
                    ell_t = Ellipse(xy=(np.mean(treat_a), np.mean(treat_b)),
                                    width=w, height=h,
                                    angle=theta, color='green')
                    ell_t.set_facecolor('green')
                    ell_t.set_alpha(0.2)
                    ax1.add_artist(ell_t)

                    ax1.set_title('Relative abundance scatter plot')

                    tmpilrdata = ilr_data[[taxa1, taxa2, condition]]
                    tmpilrdata = tmpilrdata.dropna()
                    tmpilrdata = tmpilrdata.reset_index(drop=True)
                    control_a = tmpilrdata.loc[tmpilrdata[condition] == control_label, taxa1]
                    control_b = tmpilrdata.loc[tmpilrdata[condition] == control_label, taxa2]
                    treat_a = tmpilrdata.loc[tmpilrdata[condition] == treat_label, taxa1]
                    treat_b = tmpilrdata.loc[tmpilrdata[condition] == treat_label, taxa2]
                    ax2.scatter(control_a, control_b, color='red', label=control_label)
                    ax2.scatter(treat_a, treat_b, color='green', label=treat_label)
                    ax2.set_xlabel(taxa1)
                    ax2.set_ylabel(taxa2)
                    nstd = 2
                    cov = np.cov(control_a, control_b)
                    vals, vecs = eigsorted(cov)
                    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                    w, h = 2 * nstd * np.sqrt(vals)
                    ell = Ellipse(xy=(np.mean(control_a), np.mean(control_b)),
                                  width=w, height=h,
                                  angle=theta, color='red')
                    ell.set_facecolor('red')
                    ell.set_alpha(0.2)
                    ax2.add_artist(ell)

                    nstd = 2
                    cov = np.cov(treat_a, treat_b)
                    vals, vecs = eigsorted(cov)
                    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                    w, h = 2 * nstd * np.sqrt(vals)
                    ell_t = Ellipse(xy=(np.mean(treat_a), np.mean(treat_b)),
                                    width=w, height=h,
                                    angle=theta, color='green')
                    ell_t.set_facecolor('green')
                    ell_t.set_alpha(0.2)
                    ax2.add_artist(ell_t)
                    ax2.legend()
                    ax2.set_title('ilr Relative abundance scatter plot')

                    ax3.set_title('T-squared statistics base on the ' + control_label + ' group')
                    ax3.violinplot([cc, tc])
                    ax4.set_title('T-squared statistics base on the ' + treat_label + ' group')
                    ax4.violinplot([ct, tt])
                    # set style for the axes
                    labels = [control_label, treat_label]
                    for ax in [ax3, ax4]:
                        set_axis_style(ax, labels)
                    if (cc_kde is not None) and (tc_kde is not None):
                        ax5 = plt.subplot(325)
                        minv = np.min([np.min(cc), np.min(tc)])
                        maxv = np.max([np.max(cc), np.max(tc)])
                        X_plot = np.linspace(minv, maxv, 500)
                        cc_log_dens = cc_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                        tc_log_dens = tc_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                        ax5.plot(X_plot, np.exp(cc_log_dens), color='red', label=control_label)
                        ax5.plot(X_plot, np.exp(tc_log_dens), color='green', label=treat_label)
                        ax5.legend()
                        ax5.set_title(control_label + ' based T-squared statistics kernel estimated density')
                    else:
                        ax5 = plt.subplot(325)
                        ax5.hlines(1, 1, 2)
                        ax5.set_title(control_label + ' based T-squared statistics kernel estimated density not plotted')
                    if (ct_kde is not None) and (tt_kde is not None):
                        ax6 = plt.subplot(326)
                        minv = np.min([np.min(ct), np.min(tt)])
                        maxv = np.max([np.max(ct), np.max(tt)])
                        X_plot = np.linspace(minv, maxv, 500)
                        ct_log_dens = ct_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                        tt_log_dens = tt_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                        ax6.plot(X_plot, np.exp(ct_log_dens), color='red', label=control_label)
                        ax6.plot(X_plot, np.exp(tt_log_dens), color='green', label=treat_label)
                        ax6.legend()
                        ax6.set_title(treat_label + ' based T-squared statistics kernel estimated density')
                    else:
                        ax6 = plt.subplot(326)
                        ax6.hlines(1, 1, 2)
                        ax6.set_title(treat_label + ' based T-squared statistics kernel estimated density not plotted')
                    # plt.show()
                    pdf.savefig()
                    plt.close()
    plot_scatter_and_tsquared()
    return None