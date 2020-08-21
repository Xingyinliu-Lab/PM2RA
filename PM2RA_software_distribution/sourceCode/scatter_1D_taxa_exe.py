import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import copy



def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Groups')

def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]


def plot_scatter_exe(project_dict,filter_ss,taxalist_plot,scattertype):
    scatter = scattertype
    fileplace = project_dict.get('fileplace')
    fdr = project_dict.get('fdr')
    confidenceI= [0.99, 0.95,0.90][project_dict.get('confidenceI')]
    resfilename_1D = fileplace+'/PM_scores_1D.csv'
    pm_res_1D = pd.read_csv(resfilename_1D, header=0, index_col=0)
    if filter_ss=='ON':
        if fdr=='ON':
            pm_res_1D=pm_res_1D[pm_res_1D['qvalue']<1-confidenceI]
        if fdr == 'OFF':
            pm_res_1D = pm_res_1D[pm_res_1D['pvalue'] < 1 - confidenceI]
    pm_res_1D = pm_res_1D.sort_values(by='pm', ascending=False)
    pm_res_1D=pm_res_1D[pm_res_1D['taxa']==taxalist_plot[0]]
    pm_res_1D=pm_res_1D.reset_index(drop=True)
    taxalist = list(pm_res_1D['taxa'])
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
    image_pdf = fileplace + '/taxa_difference_between_groups_withintaxa.pdf'
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
            for taxa in taxalist:
                fig= plt.figure(figsize=(20,10.5))
                tmpr = pm_res_1D[(pm_res_1D['taxa'] == taxa)]
                tmpr = tmpr.reset_index(drop=True)
                plt.suptitle(taxa + ' PM score:' + str(round(tmpr.loc[0, 'pm'], 3)))
                ax1 = plt.subplot(121)
                ax2 = plt.subplot(122)
                tmpdemodata = demodata[[taxa,condition]]
                tmpdemodata = tmpdemodata.dropna()
                tmpdemodata = tmpdemodata.reset_index(drop=True)
                control_a = tmpdemodata.loc[tmpdemodata[condition] == control_label, taxa]
                treat_a = tmpdemodata.loc[tmpdemodata[condition] == treat_label, taxa]
                control_a=list(control_a)
                treat_a = list(treat_a)
                if scatter=='Violin Chart':
                    ax1.violinplot([control_a, treat_a])
                    labels = [control_label, treat_label]
                    set_axis_style(ax1, labels)
                    ax1.set_title('Relative abundance violin plot')
                if scatter == 'Box Chart':
                    ax1.boxplot([control_a, treat_a])
                    labels = [control_label, treat_label]
                    set_axis_style(ax1, labels)
                    ax1.set_title('Relative abundance box plot')
                if scatter == 'Density Chart':
                    ax1.hist(control_a,color='red',label=control_label,density=True,alpha=0.3)
                    ax1.hist(treat_a,color='green',label=treat_label,density=True,alpha=0.3)
                    ax1.legend()
                    ax1.set_title('Relative abundance density plot')
                tmpilrdata = ilr_data[[taxa, condition]]
                tmpilrdata = tmpilrdata.dropna()
                tmpilrdata = tmpilrdata.reset_index(drop=True)
                control_a = tmpilrdata.loc[tmpilrdata[condition] == control_label, taxa]
                treat_a = tmpilrdata.loc[tmpilrdata[condition] == treat_label, taxa]
                control_a=list(control_a)
                treat_a = list(treat_a)
                if scatter == 'Violin Chart':
                    ax2.violinplot([control_a, treat_a])
                    labels = [control_label, treat_label]
                    set_axis_style(ax2, labels)
                    ax2.set_title('ilr Relative abundance violin plot')
                if scatter == 'Box Chart':
                    ax2.boxplot([control_a, treat_a])
                    labels = [control_label, treat_label]
                    set_axis_style(ax2, labels)
                    ax2.set_title('ilr Relative abundance box plot')
                if scatter == 'Density Chart':
                    ax2.hist(control_a,color='red',label=control_label,density=True,alpha=0.3)
                    ax2.hist(treat_a,color='green',label=treat_label,density=True,alpha=0.3)
                    ax2.legend()
                    ax2.set_title('ilr Relative abundance density plot')
                # plt.show()
                pdf.savefig()
                plt.close()
    plot_scatter()
    return None


def plot_t_squared_exe(project_dict,filter_ss,taxalist_plot,scattertype):
    scatter = scattertype
    fileplace = project_dict.get('fileplace')
    fdr = project_dict.get('fdr')
    confidenceI = [0.99, 0.95, 0.90][project_dict.get('confidenceI')]
    resfilename_1D = fileplace + '/PM_scores_1D.csv'
    pm_res_1D = pd.read_csv(resfilename_1D, header=0, index_col=0)
    if filter_ss == 'ON':
        if fdr == 'ON':
            pm_res_1D = pm_res_1D[pm_res_1D['qvalue'] < 1 - confidenceI]
        if fdr == 'OFF':
            pm_res_1D = pm_res_1D[pm_res_1D['pvalue'] < 1 - confidenceI]
    pm_res_1D = pm_res_1D.sort_values(by='pm', ascending=False)
    pm_res_1D = pm_res_1D[pm_res_1D['taxa'] == taxalist_plot[0]]
    pm_res_1D = pm_res_1D.reset_index(drop=True)
    taxalist = list(pm_res_1D['taxa'])
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
    image_pdf = fileplace + '/taxa_difference_between_groups_withintaxa.pdf'
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
            for taxa in taxalist:
                cc, ct, tt, tc, cc_kde, ct_kde, tt_kde, tc_kde = None, None, None, None, None, None, None, None
                sucess = False
                for i in range(threads_count):
                    tmprunlog = None
                    tmprunlog = pd.read_csv(logs_place + '/' + str(i) + 'runlog_1D.csv', index_col=0, header=0)
                    tmpr = tmprunlog[(tmprunlog['taxa'] == taxa)]
                    ct_list,ct_list,tt_list,tc_list,cc_kde_list,ct_kde_list,tt_kde_list,tc_kde_list= None, None, None, None, None, None, None, None
                    if len(tmpr) == 1:
                        tmpr = tmpr.reset_index(drop=True)
                        numpyindex = tmpr.loc[0, 'numpyindex']
                        cc_list = np.load(logs_place + '/' + str(i) + '_cc_list_1D.npy',allow_pickle=True)
                        ct_list = np.load(logs_place + '/' + str(i) + '_ct_list_1D.npy',allow_pickle=True)
                        tt_list = np.load(logs_place + '/' + str(i) + '_tt_list_1D.npy',allow_pickle=True)
                        tc_list = np.load(logs_place + '/' + str(i) + '_tc_list_1D.npy',allow_pickle=True)
                        cc_kde_list = np.load(logs_place + '/' + str(i) + '_cc_kde_list_1D.npy',allow_pickle=True)
                        ct_kde_list = np.load(logs_place + '/' + str(i) + '_ct_kde_list_1D.npy',allow_pickle=True)
                        tt_kde_list = np.load(logs_place + '/' + str(i) + '_tt_kde_list_1D.npy',allow_pickle=True)
                        tc_kde_list = np.load(logs_place + '/' + str(i) + '_tc_kde_list_1D.npy',allow_pickle=True)
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
                    fig = plt.figure(figsize=(20, 21))
                    tmpr = pm_res_1D[(pm_res_1D['taxa'] == taxa)]
                    tmpr = tmpr.reset_index(drop=True)
                    plt.suptitle(taxa + ' PM score:' + str(round(tmpr.loc[0, 'pm'], 3)))
                    ax1 = plt.subplot(221)
                    ax2 = plt.subplot(222)
                    ax1.set_title('T-squared statistics base on the ' + control_label + ' group')
                    ax1.violinplot([cc, tc])
                    ax2.set_title('T-squared statistics base on the ' + treat_label + ' group')
                    ax2.violinplot([ct, tt])
                    # set style for the axes
                    labels = [control_label, treat_label]
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
    plot_t_squared()
    return None

def plot_scatter_and_tsquared_exe(project_dict,filter_ss,taxalist_plot,scattertype):
    scatter = scattertype
    fileplace = project_dict.get('fileplace')
    fdr = project_dict.get('fdr')
    confidenceI = [0.99, 0.95, 0.90][project_dict.get('confidenceI')]
    resfilename_1D = fileplace + '/PM_scores_1D.csv'
    pm_res_1D = pd.read_csv(resfilename_1D, header=0, index_col=0)
    if filter_ss == 'ON':
        if fdr == 'ON':
            pm_res_1D = pm_res_1D[pm_res_1D['qvalue'] < 1 - confidenceI]
        if fdr == 'OFF':
            pm_res_1D = pm_res_1D[pm_res_1D['pvalue'] < 1 - confidenceI]
    pm_res_1D = pm_res_1D.sort_values(by='pm', ascending=False)
    pm_res_1D = pm_res_1D[pm_res_1D['taxa'] == taxalist_plot[0]]
    pm_res_1D = pm_res_1D.reset_index(drop=True)
    taxalist = list(pm_res_1D['taxa'])
    datafilename = project_dict.get('datafilename')
    control_label = project_dict.get('controllabel')
    treat_label = project_dict.get('treatlabel')
    condition = project_dict.get('grouplabel')
    multithreading = project_dict.get('multithreading')
    if multithreading == 'ON':
        threads_count = int(project_dict.get('threads'))
    if multithreading == 'OFF':
        threads_count = 1
    logs_place =  project_dict.get('logs_place')
    image_pdf = fileplace + '/taxa_difference_between_groups_withintaxa.pdf'
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
            for taxa in taxalist:
                cc, ct, tt, tc, cc_kde, ct_kde, tt_kde, tc_kde = None, None, None, None, None, None, None, None
                sucess = False
                for i in range(threads_count):
                    tmprunlog = None
                    tmprunlog = pd.read_csv(logs_place + '/' + str(i) + 'runlog_1D.csv', index_col=0, header=0)
                    tmpr = tmprunlog[(tmprunlog['taxa'] == taxa)]
                    ct_list,ct_list,tt_list,tc_list,cc_kde_list,ct_kde_list,tt_kde_list,tc_kde_list= None, None, None, None, None, None, None, None
                    if len(tmpr) == 1:
                        tmpr = tmpr.reset_index(drop=True)
                        numpyindex = tmpr.loc[0, 'numpyindex']
                        cc_list = np.load(logs_place + '/' + str(i) + '_cc_list_1D.npy',allow_pickle=True)
                        ct_list = np.load(logs_place + '/' + str(i) + '_ct_list_1D.npy',allow_pickle=True)
                        tt_list = np.load(logs_place + '/' + str(i) + '_tt_list_1D.npy',allow_pickle=True)
                        tc_list = np.load(logs_place + '/' + str(i) + '_tc_list_1D.npy',allow_pickle=True)
                        cc_kde_list = np.load(logs_place + '/' + str(i) + '_cc_kde_list_1D.npy',allow_pickle=True)
                        ct_kde_list = np.load(logs_place + '/' + str(i) + '_ct_kde_list_1D.npy',allow_pickle=True)
                        tt_kde_list = np.load(logs_place + '/' + str(i) + '_tt_kde_list_1D.npy',allow_pickle=True)
                        tc_kde_list = np.load(logs_place + '/' + str(i) + '_tc_kde_list_1D.npy',allow_pickle=True)
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
                    fig = plt.figure(figsize=(20, 31))
                    tmpr = pm_res_1D[(pm_res_1D['taxa'] == taxa)]
                    tmpr = tmpr.reset_index(drop=True)
                    plt.suptitle(taxa + ' PM score:' + str(round(tmpr.loc[0, 'pm'], 3)))
                    ax1 = plt.subplot(321)
                    ax2 = plt.subplot(322)
                    ax3 = plt.subplot(323)
                    ax4 = plt.subplot(324)
                    tmpdemodata = demodata[[taxa, condition]]
                    tmpdemodata = tmpdemodata.dropna()
                    tmpdemodata = tmpdemodata.reset_index(drop=True)
                    control_a = tmpdemodata.loc[tmpdemodata[condition] == control_label, taxa]
                    treat_a = tmpdemodata.loc[tmpdemodata[condition] == treat_label, taxa]
                    control_a = list(control_a)
                    treat_a = list(treat_a)
                    if scatter == 'Violin Chart':
                        ax1.violinplot([control_a, treat_a])
                        labels = [control_label, treat_label]
                        set_axis_style(ax1, labels)
                        ax1.set_title('Relative abundance violin plot')
                    if scatter == 'Box Chart':
                        ax1.boxplot([control_a, treat_a])
                        labels = [control_label, treat_label]
                        set_axis_style(ax1, labels)
                        ax1.set_title('Relative abundance box plot')
                    if scatter == 'Density Chart':
                        ax1.hist(control_a, color='red', label=control_label, density=True, alpha=0.3)
                        ax1.hist(treat_a, color='green', label=treat_label, density=True, alpha=0.3)
                        ax1.legend()
                        ax1.set_title('Relative abundance density plot')

                    tmpilrdata = ilr_data[[taxa, condition]]
                    tmpilrdata = tmpilrdata.dropna()
                    tmpilrdata = tmpilrdata.reset_index(drop=True)
                    control_a = tmpilrdata.loc[tmpilrdata[condition] == control_label, taxa]
                    treat_a = tmpilrdata.loc[tmpilrdata[condition] == treat_label, taxa]
                    control_a = list(control_a)
                    treat_a = list(treat_a)
                    if scatter == 'Violin Chart':
                        ax2.violinplot([control_a, treat_a])
                        labels = [control_label, treat_label]
                        set_axis_style(ax2, labels)
                        ax2.set_title('ilr Relative abundance violin plot')
                    if scatter == 'Box Chart':
                        ax2.boxplot([control_a, treat_a])
                        labels = [control_label, treat_label]
                        set_axis_style(ax2, labels)
                        ax2.set_title('ilr Relative abundance box plot')
                    if scatter == 'Density Chart':
                        ax2.hist(control_a, color='red', label=control_label, density=True, alpha=0.3)
                        ax2.hist(treat_a, color='green', label=treat_label, density=True, alpha=0.3)
                        ax2.legend()
                        ax2.set_title('ilr Relative abundance density plot')

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














