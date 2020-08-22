import pandas as pd
from scipy.integrate import quad
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import numpy as np
from scipy.stats import ks_2samp
from matplotlib.patches import Ellipse
import fdrcorrection
import load_save_project
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import copy
np.seterr(divide='ignore',invalid='ignore')
from scipy.spatial import distance_matrix
from PCoA import pcoa
import warnings
warnings.filterwarnings("ignore")
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


def cal_pm_score_nd(datafilename,control,treat,condition,ilr,wetherFindBestBandwidth,kernel_str,module_taxa):
    def hotellingt2d(data,data_average,data_cov):
        data_cov = np.matrix(data_cov)
        try:
            data_cov_r=np.linalg.inv(data_cov)
        except:
            data_cov_r = np.linalg.pinv(data_cov)

        cdata=data-data_average.reshape(-1,1)
        t2d=np.diag(np.matmul(np.matmul(np.transpose(cdata), data_cov_r),cdata))
        return t2d

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
        if len(tmpdemodata)<5:
            return None, None, None, None, None, None, None, None, None, None

        control_obs=tmpdemodata[tmpdemodata[condition]==control]
        control_obs=np.transpose(np.asarray(control_obs[taxalist]))
        treat_obs=tmpdemodata[tmpdemodata[condition]==treat]
        treat_obs = np.transpose(np.asarray(treat_obs[taxalist]))

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

    if ilr=='ON':
        demodata[genuslist] = np.log(demodata[genuslist])
        gmean = demodata[genuslist].mean(axis=1)
        demodata[genuslist] = demodata[genuslist].subtract(gmean, axis='rows')
    # filter low detection taxa
    genuslist = []
    for g in taxaname:
        if g != condition:
            genuslist.append(g)
    demodata=demodata[genuslist+list([condition])]
    pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde=cal_pm_socre_nD(module_taxa,demodata)
    return pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde

def plot_pm_score_nd(fileplace,datafilename,control,treat,condition,module_taxa,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde):
    if len(module_taxa)==1:
        plot_scatter_and_tsquared_1d_exe(fileplace,datafilename,control,treat,condition,module_taxa,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde)
    if len(module_taxa)==2:
        plot_scatter_and_tsquared_2d_exe(fileplace,datafilename,control,treat,condition,module_taxa,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde)
    if len(module_taxa)>2:
        plot_scatter_and_tsquared_nd_exe(fileplace,datafilename,control,treat,condition,module_taxa,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde)
    return None




def plot_scatter_and_tsquared_nd_exe(fileplace,datafilename,control,treat,condition,module_taxa,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde):
    image_pdf = fileplace + '/module_difference_between_groups.pdf'
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
    control_label=control
    treat_label=treat
    def plot_scatter_and_tsquared():
        with PdfPages(image_pdf) as pdf:
            if (cc is not None) and (ct is not None) and (tt is not None) and (tc is not None) :
                fig = plt.figure(figsize=(20, 30))
                plt.suptitle( ' PM score:' + str(round(pm_score, 3)) + ' Pvalue: ' + str(round(pvalue, 6)))
                ax1 = plt.subplot(321)
                ax2 = plt.subplot(322)
                ax3 = plt.subplot(323)
                ax4 = plt.subplot(324)

                tmpdemodata = demodata[module_taxa+[condition]]
                tmpdemodata = tmpdemodata.dropna()
                tmpdemodata = tmpdemodata.reset_index(drop=True)
                control_a = tmpdemodata.loc[tmpdemodata[condition] == control_label, module_taxa]
                treat_a = tmpdemodata.loc[tmpdemodata[condition] == treat_label, module_taxa]

                control_a = pd.DataFrame(control_a)
                index_a=control_a.columns
                control_a = pd.DataFrame(distance_matrix(control_a.values, control_a.values))
                control_a = control_a.values
                control_a = pcoa(control_a,  method="eigh", number_of_dimensions=2, inplace=False)

                treat_a = pd.DataFrame(treat_a)
                index_a=treat_a.columns
                treat_a = pd.DataFrame(distance_matrix(treat_a.values, treat_a.values))
                treat_a = treat_a.values
                treat_a = pcoa(treat_a,  method="eigh", number_of_dimensions=2, inplace=False)

                control_b = np.transpose(np.asarray(control_a['PC2']))
                control_a =  np.transpose(np.asarray(control_a['PC1']))

                treat_b = np.transpose(np.asarray(treat_a['PC2']))
                treat_a = np.transpose(np.asarray(treat_a['PC1']))

                ax1.scatter(control_a, control_b, color='red', label=control_label)
                ax1.scatter(treat_a, treat_b, color='green', label=treat_label)
                ax1.set_xlabel('PCoA1')
                ax1.set_ylabel('PCoA2')
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

                tmpilrdata = ilr_data[module_taxa+[condition]]
                tmpilrdata = tmpilrdata.dropna()
                tmpilrdata = tmpilrdata.reset_index(drop=True)
                control_a = tmpilrdata.loc[tmpilrdata[condition] == control_label, module_taxa]
                treat_a = tmpilrdata.loc[tmpilrdata[condition] == treat_label, module_taxa]

                control_a = pd.DataFrame(control_a)
                index_a = control_a.columns
                control_a = pd.DataFrame(distance_matrix(control_a.values, control_a.values))
                control_a = control_a.values
                control_a = pcoa(control_a,  method="eigh", number_of_dimensions=2, inplace=False)

                treat_a = pd.DataFrame(treat_a)
                index_a = treat_a.columns
                treat_a = pd.DataFrame(distance_matrix(treat_a.values, treat_a.values))
                treat_a = treat_a.values
                treat_a = pcoa(treat_a,  method="eigh", number_of_dimensions=2, inplace=False)

                control_b = np.transpose(np.asarray(control_a['PC2']))
                control_a =  np.transpose(np.asarray(control_a['PC1']))

                treat_b = np.transpose(np.asarray(treat_a['PC2']))
                treat_a = np.transpose(np.asarray(treat_a['PC1']))
                ax2.scatter(control_a, control_b, color='red', label=control_label)
                ax2.scatter(treat_a, treat_b, color='green', label=treat_label)
                ax2.set_xlabel('PCoA1')
                ax2.set_ylabel('PCoA2')
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



def plot_scatter_and_tsquared_2d_exe(fileplace,datafilename,control,treat,condition,module_taxa,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde):
    image_pdf = fileplace + '/module_difference_between_groups.pdf'
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
    taxa1=module_taxa[0]
    taxa2=module_taxa[1]
    control_label=control
    treat_label=treat
    def plot_scatter_and_tsquared():
        with PdfPages(image_pdf) as pdf:
            if (cc is not None) and (ct is not None) and (tt is not None) and (tc is not None):
                fig = plt.figure(figsize=(20, 30))
                plt.suptitle(taxa1 + ' ' + taxa2 +  ' PM score:' + str(round(pm_score, 3)) + ' Pvalue: ' + str(round(pvalue, 6)))
                ax1 = plt.subplot(321)
                ax2 = plt.subplot(322)
                ax3 = plt.subplot(323)
                ax4 = plt.subplot(324)
                ax5 = plt.subplot(325)
                ax6 = plt.subplot(326)
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

def plot_scatter_and_tsquared_1d_exe(fileplace,datafilename,control,treat,condition,module_taxa,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde):
    image_pdf = fileplace + '/module_difference_between_groups.pdf'
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
            if (cc is not None) and (ct is not None) and (tt is not None) and (tc is not None)and (cc_kde is not None) and (ct_kde is not None)and (tt_kde is not None)and (tc_kde is not None):
                fig = plt.figure(figsize=(20, 31))
                plt.suptitle(module_taxa[0] + ' PM score:' + str(round(pm_score, 3)) + ' Pvalue: ' + str(round(pvalue, 6)))
                ax1 = plt.subplot(321)
                ax2 = plt.subplot(322)
                ax3 = plt.subplot(323)
                ax4 = plt.subplot(324)
                ax5 = plt.subplot(325)
                ax6 = plt.subplot(326)
                tmpdemodata = demodata[[module_taxa[0], condition]]
                tmpdemodata = tmpdemodata.dropna()
                tmpdemodata = tmpdemodata.reset_index(drop=True)
                control_a = tmpdemodata.loc[tmpdemodata[condition] == control, module_taxa[0]]
                treat_a = tmpdemodata.loc[tmpdemodata[condition] == treat, module_taxa[0]]
                control_a = list(control_a)
                treat_a = list(treat_a)

                ax1.hist(control_a, color='red', label=control, density=True, alpha=0.3)
                ax1.hist(treat_a, color='green', label=treat, density=True, alpha=0.3)
                ax1.legend()
                ax1.set_title('Relative abundance density plot')

                tmpilrdata = ilr_data[[module_taxa[0], condition]]
                tmpilrdata = tmpilrdata.dropna()
                tmpilrdata = tmpilrdata.reset_index(drop=True)
                control_a = tmpilrdata.loc[tmpilrdata[condition] ==control, module_taxa[0]]
                treat_a = tmpilrdata.loc[tmpilrdata[condition] == treat, module_taxa[0]]
                control_a = list(control_a)
                treat_a = list(treat_a)

                ax2.hist(control_a, color='red', label=control, density=True, alpha=0.3)
                ax2.hist(treat_a, color='green', label=treat, density=True, alpha=0.3)
                ax2.legend()
                ax2.set_title('ilr Relative abundance density plot')

                ax3.set_title('T-squared statistics base on the ' + control + ' group')
                ax3.violinplot([cc, tc])
                ax4.set_title('T-squared statistics base on the ' + treat + ' group')
                ax4.violinplot([ct, tt])
                # set style for the axes
                labels = [control, treat]
                for ax in [ax3, ax4]:
                    set_axis_style(ax, labels)
                if (cc_kde is not None) and (tc_kde is not None):
                    ax5 = plt.subplot(325)
                    minv = np.min([np.min(cc), np.min(tc)])
                    maxv = np.max([np.max(cc), np.max(tc)])
                    X_plot = np.linspace(minv, maxv, 500)
                    cc_log_dens = cc_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                    tc_log_dens = tc_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                    ax5.plot(X_plot, np.exp(cc_log_dens), color='red', label=control)
                    ax5.plot(X_plot, np.exp(tc_log_dens), color='green', label=treat)
                    ax5.legend()
                    ax5.set_title(control + ' based T-squared statistics kernel estimated density')
                else:
                    ax5 = plt.subplot(325)
                    ax5.hlines(1, 1, 2)
                    ax5.set_title(control + ' based T-squared statistics kernel estimated density not plotted')
                if (ct_kde is not None) and (tt_kde is not None):
                    ax6 = plt.subplot(326)
                    minv = np.min([np.min(ct), np.min(tt)])
                    maxv = np.max([np.max(ct), np.max(tt)])
                    X_plot = np.linspace(minv, maxv, 500)
                    ct_log_dens = ct_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                    tt_log_dens = tt_kde.score_samples(np.array([X_plot]).reshape(-1, 1))
                    ax6.plot(X_plot, np.exp(ct_log_dens), color='red', label=control)
                    ax6.plot(X_plot, np.exp(tt_log_dens), color='green', label=treat)
                    ax6.legend()
                    ax6.set_title(treat + ' based T-squared statistics kernel estimated density')
                else:
                    ax6 = plt.subplot(326)
                    ax6.hlines(1, 1, 2)
                    ax6.set_title(treat + ' based T-squared statistics kernel estimated density not plotted')
                # plt.show()
                pdf.savefig()
                plt.close()
    plot_scatter_and_tsquared()
    return None


