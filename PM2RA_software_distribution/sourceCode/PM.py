# -*- coding: utf-8 -*-
from PM_gui import Ui_MainWindow
import sys
from PyQt5.QtGui import QIntValidator,QDoubleValidator,QDesktopServices
from PyQt5.QtCore import QDir,QUrl,QThread,pyqtSignal,QTimer
from PyQt5.QtWidgets import QAction, QFileDialog, QMainWindow,QApplication,QListWidgetItem
import os
import re
import pandas as pd
import platform
import multiprocessing
import load_save_project
import datetime
import time
import scatter_1D_exe
import network_plot_exe
import scatter_exe
import network_plot_ontaxa_exe
import scatter_1D_taxa_exe
import scatter_taxa_exe
import pm_score_nD_exe
from time import strftime
import pm_score_2D_exe
import subprocess
import warnings
warnings.filterwarnings("ignore")

path = getattr(sys, '_MEIPASS', os.getcwd())
os.chdir(path)


s = platform.uname()
os_p = s[0]

cpu_count = multiprocessing.cpu_count()


class query_window(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.lineEdit_4.setValidator(QIntValidator())
        self.ui.lineEdit_6.setValidator(QIntValidator())
        self.ui.lineEdit_10.setValidator(QIntValidator())
        self.ui.lineEdit_12.setValidator(QIntValidator())
        self.ui.lineEdit_9.setValidator(QDoubleValidator())
        self.ui.pushButton.clicked.connect(self.reset_all_input)
        self.ui.pushButton_14.clicked.connect(self.reset_all_input14)

        self.ui.progressBar.setDisabled(True)
        self.ui.progressBar_2.setDisabled(True)

        self.ui.checkBox.setChecked(True)
        self.ui.checkBox_4.setChecked(True)
        self.ui.checkBox_6.setChecked(True)

        self.ui.checkBox.clicked.connect(self.checkbox_pm_plot)
        self.ui.checkBox_4.clicked.connect(self.checkbox_pm_plot_4)
        self.ui.checkBox_6.clicked.connect(self.checkbox_pm_plot_6)

        self.ui.pushButton_3.clicked.connect(self.getFiles)
        self.ui.pushButton_4.clicked.connect(self.getFiles4)
        self.ui.pushButton_9.clicked.connect(self.getFiles9)
        self.ui.pushButton_16.clicked.connect(self.getFiles16)
        self.ui.pushButton_13.clicked.connect(self.getFiles13)
        self.ui.pushButton_18.clicked.connect(self.getFiles18)
        self.ui.pushButton_19.clicked.connect(self.getFiles19)

        self.ui.pushButton_17.clicked.connect(self.plot_pm_1d)
        self.ui.pushButton_2.clicked.connect(self.run_pm_2d)
        self.ui.pushButton_5.clicked.connect(self.run_pm_plot)
        self.ui.pushButton_8.clicked.connect(self.run_pm_plot_on_certain_taxa)
        self.ui.pushButton_12.clicked.connect(self.run_module_pm)
        self.ui.pushButton_11.clicked.connect(self.right_list_module)
        self.ui.pushButton_10.clicked.connect(self.left_list_module)
        self.ui.pushButton_6.clicked.connect(self.right_list_pm_network)
        self.ui.pushButton_7.clicked.connect(self.left_list_pm_network)
        self.ui.pushButton_15.clicked.connect(self.select_taxa)
        self.ui.comboBox_4.setCurrentIndex(1)

        self.ui.comboBox_11.setEnabled(False)
        self.ui.comboBox_25.setEnabled(False)
        self.ui.comboBox_10.setEnabled(False)

        self.ui.comboBox_7.currentIndexChanged.connect(self.mtchange)

        self.openmanual = QAction("Open Manual")
        self.openmanual.setShortcut('Ctrl+m')
        self.openmanual.triggered.connect(self.open_manual)
        self.ui.menu.addAction(self.openmanual)

        self.aboutpm = QAction("About PM2RA")
        self.aboutpm.setShortcut('Ctrl+a')
        self.aboutpm.triggered.connect(self.about_pm)
        self.ui.menu.addAction(self.aboutpm)

        self.exit = QAction("Exit Application")
        self.exit.setShortcut('Ctrl+q')
        self.exit.triggered.connect(self.exit_app)
        self.ui.menu.addAction(self.exit)

        self.ui.checkBox_8.setDisabled(True)
        self.ui.checkBox_7.setDisabled(True)
        self.ui.checkBox_6.setChecked(False)
        self.ui.checkBox_6.setDisabled(True)


    def checkbox_pm_plot_4(self):
        if self.ui.checkBox_4.isChecked():
            self.ui.comboBox_21.setEnabled(True)
        if not self.ui.checkBox_4.isChecked():
            self.ui.comboBox_21.setDisabled(True)

    def checkbox_pm_plot_6(self):
        if self.ui.checkBox_6.isChecked():
            self.ui.comboBox_18.setEnabled(True)
            self.ui.comboBox_19.setEnabled(True)
            self.ui.comboBox_24.setEnabled(True)
        if not self.ui.checkBox_6.isChecked():
            self.ui.comboBox_18.setDisabled(True)
            self.ui.comboBox_19.setDisabled(True)
            self.ui.comboBox_24.setDisabled(True)

    def checkbox_pm_plot(self):
        if self.ui.checkBox.isChecked():
            self.ui.comboBox_13.setEnabled(True)
            self.ui.comboBox_14.setEnabled(True)
            self.ui.comboBox_22.setEnabled(True)
        if not self.ui.checkBox.isChecked():
            self.ui.comboBox_13.setDisabled(True)
            self.ui.comboBox_14.setDisabled(True)
            self.ui.comboBox_22.setDisabled(True)


    def plot_pm_1d(self, i):
        self.ui.textEdit.setText('Logs...')
        QApplication.processEvents()
        tmpdatafilename = self.ui.lineEdit_11.text()
        try:
            project_dict = load_save_project.load_dict(tmpdatafilename)
            fileplace = project_dict.get('fileplace')
            if not os.path.exists(fileplace + '/PM_scores_1D.csv'):
                self.ui.textEdit.setText(
                    'Can not find PM_scores file in the project folder. Please check the input pm file.')
                return None
            if os.path.exists(fileplace + '/PM_scores_1D.csv'):
                pm_res_1D = pd.read_csv(fileplace + '/PM_scores_1D.csv', header=0, index_col=0)
                if 'taxa' not in pm_res_1D.columns:
                    self.ui.textEdit.setText(
                        'Can not parse PM_scores file in the project folder. Please check the input pm file.')
                    return None
                filter_ss=self.ui.comboBox_17.currentText()
                topnum=int( self.ui.lineEdit_18.text())
                plotscatter=self.ui.checkBox_4.isChecked()
                plotts = self.ui.checkBox_5.isChecked()
                scattertype = self.ui.comboBox_21.currentText()
                # print(filter_ss,topnum,plotscatter,plotts,scattertype)
                if (not plotscatter) and (not plotts):
                    self.ui.textEdit.setText(
                        'Please check the plot type.')
                    return None
                plottype=''
                if plotts and plotscatter:
                    scatter_1D_exe.plot_scatter_and_tsquared_exe(project_dict,filter_ss,topnum,scattertype)
                    plottype='Taxa differenct '+scattertype+', T squared statsistics plot'
                if plotscatter and (not plotts):
                    scatter_1D_exe.plot_scatter_exe(project_dict,filter_ss,topnum,scattertype)
                    plottype = 'Taxa differenct ' + scattertype
                if (not plotscatter) and plotts:
                    scatter_1D_exe.plot_t_squared_exe(project_dict, filter_ss, topnum, scattertype)
                    plottype = 'Taxa T squared statsistics plot'
                setting_str=project_dict.get('loggs')
                if os.path.exists(fileplace + '/taxa_difference_between_groups.pdf'):
                    sharp_str = '######################\r######################\r\r'
                    currenttime = str(datetime.datetime.now()) + '\r'
                    if os_p == 'Windows':
                        imageplace=fileplace + '\\taxa_difference_between_groups.pdf'
                    else:
                        imageplace = fileplace + '/taxa_difference_between_groups.pdf'
                    setting_str =setting_str+ sharp_str + currenttime + 'Taxa difference plotting:\r'+ \
                                '\tPlot the top ' +str(topnum)+' taxa with the highest scores\r'+ \
                                 '\tStatistically significant taxa filter: ' + str(filter_ss) + '\r' + \
                                 '\tImage type: ' + plottype+ '\r'+ \
                                  '\tImage file stored in ' +  imageplace+ '\r'
                    self.ui.textEdit.setText(setting_str)
                    project_dict['loggs']=setting_str
                    load_save_project.save_dict(
                        fileplace + '/project.pm', project_dict)
        except:
            self.ui.textEdit.setText('Project loads failed. Please check the input pm file.')
        return None

    def run_pm_plot(self):
        self.ui.textEdit.setText('Logs...')
        QApplication.processEvents()
        tmpdatafilename = self.ui.lineEdit_11.text()
        try:
            project_dict = load_save_project.load_dict(tmpdatafilename)
            fileplace = project_dict.get('fileplace')
            if not os.path.exists(fileplace + '/PM_scores.csv'):
                self.ui.textEdit.setText(
                    'Can not find PM_scores file in the project folder. Please check the input pm file.')
                return None
            if os.path.exists(fileplace + '/PM_scores.csv'):
                pm_res = pd.read_csv(fileplace + '/PM_scores.csv', header=0, index_col=0)
                if 'taxa1' not in pm_res.columns:
                    self.ui.textEdit.setText(
                        'Can not parse PM_scores file in the project folder. Please check the input pm file.')
                    return None
                plotpm = self.ui.checkBox.isChecked()
                plotscatter = self.ui.checkBox_2.isChecked()
                plotts= self.ui.checkBox_3.isChecked()
                topnum = int(self.ui.lineEdit_12.text())
                networktype= self.ui.comboBox_15.currentText()
                edgecolormap= self.ui.comboBox_13.currentText()
                nodecolormap = self.ui.comboBox_14.currentText()
                networklayout=self.ui.comboBox_22.currentText()
                filter_ss = self.ui.comboBox_16.currentText()
                if (not plotpm) and (not plotscatter)and (not plotts):
                    self.ui.textEdit.setText(
                        'Please check the plot type.')
                    return None
                if plotpm:
                    network_plot_exe.plotPMnetwork(project_dict, topnum, networktype, edgecolormap, nodecolormap, networklayout, filter_ss)
                    setting_str = project_dict.get('loggs')
                    if os.path.exists(fileplace + '/PM_network.pdf'):
                        sharp_str = '######################\r######################\r\r'
                        currenttime = str(datetime.datetime.now()) + '\r'
                        if os_p == 'Windows':
                            imageplace = fileplace + '\\PM_network.pdf'
                        else:
                            imageplace = fileplace + '/PM_network.pdf'
                        setting_str = setting_str + sharp_str + currenttime + 'PM network plotting:\r' + \
                                      '\tPlot the top ' + str(topnum) + ' edges with the highest '+networktype +'\r' + \
                                      '\tStatistically significant taxa filter: ' + str(filter_ss) + '\r' + \
                                      '\tEdge color map: ' + edgecolormap + '\r' + \
                                      '\tNode color map: ' + nodecolormap + '\r' + \
                                      '\tImage file stored in ' + imageplace + '\r'
                        self.ui.textEdit.setText(setting_str)
                        QApplication.processEvents()

                plottype = ''
                if plotts and plotscatter:
                    scatter_exe.plot_scatter_and_tsquared_exe(project_dict, topnum, networktype, filter_ss)
                    plottype = 'Taxa relaltionship difference and T squared statsistics plot'
                if plotscatter and (not plotts):
                    scatter_exe.plot_scatter_exe(project_dict, topnum, networktype, filter_ss)
                    plottype = 'Taxa relaltionship difference chart'
                if (not plotscatter) and plotts:
                    scatter_exe.plot_t_squared_exe(project_dict, topnum, networktype, filter_ss)
                    plottype = 'Taxa relationship T squared statsistics plot'
                setting_str = project_dict.get('loggs')
                if os.path.exists(fileplace + '/taxa_relationship_difference_between_groups.pdf'):
                    sharp_str = '######################\r######################\r\r'
                    currenttime = str(datetime.datetime.now()) + '\r'
                    if os_p == 'Windows':
                        imageplace = fileplace + '\\taxa_relationship_difference_between_groups.pdf'
                    else:
                        imageplace = fileplace + '/taxa_relationship_difference_between_groups.pdf'
                    setting_str = setting_str + sharp_str + currenttime + 'Taxa relationship difference plotting:\r' + \
                                  '\tPlot the taxa relationship in top ' + str(topnum) + ' edges with the highest ' + networktype + '\r' + \
                                  '\tStatistically significant taxa filter: ' + str(filter_ss) + '\r' + \
                                  '\tImage file stored in ' + imageplace + '\r'
                    self.ui.textEdit.setText(setting_str)
                    project_dict['loggs'] = setting_str
                    load_save_project.save_dict(
                        fileplace + '/project.pm', project_dict)

        except:
            self.ui.textEdit.setText('Project loads failed. Please check the input pm file.')
        return None

    def run_pm_plot_on_certain_taxa(self):
        self.ui.textEdit.setText('Logs...')
        QApplication.processEvents()
        tmpdatafilename = self.ui.lineEdit_11.text()
        try:
            project_dict = load_save_project.load_dict(tmpdatafilename)
            fileplace = project_dict.get('fileplace')
            if (not os.path.exists(fileplace + '/PM_scores.csv')) or (not os.path.exists(fileplace + '/PM_scores_1D.csv')):
                self.ui.textEdit.setText(
                    'Can not find PM_scores file in the project folder. Please check the input pm file.')
                return None
            if os.path.exists(fileplace + '/PM_scores.csv'):
                pm_res = pd.read_csv(fileplace + '/PM_scores.csv', header=0, index_col=0)
                if 'taxa1' not in pm_res.columns:
                    self.ui.textEdit.setText(
                        'Can not parse PM_scores file in the project folder. Please check the input pm file.')
                    return None
                plotpm = self.ui.checkBox_6.isChecked()
                plotscatter = self.ui.checkBox_7.isChecked()
                plotts = self.ui.checkBox_8.isChecked()
                itemcount = self.ui.listWidget_2.count()
                taxalist_plot=[]
                for itemindex in range(itemcount):
                    taxalist_plot.append(self.ui.listWidget_2.item(itemindex).text())
                if len(taxalist_plot)==0:
                    self.ui.textEdit.setText(
                        'Please select the taxa to be analyzed.')
                    return None
                networktype = self.ui.comboBox_23.currentText()
                edgecolormap = self.ui.comboBox_18.currentText()
                nodecolormap = self.ui.comboBox_19.currentText()
                networklayout = self.ui.comboBox_24.currentText()
                filter_ss = self.ui.comboBox_20.currentText()
                if (not plotpm) and (not plotscatter) and (not plotts):
                    self.ui.textEdit.setText(
                        'Please check the plot type.')
                    return None
                if plotpm:
                    if len(taxalist_plot)<3:
                        self.ui.textEdit.setText(
                            'The taxa number is too less for network plotting.')
                        return None
                    network_plot_ontaxa_exe.plotPMnetwork(project_dict, taxalist_plot, networktype, edgecolormap, nodecolormap,
                                                   networklayout, filter_ss)
                    setting_str = project_dict.get('loggs')
                    if os.path.exists(fileplace + '/PM_network_withintaxa.pdf'):
                        sharp_str = '######################\r######################\r\r'
                        currenttime = str(datetime.datetime.now()) + '\r'
                        if os_p == 'Windows':
                            imageplace = fileplace + '\\PM_network.pdf'
                        else:
                            imageplace = fileplace + '/PM_network.pdf'
                        setting_str = setting_str + sharp_str + currenttime + 'PM network plotting:\r' + \
                                      '\tPlot ' + '\t'.join(taxalist_plot) + ' \r' + \
                                      '\tStatistically significant taxa filter: ' + str(filter_ss) + '\r' + \
                                      '\tEdge color map: ' + edgecolormap + '\r' + \
                                      '\tNode color map: ' + nodecolormap + '\r' + \
                                      '\tImage file stored in ' + imageplace + '\r'
                        self.ui.textEdit.setText(setting_str)
                        QApplication.processEvents()

                plottype = ''
                if plotts and plotscatter:
                    if len(taxalist_plot)==1:
                        scatter_1D_taxa_exe.plot_scatter_and_tsquared_exe(project_dict,filter_ss,taxalist_plot,scattertype='Density Chart')
                    if len(taxalist_plot) > 1:
                        scatter_taxa_exe.plot_scatter_and_tsquared_exe(project_dict, taxalist_plot, networktype, filter_ss)
                    plottype = 'Taxa relaltionship difference and T squared statsistics plot'
                if plotscatter and (not plotts):
                    if len(taxalist_plot)==1:
                        scatter_1D_taxa_exe.plot_scatter_exe(project_dict,filter_ss,taxalist_plot,scattertype='Density Chart')
                    if len(taxalist_plot) > 1:
                        scatter_taxa_exe.plot_scatter_exe(project_dict, taxalist_plot, networktype, filter_ss)
                    plottype = 'Taxa relaltionship difference chart'
                if (not plotscatter) and plotts:
                    if len(taxalist_plot) == 1:
                        scatter_1D_taxa_exe.plot_t_squared_exe(project_dict, filter_ss, taxalist_plot,scattertype='Density Chart')
                    if len(taxalist_plot) > 1:
                        scatter_taxa_exe.plot_t_squared_exe(project_dict, taxalist_plot, networktype, filter_ss)
                    plottype = 'Taxa relationship T squared statsistics plot'
                setting_str = project_dict.get('loggs')
                if len(taxalist_plot) > 1:
                    if os.path.exists(fileplace + '/taxa_relationship_difference_between_groups.pdf'):
                        sharp_str = '######################\r######################\r\r'
                        currenttime = str(datetime.datetime.now()) + '\r'
                        if os_p == 'Windows':
                            imageplace = fileplace + '\\taxa_relationship_difference_between_groups_withintaxa.pdf'
                        else:
                            imageplace = fileplace + '/taxa_relationship_difference_between_groups_withintaxa.pdf'
                        setting_str = setting_str + sharp_str + currenttime + 'Taxa relationship difference plotting:\r' + \
                                      '\tPlot ' + '\t'.join(taxalist_plot) + ' \r' + \
                                      '\tStatistically significant taxa filter: ' + str(filter_ss) + '\r' + \
                                      '\tImage file stored in ' + imageplace + '\r'
                        self.ui.textEdit.setText(setting_str)
                        project_dict['loggs'] = setting_str
                        load_save_project.save_dict(
                            fileplace + '/project.pm', project_dict)
                if len(taxalist_plot) == 1:
                    if os.path.exists(fileplace + '/taxa_difference_between_groups_withintaxa.pdf'):
                        sharp_str = '######################\r######################\r\r'
                        currenttime = str(datetime.datetime.now()) + '\r'
                        if os_p == 'Windows':
                            imageplace = fileplace + '\\taxa_difference_between_groups_withintaxa.pdf'
                        else:
                            imageplace = fileplace + '/taxa_difference_between_groups_withintaxa.pdf'
                        setting_str = setting_str + sharp_str + currenttime + 'Taxa relationship difference plotting:\r' + \
                                      '\tPlot ' + '\t'.join(taxalist_plot) + ' \r' + \
                                      '\tStatistically significant taxa filter: ' + str(filter_ss) + '\r' + \
                                      '\tImage file stored in ' + imageplace + '\r'
                        self.ui.textEdit.setText(setting_str)
                        project_dict['loggs'] = setting_str
                        load_save_project.save_dict(
                            fileplace + '/project.pm', project_dict)

        except:
            self.ui.textEdit.setText('Project loads failed. Please check the input pm file.')
        return None

    def pm_run_pbar(self, i):
        self.ui.progressBar.setValue(i)
        QApplication.processEvents()

    def pm_run_pbar_1D(self, i):
        self.ui.progressBar_2.setValue(i)
        QApplication.processEvents()

    def getFiles(self):
        dig = QFileDialog()
        dig.setNameFilters(["csv data file(*.csv)"])
        dig.setFileMode(QFileDialog.ExistingFile)
        dig.setFilter(QDir.Files)
        if dig.exec_():
            filenames = dig.selectedFiles()
            if os_p == 'Windows':
                tmpdatafilename = filenames[0].replace('/', '\\')
            else:
                tmpdatafilename = filenames[0]
            self.ui.lineEdit_3.setText(tmpdatafilename)
    def getFiles4(self):
        dig = QFileDialog()
        dig.setNameFilters(["pm project file(*.pm)"])
        dig.setFileMode(QFileDialog.ExistingFile)
        dig.setFilter(QDir.Files)
        if dig.exec_():
            filenames = dig.selectedFiles()
            if os_p == 'Windows':
                tmpdatafilename = filenames[0].replace('/', '\\')
            else:
                tmpdatafilename = filenames[0]
            self.ui.lineEdit_11.setText(tmpdatafilename)
            try:

                project_dict = load_save_project.load_dict(tmpdatafilename)
                self.ui.lineEdit_2.setText(project_dict.get('fileplace'))
                self.ui.lineEdit_3.setText(project_dict.get('datafilename'))
                self.ui.lineEdit_4.setText(project_dict.get('minimum_taxa_detection_samples'))
                self.ui.lineEdit_9.setText(project_dict.get('minimum_taxa_median_abundance'))
                self.ui.lineEdit_6.setText(project_dict.get('minimum_coexist_taxa'))
                self.ui.comboBox_4.setCurrentIndex(project_dict.get('confidenceI'))
                self.ui.comboBox_8.setCurrentText(project_dict.get('fdr'))
                self.ui.comboBox_5.setCurrentText(project_dict.get('ilr'))
                self.ui.comboBox_6.setCurrentText(project_dict.get('autoband'))
                self.ui.comboBox_7.setCurrentText(project_dict.get('multithreading'))
                if project_dict.get('multithreading')=='OFF':
                    self.ui.label_17.setEnabled(False)
                    self.ui.lineEdit_10.setEnabled(False)
                if project_dict.get('multithreading')=='ON':
                    self.ui.label_17.setEnabled(True)
                    self.ui.lineEdit_10.setEnabled(True)
                self.ui.lineEdit_5.setText(project_dict.get('grouplabel'))
                self.ui.lineEdit_8.setText(project_dict.get('controllabel'))
                self.ui.lineEdit_7.setText(project_dict.get('treatlabel'))
                self.ui.lineEdit_10.setText(project_dict.get('threads'))
                self.ui.textEdit.setText(project_dict.get('loggs'))
                self.ui.listWidget.clear()
                self.ui.listWidget_2.clear()
                fileplace=project_dict.get('fileplace')
                if not os.path.exists(fileplace+'/PM_scores.csv'):
                    self.ui.textEdit.setText('Can not find PM_scores file in the project folder. Please check the input pm file.')
                if os.path.exists(fileplace + '/PM_scores.csv'):
                    pm_res=pd.read_csv(fileplace + '/PM_scores.csv',header=0,index_col=0)
                    try:
                        taxalist=list(set(list(pm_res['taxa1'])+list(pm_res['taxa2'])))
                        self.ui.listWidget.addItems(taxalist)
                    except:
                        self.ui.textEdit.setText(
                            'Can not parse PM_scores file in the project folder. Please check the input pm file.')
            except:
                self.ui.textEdit.setText('Project loads failed. Please check the input pm file.')

    def getFiles9(self):
        dig = QFileDialog()
        dig.setNameFilters(["csv data file(*.csv)"])
        dig.setFileMode(QFileDialog.ExistingFile)
        dig.setFilter(QDir.Files)
        if dig.exec_():
            filenames = dig.selectedFiles()
            if os_p == 'Windows':
                tmpdatafilename = filenames[0].replace('/', '\\')
            else:
                tmpdatafilename = filenames[0]
            self.ui.lineEdit_13.setText(tmpdatafilename)
            self.ui.listWidget_4.clear()
            self.ui.listWidget_3.clear()

    def getFiles13(self):
        directory1 = QFileDialog.getExistingDirectory(self, "select folder to storage all analysis results", "/")
        # print(directory1)
        if os_p == 'Windows':
            tmpdatafilename = directory1.replace('/', '\\')
        else:
            tmpdatafilename = directory1
        self.ui.lineEdit_2.setText(tmpdatafilename)

    def getFiles18(self):
        directory1 = QFileDialog.getExistingDirectory(self, "select folder to storage mudule analysis results", "/")
        # print(directory1)
        if os_p == 'Windows':
            tmpdatafilename = directory1.replace('/', '\\')
        else:
            tmpdatafilename = directory1
        self.ui.lineEdit_14.setText(tmpdatafilename)


    def getFiles19(self):
        self.ui.listWidget_3.reset()
        self.ui.listWidget_4.reset()
        dig = QFileDialog()
        dig.setNameFilters(["pm project file(*.pm)"])
        dig.setFileMode(QFileDialog.ExistingFile)
        dig.setFilter(QDir.Files)
        if dig.exec_():
            filenames = dig.selectedFiles()
            if os_p == 'Windows':
                tmpdatafilename = filenames[0].replace('/', '\\')
            else:
                tmpdatafilename = filenames[0]
            try:
                project_dict = load_save_project.load_dict(tmpdatafilename)
                self.ui.lineEdit_14.setText(project_dict.get('fileplace'))
                self.ui.lineEdit_13.setText(project_dict.get('datafilename'))
                self.ui.comboBox_10.setCurrentText(project_dict.get('ilr'))
                self.ui.comboBox_11.setCurrentText(project_dict.get('autoband'))
                self.ui.lineEdit_17.setText(project_dict.get('grouplabel'))
                self.ui.lineEdit_16.setText(project_dict.get('controllabel'))
                self.ui.lineEdit_15.setText(project_dict.get('treatlabel'))
                self.ui.textEdit.setText(project_dict.get('loggs'))
                self.ui.listWidget_3.clear()
                self.ui.listWidget_4.clear()
            except:
                self.ui.textEdit.setText('Project loads failed. Please check the input pm file.')

    def getFiles16(self):
        dig = QFileDialog()
        dig.setNameFilters(["pm project file(*.pm)"])
        dig.setFileMode(QFileDialog.ExistingFile)
        dig.setFilter(QDir.Files)
        if dig.exec_():
            filenames = dig.selectedFiles()
            if os_p == 'Windows':
                tmpdatafilename = filenames[0].replace('/', '\\')
            else:
                tmpdatafilename = filenames[0]
            try:
                project_dict = load_save_project.load_dict(tmpdatafilename)
                self.ui.lineEdit_2.setText(project_dict.get('fileplace'))
                self.ui.lineEdit_3.setText(project_dict.get('datafilename'))
                self.ui.lineEdit_4.setText(project_dict.get('minimum_taxa_detection_samples'))
                self.ui.lineEdit_9.setText(project_dict.get('minimum_taxa_median_abundance'))
                self.ui.lineEdit_6.setText(project_dict.get('minimum_coexist_taxa'))
                self.ui.comboBox_4.setCurrentIndex(project_dict.get('confidenceI'))
                self.ui.comboBox_8.setCurrentText(project_dict.get('fdr'))
                self.ui.comboBox_5.setCurrentText(project_dict.get('ilr'))
                self.ui.comboBox_6.setCurrentText(project_dict.get('autoband'))
                self.ui.comboBox_7.setCurrentText(project_dict.get('multithreading'))
                if project_dict.get('multithreading')=='OFF':
                    self.ui.label_17.setEnabled(False)
                    self.ui.lineEdit_10.setEnabled(False)
                if project_dict.get('multithreading')=='ON':
                    self.ui.label_17.setEnabled(True)
                    self.ui.lineEdit_10.setEnabled(True)
                self.ui.lineEdit_5.setText(project_dict.get('grouplabel'))
                self.ui.lineEdit_8.setText(project_dict.get('controllabel'))
                self.ui.lineEdit_7.setText(project_dict.get('treatlabel'))
                self.ui.lineEdit_10.setText(project_dict.get('threads'))
                self.ui.textEdit.setText(project_dict.get('loggs'))
            except:
                self.ui.textEdit.setText('Project loads failed. Please check the input pm file.')

    def reset_all_input(self):

        self.ui.progressBar.reset()
        self.ui.progressBar.setEnabled(False)
        self.ui.progressBar_2.reset()
        self.ui.progressBar_2.setEnabled(False)
        self.ui.comboBox_4.setCurrentIndex(1)
        self.ui.comboBox_8.setCurrentIndex(0)
        self.ui.comboBox_6.setCurrentIndex(0)
        self.ui.comboBox_7.setCurrentIndex(0)

        self.ui.lineEdit_2.setText('your project name')
        self.ui.lineEdit_3.setText('D:\\data\\abundance_to_be_analyzed.csv')
        self.ui.lineEdit_4.setText('5')
        self.ui.lineEdit_5.setText('condition')
        self.ui.lineEdit_6.setText('5')
        self.ui.lineEdit_7.setText('treat')
        self.ui.lineEdit_8.setText('control')
        self.ui.lineEdit_9.setText('0')
        self.ui.lineEdit_10.setText('4')
        self.ui.textEdit.setText('Logs..')


    def select_taxa(self):
        # to be edited
        if os_p == 'Windows':
            datafilename = self.ui.lineEdit_13.text().replace('/', '\\')
        else:
            datafilename = self.ui.lineEdit_13.text()
        if not os.path.exists(datafilename):
            setting_str='input data file not exists. Please check your input configuration.Quit the analysis'
            self.ui.textEdit.setText(setting_str)
            return None
        try:
            demodata = pd.read_csv(datafilename, header=0, index_col=False)
        except:
            setting_str = 'Can not parse the input data as csv data file.Please check your input configuration. Quit the analysis.'
            self.ui.textEdit.setText(setting_str)
            return None
        condition=self.ui.lineEdit_17.text()
        taxaname = demodata.columns
        if condition not in taxaname:
            setting_str = 'The group indicator column are not found in the datafile. Please check your input configuration. Quit the analysis.'
            self.ui.textEdit.setText(setting_str)
            return None
        genuslist=[]
        for g in taxaname:
            if g != condition:
                tmpgdata=demodata[g]
                obs=sum(tmpgdata>0)
                if obs>12:
                    genuslist.append(g)
        self.ui.listWidget_4.clear()
        self.ui.listWidget_3.clear()
        self.ui.listWidget_3.addItems(genuslist)

        return None

    def reset_all_input14(self):
        # to be edited
        self.ui.listWidget_4.clear()
        self.ui.listWidget_3.clear()

        self.ui.comboBox_10.setCurrentIndex(0)
        self.ui.comboBox_11.setCurrentIndex(0)
        self.ui.comboBox_25.setCurrentIndex(0)

        self.ui.lineEdit_14.setText('D:\\data')
        self.ui.lineEdit_13.setText('D:\\data\\abundance_to_be_analyzed.csv')
        self.ui.lineEdit_17.setText('condition')
        self.ui.lineEdit_16.setText('control')
        self.ui.lineEdit_15.setText('treat')

        self.ui.textEdit.setText('Logs..')

        return None


    def mtchange(self):
        if self.ui.comboBox_7.currentText() =='ON':
            self.ui.label_17.setEnabled(True)
            self.ui.lineEdit_10.setEnabled(True)
        if self.ui.comboBox_7.currentText() =='OFF':
            self.ui.label_17.setEnabled(False)
            self.ui.lineEdit_10.setEnabled(False)
        return None

    def right_list_module(self):
        try:
            checked=self.ui.listWidget_3.currentItem().text()
            listItems = self.ui.listWidget_3.selectedItems()
            for item in listItems:
                self.ui.listWidget_3.takeItem(self.ui.listWidget_3.row(item))
            listWidgetItem = QListWidgetItem(checked)
            self.ui.listWidget_4.addItem(listWidgetItem)
            itemcount=self.ui.listWidget_4.count()
            if itemcount>0:
                self.ui.comboBox_11.setEnabled(True)
                self.ui.comboBox_25.setEnabled(True)
                self.ui.comboBox_10.setEnabled(True)
        except:
            pass
        return None

    def left_list_module(self):
        try:
            checked=self.ui.listWidget_4.currentItem().text()
            listItems = self.ui.listWidget_4.selectedItems()
            for item in listItems:
                self.ui.listWidget_4.takeItem(self.ui.listWidget_4.row(item))
            listWidgetItem = QListWidgetItem(checked)
            self.ui.listWidget_3.addItem(listWidgetItem)
            itemcount=self.ui.listWidget_4.count()
            if itemcount==0:
                self.ui.comboBox_11.setEnabled(False)
                self.ui.comboBox_25.setEnabled(False)
                self.ui.comboBox_10.setEnabled(False)
        except:
            pass
        return None
    def right_list_pm_network(self):
        # to be edited
        try:
            checked=self.ui.listWidget.currentItem().text()
            listItems = self.ui.listWidget.selectedItems()
            for item in listItems:
                self.ui.listWidget.takeItem(self.ui.listWidget.row(item))
            listWidgetItem = QListWidgetItem(checked)
            self.ui.listWidget_2.addItem(listWidgetItem)
            listItemsCount = self.ui.listWidget_2.count()
            if listItemsCount>0:
                self.ui.checkBox_8.setEnabled(True)
                self.ui.checkBox_7.setEnabled(True)
            if listItemsCount>2:
                self.ui.checkBox_6.setEnabled(True)
        except:
            pass
        return None

    def left_list_pm_network(self):
        # to be edited
        try:
            checked=self.ui.listWidget_2.currentItem().text()
            listItems = self.ui.listWidget_2.selectedItems()
            for item in listItems:
                self.ui.listWidget_2.takeItem(self.ui.listWidget_2.row(item))
            listWidgetItem = QListWidgetItem(checked)
            self.ui.listWidget.addItem(listWidgetItem)
            listItemsCount = self.ui.listWidget_2.count()
            if listItemsCount==0:
                self.ui.checkBox_8.setDisabled(True)
                self.ui.checkBox_7.setDisabled(True)
            if listItemsCount<=2:
                self.ui.checkBox_6.setChecked(False)
                self.ui.checkBox_6.setDisabled(True)

        except:
            pass
        return None
    def about_pm(self):
        file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "source/About_PM.html"))
        local_url = QUrl.fromLocalFile(file_path)
        QDesktopServices.openUrl(local_url)

    def open_manual(self):
        file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "source/manual.html"))
        local_url = QUrl.fromLocalFile(file_path)
        QDesktopServices.openUrl(local_url)

    def exit_app(self):
        self.close()

    def run_pm_2d(self):
        if os_p == 'Windows':
            datafilename = self.ui.lineEdit_3.text().replace('/', '\\')
        else:
            datafilename = self.ui.lineEdit_3.text()

        fileplace = self.ui.lineEdit_2.text()

        sharp_str = '######################\r######################\r\r'
        currenttime = str(datetime.datetime.now()) + '\r'
        setting_str = sharp_str+currenttime+'Setting:\r'
        setting_str = setting_str+'\tProject analysis results storage folder: '+self.ui.lineEdit_2.text()+'\r' + \
                      '\tData file: ' + datafilename + '\r'+ \
                      '\tMinimum taxa detection samples in each group: ' + self.ui.lineEdit_4.text() + '\r' + \
                      '\tMinimum taxa median abundance in control group: ' + self.ui.lineEdit_9.text() + '\r' + \
                      '\tMinimum coexist taxa number: ' + self.ui.lineEdit_6.text() + '\r' + \
                    '\tConfidence Interval: '+self.ui.comboBox_4.currentText()+'\r' + \
                    '\tFDR adjustment: ' + self.ui.comboBox_8.currentText() + '\r' + \
                      '\tilr transform: ' + self.ui.comboBox_5.currentText() + '\r' + \
                      '\tAutomated Best Bandwidth:: ' + self.ui.comboBox_6.currentText() + '\r' +\
                    '\tColumn name of group indicator: ' + self.ui.lineEdit_5.text() + '\r' + \
                    '\tLabel of control group: ' + self.ui.lineEdit_8.text() + '\r' + \
                    '\tLabel of treatment group: ' + self.ui.lineEdit_7.text() + '\r'
        if self.ui.comboBox_7.currentText()=='ON':
            setting_str=setting_str+'\tMulti threading: '+self.ui.comboBox_7.currentText()+'\r' + \
                        '\tThreads: ' + self.ui.lineEdit_10.text() + '\r'
        if self.ui.comboBox_7.currentText() == 'OFF':
            setting_str=setting_str+'\tMulti threading: '+self.ui.comboBox_7.currentText()+ '\r'
        QApplication.processEvents()
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        isExists = os.path.exists(self.ui.lineEdit_3.text())
        if not isExists:
            setting_str = setting_str + '\tThe data file not existed. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if isExists:
            try:
                relativedata = pd.read_csv(datafilename, header=0, index_col=False)
            except:
                setting_str = setting_str + \
                    '\tCan not read into the data file. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None
        # check the parameter
        try:
            group = relativedata[self.ui.lineEdit_5.text()]
        except:
            setting_str = setting_str + '\tCan not find the group indicator column in the data file. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        controlnum = len(
            relativedata[relativedata[self.ui.lineEdit_5.text()] == self.ui.lineEdit_8.text()])
        treatnum = len(
            relativedata[relativedata[self.ui.lineEdit_5.text()] == self.ui.lineEdit_7.text()])
        if controlnum == 0:
            setting_str = setting_str + '\tData file do not contain the control group data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum == 0:
            setting_str = setting_str + '\tData file do not contain the treat group data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if controlnum < 5:
            setting_str = setting_str + \
                '\tThe sample size of control group data is to small. It is preferably larger than 5. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum < 5:
            setting_str = setting_str + \
                '\tThe sample size of control group data is to small. It is preferably larger than 5. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        minimum_taxa_detection_num = int(self.ui.lineEdit_4.text())
        if minimum_taxa_detection_num < 2:
            setting_str = setting_str + \
                '\tThe minimum taxa detection num should be larger than 2. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if controlnum < minimum_taxa_detection_num:
            setting_str = setting_str + \
                '\tThe minimum taxa detection num should be smaller than the sample size of control group ('+str(
                    controlnum)+'). Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum < minimum_taxa_detection_num:
            setting_str = setting_str + \
                '\tThe minimum taxa detection num should be smaller than the sample size of treat group ('+str(
                    treatnum)+'). Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        minimum_coexist_taxa_num = int(self.ui.lineEdit_6.text())
        if minimum_coexist_taxa_num < 0:
            setting_str = setting_str + \
                '\tThe minimum coexist taxa num should be positive. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if controlnum < minimum_coexist_taxa_num:
            setting_str = setting_str + \
                '\tThe minimum coexist taxa num  should be smaller than the sample size of control group ('+str(
                    controlnum)+'). Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum < minimum_coexist_taxa_num:
            setting_str = setting_str + \
                '\tThe minimum coexist taxa num should be smaller than the sample size of treat group ('+str(
                    treatnum)+'). Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        tmptaxalist = list(relativedata.columns)
        taxalist = []
        all_float_dtype=True
        countneg=0
        for g in tmptaxalist:
            if g != self.ui.lineEdit_5.text():
                taxalist.append(g)
                if not pd.api.types.is_numeric_dtype(relativedata[g]):
                    all_float_dtype=False
                else:
                    countneg=countneg+sum(relativedata[g]<0)
        if countneg>0:
            setting_str = setting_str + \
                '\tThe data file contains negative data. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if not all_float_dtype:
            setting_str = setting_str + \
                '\tSome columns of the data file are not abundance data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if len(taxalist) == 0:
            setting_str = setting_str + \
                '\tNo taxa are found in the data file. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        setting_str = setting_str + '\tThe data file contains ' + \
            str(len(taxalist))+' taxa' + '\r'
        self.ui.textEdit.setText(setting_str)
        setting_str = setting_str + '\tThe '+self.ui.lineEdit_8.text()+' group contains ' + \
            str(controlnum) + ' samples' + '\r'
        self.ui.textEdit.setText(setting_str)
        setting_str = setting_str + '\tThe ' + self.ui.lineEdit_7.text() + ' group contains ' + \
            str(treatnum) + ' samples' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        # print(self.ui.lineEdit.text())
        # print(self.ui.comboBox_2.currentText())
        minimum_taxa_abundance_median=float(self.ui.lineEdit_9.text())
        if minimum_taxa_abundance_median < 0:
            setting_str = setting_str + \
                '\tThe parameter minimum taxa median abundance in control group should be non-negative. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        if self.ui.comboBox_7.currentText()=='ON':
            if int( self.ui.lineEdit_10.text())>cpu_count:
                setting_str = setting_str + \
                              '\tThe threads should be set less than counts of logical processors ('+ str(cpu_count)+'). Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None

        isExists = os.path.exists(fileplace)
        if not isExists:
            setting_str = setting_str+'\tThe project analysis results storage folder not exist. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        try:
            timestr=strftime("%Y%m%d%H%M%S")
            if os_p == 'Windows':
                logs_place = fileplace+"\\logs"+timestr
            else:
                logs_place = fileplace + "/logs"+timestr
            if not os.path.exists(logs_place):
                os.makedirs(logs_place)
            setting_str = setting_str + '\tProject run log are stored in "' + \
                logs_place + '"\r'
            self.ui.textEdit.setText(setting_str)
        except:
            setting_str = setting_str + \
                '\tCreating project running logs storage folder ("' + \
                logs_place+') failed. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        currenttime = str(datetime.datetime.now()) + '\r'
        setting_str = setting_str + sharp_str+currenttime + 'Start analysis' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()

        self.ui.progressBar.reset()
        self.ui.progressBar.setEnabled(True)
        self.ui.progressBar_2.reset()
        self.ui.progressBar_2.setEnabled(True)
        QApplication.processEvents()

        threads_count = int(self.ui.lineEdit_10.text())
        control_label = self.ui.lineEdit_8.text()
        treat_label = self.ui.lineEdit_7.text()
        group_label = self.ui.lineEdit_5.text()
        fdr = self.ui.comboBox_8.currentText()
        ilr = self.ui.comboBox_5.currentText()
        confidence_level = [0.99, 0.95,0.90][self.ui.comboBox_4.currentIndex()]
        autoBand = self.ui.comboBox_6.currentText()
        kernel_str=self.ui.comboBox_12.currentText()
        # fileplace, logs_place, datafilename, minimum_coexist_taxa_num, minimum_taxa_detection_num, minimum_taxa_abundance_median,

        if self.ui.comboBox_7.currentText() == 'ON':
            # this cmd_str should be modified to ./python.exe in the distributed exe
            cmd_str="python pm_score_2D_multi_thread_exe.py "+fileplace+' '+logs_place+' '+datafilename+' '+ \
                      str(minimum_coexist_taxa_num)+' '+str(minimum_taxa_detection_num)+ \
                      ' ' + str(minimum_taxa_abundance_median) +' '+control_label+' '+treat_label \
                    +' '+group_label+' '+str(threads_count)+' '+fdr+' '+ilr+' '+str(confidence_level)+' '+autoBand+' '+kernel_str
            print(cmd_str)
            self.ui.thread_1 = Worker(logs_place,threads_count)
            self.ui.thread_1.progressBarValue.connect(self.pm_run_pbar)
            self.ui.thread_1.start()
            QApplication.processEvents()
            self.ui.thread_3 = Worker3(logs_place,threads_count)
            self.ui.thread_3.progressBarValue.connect(self.pm_run_pbar_1D)
            self.ui.thread_3.start()
            QApplication.processEvents()
            self.ui.thread_2 = Worker2(cmd_str)
            self.ui.thread_2.start()
            self.ui.thread_2.trigger.connect(self.finish_pm)
            QApplication.processEvents()

        if self.ui.comboBox_7.currentText() == 'OFF':
            cmd_str="python pm_score_2D_exe.py "+fileplace+' '+logs_place+' '+datafilename+' '+ \
                      str(minimum_coexist_taxa_num)+' '+str(minimum_taxa_detection_num)+ \
                      ' ' + str(minimum_taxa_abundance_median) +' '+control_label+' '+treat_label \
                    +' '+group_label+' '+fdr+' '+ilr+' '+str(confidence_level)+' '+autoBand+' '+kernel_str
            print(cmd_str)
            self.ui.thread_1 = Worker(logs_place,0)
            self.ui.thread_1.progressBarValue.connect(self.pm_run_pbar)
            self.ui.thread_1.start()
            QApplication.processEvents()
            self.ui.thread_3 = Worker3(logs_place,threads_count)
            self.ui.thread_3.progressBarValue.connect(self.pm_run_pbar_1D)
            self.ui.thread_3.start()
            QApplication.processEvents()
            self.ui.thread_4 = Worker4(cmd_str,fileplace,logs_place,datafilename,minimum_coexist_taxa_num,minimum_taxa_detection_num,minimum_taxa_abundance_median,control_label,treat_label,group_label,fdr,ilr,confidence_level,autoBand,kernel_str)
            self.ui.thread_4.start()
            self.ui.thread_4.trigger.connect(self.finish_pm)
            QApplication.processEvents()



    def run_module_pm(self):
        itemcount=self.ui.listWidget_4.count()
        if itemcount==0:
            setting_str = '\tPlease select the taxa to be analyzed. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        if os_p == 'Windows':
            datafilename = self.ui.lineEdit_13.text().replace('/', '\\')
        else:
            datafilename = self.ui.lineEdit_13.text()

        fileplace=self.ui.lineEdit_14.text()
        sharp_str = '######################\r######################\r\r'
        currenttime = str(datetime.datetime.now()) + '\r'
        setting_str = sharp_str + currenttime + 'Setting:\r'
        setting_str = setting_str + '\tModule analysis results storage folder: ' + self.ui.lineEdit_14.text() + '\r' + \
                      '\tData file: ' + datafilename + '\r' + \
                      '\tKernel density estimation: ' + self.ui.comboBox_25.currentText() + '\r' + \
                      '\tilr transform: ' + self.ui.comboBox_10.currentText() + '\r' + \
                      '\tAutomated Best Bandwidth:: ' + self.ui.comboBox_11.currentText() + '\r' + \
                      '\tColumn name of group indicator: ' + self.ui.lineEdit_17.text() + '\r' + \
                      '\tLabel of control group: ' + self.ui.lineEdit_16.text() + '\r' + \
                      '\tLabel of treatment group: ' + self.ui.lineEdit_15.text() + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        isExists = os.path.exists(self.ui.lineEdit_13.text())
        if not isExists:
            setting_str = setting_str + '\tThe data file not existed. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if isExists:
            try:
                relativedata = pd.read_csv(datafilename, header=0, index_col=False)
            except:
                setting_str = setting_str + \
                              '\tCan not read into the data file. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None
        # check the parameter
        try:
            group = relativedata[self.ui.lineEdit_17.text()]
        except:
            setting_str = setting_str + '\tCan not find the group indicator column in the data file. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        controlnum = len(
            relativedata[relativedata[self.ui.lineEdit_17.text()] == self.ui.lineEdit_16.text()])
        treatnum = len(
            relativedata[relativedata[self.ui.lineEdit_17.text()] == self.ui.lineEdit_15.text()])
        if controlnum == 0:
            setting_str = setting_str + '\tData file do not contain the control group data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum == 0:
            setting_str = setting_str + '\tData file do not contain the treat group data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if controlnum < 5:
            setting_str = setting_str + \
                          '\tThe sample size of control group data is to small. It is preferably larger than 5. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if treatnum < 5:
            setting_str = setting_str + \
                          '\tThe sample size of control group data is to small. It is preferably larger than 5. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        tmptaxalist = list(relativedata.columns)
        taxalist = []
        all_float_dtype = True
        countneg = 0
        for g in tmptaxalist:
            if g != self.ui.lineEdit_17.text():
                taxalist.append(g)
                if not pd.api.types.is_numeric_dtype(relativedata[g]):
                    all_float_dtype = False
                else:
                    countneg = countneg + sum(relativedata[g] < 0)
        if countneg > 0:
            setting_str = setting_str + \
                          '\tThe data file contains negative data. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if not all_float_dtype:
            setting_str = setting_str + \
                          '\tSome columns of the data file are not abundance data. Please check you data file. Make sure it follows the data template. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if len(taxalist) == 0:
            setting_str = setting_str + \
                          '\tNo taxa are found in the data file. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        setting_str = setting_str + '\tThe data file contains ' + \
                      str(len(taxalist)) + ' taxa' + '\r'
        self.ui.textEdit.setText(setting_str)
        setting_str = setting_str + '\tThe ' + self.ui.lineEdit_8.text() + ' group contains ' + \
                      str(controlnum) + ' samples' + '\r'
        self.ui.textEdit.setText(setting_str)
        setting_str = setting_str + '\tThe ' + self.ui.lineEdit_7.text() + ' group contains ' + \
                      str(treatnum) + ' samples' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        # print(self.ui.lineEdit.text())
        # print(self.ui.comboBox_2.currentText())
        minimum_taxa_abundance_median = float(self.ui.lineEdit_9.text())
        if minimum_taxa_abundance_median < 0:
            setting_str = setting_str + \
                          '\tThe parameter minimum taxa median abundance in control group should be non-negative. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        if self.ui.comboBox_7.currentText() == 'ON':
            if int(self.ui.lineEdit_10.text()) > cpu_count:
                setting_str = setting_str + \
                              '\tThe threads should be set less than counts of logical processors (' + str(
                    cpu_count) + '). Quit the analysis.' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None

        isExists = os.path.exists(fileplace)
        if not isExists:
            setting_str = setting_str + '\tModule analysis results storage folder not exist. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None

        currenttime = str(datetime.datetime.now()) + '\r'
        setting_str = setting_str + sharp_str + currenttime + 'Start analysis' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()

        control_label = self.ui.lineEdit_16.text()
        treat_label = self.ui.lineEdit_15.text()
        group_label = self.ui.lineEdit_17.text()
        ilr = self.ui.comboBox_10.currentText()
        autoBand = self.ui.comboBox_11.currentText()
        kernel_str = self.ui.comboBox_25.currentText()

        taxalist_plot = []
        for itemindex in range(itemcount):
            taxalist_plot.append(self.ui.listWidget_4.item(itemindex).text())
        # fileplace, logs_place, datafilename, minimum_coexist_taxa_num, minimum_taxa_detection_num, minimum_taxa_abundance_median,

        pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde=pm_score_nD_exe.cal_pm_score_nd(datafilename, control_label, treat_label, group_label, ilr, autoBand, kernel_str, taxalist_plot)
        if pm_score is None:
            setting_str = setting_str + \
                          '\tThe sample size for the selected taxa is smaller than 6. Quit the analysis.' + '\r'
            self.ui.textEdit.setText(setting_str)
            return None
        if pm_score is not None:
            pm_score_nD_exe.plot_pm_score_nd(fileplace,datafilename,control_label,treat_label,group_label,taxalist_plot,pm_score,pvalue,cc,ct,tt,tc,cc_kde,ct_kde,tt_kde,tc_kde)
            if os.path.exists(fileplace+'\module_difference_between_groups.pdf'):
                setting_str = setting_str + \
                              '\tThe generated image file is stored in .' + fileplace+'\module_difference_between_groups.pdf' + '\r' \
                              '\tPM score is  ' + str(round(pm_score,3)) + '\r'\
                              '\tP value is  ' + str(round(pvalue, 3)) + '\r'
                self.ui.textEdit.setText(setting_str)
                return None
            if not os.path.exists(fileplace+'\module_difference_between_groups.pdf'):
                setting_str = setting_str + \
                              '\tImage generaation failed. Please check your input datafile' + '\r'
                self.ui.textEdit.setText(setting_str)
                return None



    def finish_pm(self,cmd_str):
        cmd_str=cmd_str.split(' ')
        fileplace=cmd_str[2]
        logs_place=cmd_str[3]
        # print(logs_place)
        sucess=False
        if os.path.exists(logs_place + '/task_finished.json'):
            try:
                finishedtask = load_save_project.load_dict(logs_place + '/task_finished.json')
                finishedtask = float(finishedtask.get('success'))
                sucess=finishedtask
                # print(sucess)
            except:
                pass
        currenttime = str(datetime.datetime.now()) + '\r'
        sharp_str = '######################\r######################\r\r'
        setting_str=self.ui.textEdit.toPlainText()
        if sucess:
            setting_str = setting_str + sharp_str+currenttime + 'Analysis finished' + '\r'
        if not sucess:
            setting_str = setting_str + sharp_str + currenttime + 'Analysis failed, please check your input configuration.' + '\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()
        project_dict = {
                        'datafilename': self.ui.lineEdit_3.text(),
                        'minimum_taxa_detection_samples': self.ui.lineEdit_4.text(),
                        'minimum_taxa_median_abundance': self.ui.lineEdit_9.text(),
                        'minimum_coexist_taxa': self.ui.lineEdit_6.text() ,
                        'confidenceI': self.ui.comboBox_4.currentIndex(),
                        'fdr': self.ui.comboBox_8.currentText(),
                        'ilr': self.ui.comboBox_5.currentText(),
                        'autoband': self.ui.comboBox_6.currentText(),
                        'multithreading': self.ui.comboBox_7.currentText(),
                        'grouplabel': self.ui.lineEdit_5.text(),
                        'controllabel': self.ui.lineEdit_8.text(),
                        'treatlabel': self.ui.lineEdit_7.text(),
                        'threads': self.ui.lineEdit_10.text(),
                        'fileplace': fileplace,
                        'logs_place': logs_place,
                        'loggs': setting_str

                        }
        load_save_project.save_dict(
            fileplace + '/project.pm', project_dict)
        if os_p == 'Windows':
            logfileplace = (
                    fileplace + '/project.pm').replace('/', '\\')
        else:
            logfileplace = (fileplace + '/project.pm')
        setting_str = setting_str + 'The project configuration and logs are stored in "' + logfileplace + '."\r'
        self.ui.textEdit.setText(setting_str)
        QApplication.processEvents()



class Worker(QThread):
    progressBarValue = pyqtSignal(int)  # 更新进度条
    def __init__(self,logs_place,threads_count):
        super(Worker, self).__init__()
        self.logs_place = logs_place
        self.threads_count = threads_count
    def run(self):
        count=0
        while (1):
            time.sleep(1)
            # print(self.logs_place)
            try:
                totaltask = load_save_project.load_dict(self.logs_place + '/totaltask.json')
                totaltask = float(totaltask.get('totalcomb'))
                if self.threads_count==0:
                    finishedtask=load_save_project.load_dict(self.logs_place+'/0_task.json')
                    finishedtask = float(finishedtask.get('hasfinished'))
                    count=int(finishedtask/totaltask*100)
                if self.threads_count!=0:
                    allfinished=0
                    for i in range(self.threads_count):
                        finishedtask=0
                        try:
                            finishedtask=load_save_project.load_dict(self.logs_place+'/'+str(i)+'_task.json')
                            finishedtask = float(finishedtask.get('hasfinished'))
                        except:
                            pass
                        allfinished=allfinished+finishedtask
                    count=int(allfinished/totaltask*100)
            except:
                pass
            if os.path.exists(self.logs_place+'/task_finished.json'):
                count=100
            # print(count)
            self.progressBarValue.emit(count)
            if count>99:
                break


class Worker3(QThread):
    progressBarValue = pyqtSignal(int)  # 更新进度条
    def __init__(self,logs_place,threads_count):
        super(Worker3, self).__init__()
        self.logs_place = logs_place
        self.threads_count = threads_count
    def run(self):
        count=0
        while (1):
            time.sleep(1)
            # print(self.logs_place)
            try:
                totaltask = load_save_project.load_dict(self.logs_place + '/totaltask_1D.json')
                totaltask = float(totaltask.get('totalcomb'))
                if self.threads_count==0:
                    finishedtask=load_save_project.load_dict(self.logs_place+'/0_task_1D.json')
                    finishedtask = float(finishedtask.get('hasfinished'))
                    count=int(finishedtask/totaltask*100)
                if self.threads_count!=0:
                    allfinished=0
                    for i in range(self.threads_count):
                        finishedtask=0
                        try:
                            finishedtask=load_save_project.load_dict(self.logs_place+'/'+str(i)+'_task_1D.json')
                            finishedtask = float(finishedtask.get('hasfinished'))
                        except:
                            pass
                        allfinished=allfinished+finishedtask
                    count=int(allfinished/totaltask*100)
            except:
                pass
            if os.path.exists(self.logs_place+'/task_finished_1D.json'):
                count=100
            # print(count)
            self.progressBarValue.emit(count)
            if count>99:
                break

class Worker2(QThread):
    trigger = pyqtSignal(str)  # 更新进度条
    def __init__(self,cmd_str):
        super(Worker2, self).__init__()
        self.cmd_str = cmd_str
    def run(self):
        # subprocess.call(self.cmd_str,shell=False)
        subprocess.call(self.cmd_str,shell=True)
        self.trigger.emit(self.cmd_str)#

class Worker4(QThread):
    trigger = pyqtSignal(str)  # 更新进度条
    def __init__(self,cmd_str,fileplace,logs_place,datafilename,minimum_coexist_taxa_number,minimum_taxa_number,minimum_taxa_prevlance,control,treat,condition,fdr,ilr,confidence_level,wetherFindBestBandwidth,kernel_str):
        super(Worker4, self).__init__()
        self.cmd_str=cmd_str
        self.fileplace = fileplace
        self.logs_place = logs_place
        self.datafilename = datafilename
        self.minimum_coexist_taxa_number = minimum_coexist_taxa_number
        self.minimum_taxa_number = minimum_taxa_number
        self.minimum_taxa_prevlance = minimum_taxa_prevlance
        self.control = control
        self.treat = treat
        self.condition = condition
        self.fdr = fdr
        self.ilr = ilr
        self.confidence_level = confidence_level
        self.wetherFindBestBandwidth = wetherFindBestBandwidth
        self.kernel_str = kernel_str

    def run(self):
        # subprocess.call(self.cmd_str,shell=False)
        pm_score_2D_exe.pm_score_2d(self.fileplace, self.logs_place, self.datafilename, self.minimum_coexist_taxa_number, self.minimum_taxa_number, self.minimum_taxa_prevlance, self.control, self.treat, self.condition, self.fdr, self.ilr, self.confidence_level, self.wetherFindBestBandwidth, self.kernel_str)
        self.trigger.emit(self.cmd_str)#


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = query_window()
    window.show()
    sys.exit(app.exec_())
