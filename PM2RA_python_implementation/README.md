# Character User Interface PM2RA@XYL_Lab

This is the Character User Interface version (Python implemented) of PM2RA@XYL_Lab.
The multi-threading version and single-threading version are both provided.

## Installation recommendations
It is recommended to install CUI-PM2RA@XYL_Lab as follows
### using conda
```python
conda create -n PM2RA python=3.6
source activate PM2RA
conda install numpy==1.19.1
conda install pandas==1.1.0
conda install scipy==1.5.2
conda install scikit-learn==0.23.2
```
### using pipenv
```python
pipenv install --python 3.6.7
pipenv shell
pip install numpy==1.19.1
pip install pandas==1.1.0
pip install scipy==1.5.2
pip install scikit-learn==0.23.2
```
Note: only python 3.6.7 was tested.


## Instructions for use
### Arguments for single-processing version PM2RA
The *PM_Analysis.py* is the entrance of the single-processing analysis. It takes following arguments.
> 1. fileplace. The place to store analysis result.
> 2. logs_place. The place to store analysis process logs. It is recommend using subfolder under fileplace. Please make sure the fileplace folder and logs_place folder do exist in your system. 
> 3. datafilename. The absolute path of the csv file to be analyzed. Both the datafile name and its relative or absolute path should be provided.
> 4. minimum_coexist_taxa_number. PM2RA analyzes the taxa relationship alternations between groups. Only samples with paired taxa are pipelined into the analysis. Take taxa1 and taxa2 as an example, only taxa1's abundance and taxa2's abundance are both larger than 0, the sample will be into the analysis. This parameter minimum_coexist_taxa_number specifies the minimum sample size of each group with detected taxa1 and taxa2. The relationship alternation of taxas will not be analyzed if the sample size do not meet this requirement.
> 5. minimum_taxa_number. Minimum taxa detection samples in each group. The value should be an integer. If a taxa is only detected in a few samples less than the set value in either the control group or the treatment group, it will be filtered out. It is suggested the value is set larger than 6. 
> 6. minimum_taxa_prevalence. Minimum taxa median abundance in control group.The value should be in [0,1). If the median relative abundance of taxa is less than the set value, the taxa will be filtered out. This is another way to control the prevalence of taxa. One can specify this value to 0 if one do not want to filter out low abundance taxa.
> 7. control_label. Label of control group. Specify the control group label. For the demo project in this folder, it should be *H2029*.
> 8. treat_label. Label of treatment group. Specify the treatment group label. For the demo project in this folder, it should be *crc*.
> 9. group_label. Column name of group indicator. Specify the column name of group indicator. This should be consistent with the datafile. For the demo project in this folder, it should be *condition*.
> 10. fdr. Whether turn FDR adjustment on. Benjamini/Hochberg FDR p value adjustment are provided in the software. One can choose turning it on or off. It is suggested if the sample size in each group are larger than 50, the FDR adjust is on. Valid values include "ON" and "OFF".
> 11. ilr. Whether performing isometric log-ratio transformation on relative abundance data. Valid values include "ON" and "OFF".
> 12. confidence_level. The confidence interval for statisitical significantly altered relationship identification. The value should be in (0.7,1). 
> 13. AutoBandwidth. Whether turn auto optimization of bandwidth parameter on. The bandwidth is used in kernel density estimation. The larger bandwidth is, the estimated distribution is more smooth. On the other hand, a larger bandwidth will cause information missing. If the autobandwidth is off, the program will use 0.1 as default para. Turnning off this feature alao can save about one-third computing time.
> 14. kernel_str. Kernel used in the kernel density estimation. Valid kernels includes 'gaussian', 'tophat', 'epanechnikov', 'exponential', 'linear' and 'cosine'. The 'gaussian' and 'epanechnikov' are recommended.


### Arguments for multi-processing version PM2RA
The *PM_Analysis_multiprocessing.py* is the entrance of the multi-threading analysis. It takes following arguments.
> 1. fileplace. The place to store analysis result.
> 2. logs_place. The place to store analysis process logs. It is recommend using subfolder under fileplace. Please make sure the fileplace folder and logs_place folder do exist in your system. 
> 3. datafilename. The absolute path of the csv file to be analyzed. Both the datafile name and its relative or absolute path should be provided.
> 4. minimum_coexist_taxa_number. PM2RA analyzes the taxa relationship alternations between groups. Only samples with paired taxa are pipelined into the analysis. Take taxa1 and taxa2 as an example, only taxa1's abundance and taxa2's abundance are both larger than 0, the sample will be into the analysis. This parameter minimum_coexist_taxa_number specifies the minimum sample size of each group with detected taxa1 and taxa2. The relationship alternation of taxas will not be analyzed if the sample size do not meet this requirement.
> 5. minimum_taxa_number. Minimum taxa detection samples in each group. The value should be an integer. If a taxa is only detected in a few samples less than the set value in either the control group or the treatment group, it will be filtered out. It is suggested the value is set larger than 6. 
> 6. minimum_taxa_prevalence. Minimum taxa median abundance in control group.The value should be in [0,1). If the median relative abundance of taxa is less than the set value, the taxa will be filtered out. This is another way to control the prevalence of taxa. One can specify this value to 0 if one do not want to filter out low abundance taxa.
> 7. control_label. Label of control group. Specify the control group label. For the demo project in this folder, it should be *Healthy*.
> 8. treat_label. Label of treatment group. Specify the treatment group label. For the demo project in this folder, it should be *CD*.
> 9. group_label. Column name of group indicator. Specify the column name of group indicator. This should be consistent with the datafile. For the demo project in this folder, it should be *condition*.
> 10. processers. Specifiy the threads used in the analysis. The number should not exceed the number of logical processors.
> 11. fdr. Whether turn FDR adjustment on. Benjamini/Hochberg FDR p value adjustment are provided in the software. One can choose turning it on or off. It is suggested if the sample size in each group are larger than 50, the FDR adjust is on. Valid values include "ON" and "OFF".
> 12. ilr. Whether performing isometric log-ratio transformation on relative abundance data. Valid values include "ON" and "OFF".
> 13. confidence_level. The confidence interval for statisitical significantly altered relationship identification. The value should be in (0.7,1). 
> 14. AutoBandwidth. Whether turn auto optimization of bandwidth parameter on. The bandwidth is used in kernel density estimation. The larger bandwidth is, the estimated distribution is more smooth. On the other hand, a larger bandwidth will cause information missing. If the autobandwidth is off, the program will use 0.1 as default para. Turnning off this feature alao can save about one-third computing time.
> 15. kernel_str. Kernel used in the kernel density estimation. Valid kernels includes 'gaussian', 'tophat', 'epanechnikov', 'exponential', 'linear' and 'cosine'. The 'gaussian' and 'epanechnikov' are recommended.


### Outputs

PM2RA@XYL_Lab will output the analysis result in the specified folder. It includes

> 1. The *PM_scores_1D.csv*. This file contains a table listing the 1-dimensional taxa changes (in the column *pm*), the p value of PM(in the column *pvalue*), the q value(FDR adjusted pvalue) (in the column *qvalue*). If FDR adjust is off, q value will be not provided. 1-dimensional pm scores describe the taxa abundance changes between groups.
> 2. The *PM_scores.csv*. This file contains a table list the 2-dimenesional relationship alternations. The column *taxa1* and *taxa2* specify the corresponding 2 dimenesions. *pm1*, *pvalue1*, *qvalue1* and *pm2*, *pvalue2*, *qvalue2* are *taxa1*'s and *taxa2*'s 1D PM score respectively. The *PM_scores.csv* provides two types of 2-dimenesional PM scores. One is *raw_pm_2d* which brings the changes of taxa1, the changes of taxa2, and the covariance changes of taxa1 and taxa2. The other one is *co_PM_2D*. *co_PM_2D* removes out the effects of the changes of taxa1, the changes of taxa2. The *pvalue* and *qvalue* are corresponding statisitical testing results.


## Examples
A demodata (demodata_small.csv) can be found in this folder. To analysis the demo project, users can using following command.
### Single-processing PM2RA in linux platform
```python
source activate PM2RA # use pipenv shell in python env
python PM_Analysis.py demo_data/ demo_data/logs/ demo_data/demodata_small.csv 15 10 0 H2029 crc condition ON ON 0.95 ON gaussian 
```

### Multi-processing PM2RA in linux platform
```python
source activate PM2RA # use pipenv shell in python env
python PM_Analysis_multiprocessing.py demo_data/ demo_data/logs/ demo_data/demodata_small.csv 15 10 0 H2029 crc condition 6 ON ON 0.95 ON gaussian 
```