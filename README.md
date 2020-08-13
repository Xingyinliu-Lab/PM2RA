# PM2RA
Profile Monitoring for Microbial Community Alteration

Version 1.0.2

Updated date: 2019.06.17

Cite us:

Author and License<br>
===
Author: Zhi Liu<br>

Email: liuzhi@njum.edu.cn<br>

DEPENDENCIES
===
R 3.5.1 *(Other version not tested)*<br>
R packages:<br>
*gtools,MASS,sfsmisc,DMwR，pracma，qcc,stringr,ggplot2,igraph,RColorBrewer,pheatmap,reshape2,parallel.*<br>

Script tested on Mac and Linux system.<br>

*The average computation time of PM2RA for a dataset contains 100 features is 30 minutes (running on a Linux with 8 processors). 

USAGE<br>
-----
This PM2RA script is developed with two functional modules.<br>
1. 2D scanning(2DScan):  scanning for pairwise association change among the microbial community between two conditions, and the two-component community alterations form a community alteration(CA) network, where the nodes represent microbiota and the edges represent the measure of the community alteration (or association changes) of the two connected microbiota by PM score. Hub microbiotas in the CA network are the ones with extensively altered associations between two compared conditions.<br>

2. Multiple-dimensional search(MultiSearch):  where the PM score of any defined sub-community with two or more microbiotas could be calculated.<br>

2DScan usage<br>
-----

**Seven** parameters should be set for 2DScan function<br>
`Rscript PM2RA.R abundance.table.csv ConditionA-Label ConditionB-Label prevalence.filter project.name 2DScan Test-Type`<br>
1.  abundance table: the abundance table of microbiota in csv format, where the first colunm indicates the labels for samples from two conditions.<br><br>
**Input Example:**<br>
* Missing value is not allowed<br><br>
label	Bifidobacterium_animalis	Bifidobacterium_breve	Bacteroides_ovatus	Parabacteroides_distasonis<br>
H	1.64E-05	0.000212566	0.318188189	0.005378391<br>
H	0.000106283	0.001448253	0.142495246	0.011803262<br>
H	0.000205558	0.000363231	0.022839183	0.023493233<br>
H	5.72E-05	0.001193641	0.053880851	0.082181349<br>
H	0.000629523	0.001410879	0.027530822	0.030052417<br>
H	0.000157673	0.002483053	0.005934333	0.009705631<br>
H	0.000196215	0.000908662	0.016441175	0.018664944<br>
D	0.00030133	0.001388688	0.019018832	0.010586262<br>
D	0.000416957	0.000828074	0.045962177	0.043511827<br>
D	0.000153001	0.001358321	0.042005176	0.030979767<br>
D	0.00013665	0.003584426	0.032091651	0.039965943<br>
D	0.000140154	0.000798875	0.041534494	0.01880393<br>


2.  ConditionA-Label: label for samples from conditionA(e.g. H,disease,A)<br>
3.  ConditionB-Label: label for samples from conditionB(e.g. D,control,B)<br>
4.  prevalence filter: remove microbiotas with prevalence less than the prevalence filter(e.g. 0.1)<br>
5.  project.name: prefix attached to results files<br>
6.  function type: 2DScan<br>
7.  Statistic method: Kolmogorov-Smirnov test(ks,default) or Wilcoxon rank-sum test(wc)

**Run Example**<br>
`Rscript PM2RA.R Example.csv H D 0.1 2DscanExample 2DScan ks`<br>


MultiSearch usage<br>
-----

**Eight** parameters should be set for MultiSearch function. Except for those set for 2DScan, a microbiota list should be provided for MultiSearch.<br>

`Rscript PM2RA.R abundance.table.csv ConditionA-Label ConditionB-Label prevalence.filter project.name MultiSearch Test-type micriobita.list`<br>
1.  abundance table: the abundance table of microbiota in csv format, where the first colunm indicates the labels for samples from two conditions.<br><br>
**Input Example:**<br>
label	Bifidobacterium_animalis	Bifidobacterium_breve	Bacteroides_ovatus	Parabacteroides_distasonis<br>
H	1.64E-05	0.000212566	0.318188189	0.005378391<br>
H	0.000106283	0.001448253	0.142495246	0.011803262<br>
H	0.000205558	0.000363231	0.022839183	0.023493233<br>
H	5.72E-05	0.001193641	0.053880851	0.082181349<br>
H	0.000629523	0.001410879	0.027530822	0.030052417<br>
H	0.000157673	0.002483053	0.005934333	0.009705631<br>
H	0.000196215	0.000908662	0.016441175	0.018664944<br>
D	0.00030133	0.001388688	0.019018832	0.010586262<br>
D	0.000416957	0.000828074	0.045962177	0.043511827<br>
D	0.000153001	0.001358321	0.042005176	0.030979767<br>
D	0.00013665	0.003584426	0.032091651	0.039965943<br>
D	0.000140154	0.000798875	0.041534494	0.01880393<br>

2.  ConditionA-Label: label for samples from conditionA(e.g. H,disease,A)<br>
3.  ConditionB-Label: label for samples from conditionB(e.g. D,control,B)<br>
4.  prevalence filter: remove microbiotas with prevalence less than the prevalence filter(e.g. 0.1)<br>
5.  project.name: prefix attached to results files<br>
6.  function type: MultiSearch<br>
7.  Statistic method: Kolmogorov-Smirnov test(ks,default) or Wilcoxon rank-sum test(wc)
8.  microbita.list: a plain text file provides a list of micriobita，where each line contains one microbiota.<br>

**Run Example**<br>
`Rscript PM2RA.R Example.csv H D 0.1 MultiSearchExample MultiSearch ks list.txt`<br>


RELEASE NOTES
==
v1.0.2, 2019.06.17, adding Wilcoxon rank-sum test as an alternative test.<br>
v1.0.0, 2018.12.01

License
==
This project is covered under the Apache 2.0 License.


