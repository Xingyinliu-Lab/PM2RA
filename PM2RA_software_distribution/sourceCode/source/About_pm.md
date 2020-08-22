<center>About PM2RA@XYL_Lab</center>
==========================
PM2RA@XYL_Lab: Detecting and Quantifying Relationship Alterations in Microbial Community

The human microbiome is a complex bacterial community that forms sub-communities based on shared niche specializations and specific interactions between individual microbes. The dysbiosis of gut microbiota relationship is associated with the pathogenesis of human disease. This software, *PM2RA@XYL_Lab* examines the relationship alteration (RA) in the microbiome between groups and can provide additional hints about the pathogenesis of human disease. For more details, please refer to:.


## <font color=#8A2BE2F>PM2RA Framework</font>
PM2RA is specifically designed to quantify the relationship alternation involving two or more microbes (“sub-community”) under different conditions. The basic idea of PM2RA analysis is to project the abundance data of two or more microbes under two conditions into the same space via Hoteling’s $T^{2}$ statistics, and compare the difference in the distribution of $T^{2}$ statistics to represent the RA between two conditions. We developed a new scoring scheme called PM (profile monitoring) score to quantify the RA of each sub-community under different conditions in five steps. The more the sub-community alters, the bigger the PM score is. Next, we build a RA network in which edges denote the corresponding PM score.

<center>

![fig1](md_source/fig1.svg)
</center>
<center>Fig 1 PM2RA framwork</center>  
<br />
The PM score of any defined sub-community that is referring to two or more microbes can be calculated with PM2RA. We proposed two main functions of PM2RA.  Every two microbes comprise a sub-community. After traversing all sub-communities, a weighted network is built to visualize the overall RAs. In the relationship alternation network G=(V, E), where V is the set of vertices representing microbes and E is the set of edges denoting the relationship alternation between the two conditions. The edge width and vertices size denotes PM score and topological degree, respectively. In this pairwise network, hub microbes, which have extensively altered associations between two compared conditions, can be identified. 


### <font color=#8A2BE2F>Input of PM2RA@XYL_Lab</font>

>1. One can feed either microbiome absoulte abundance table or relative abundance table to the software. The PM2RA uses relative abundance data as its input. If raw counts from sequencing is fed, it will convert the table to relative abundance automatically.
>2. PM2RA use kernel distribution estimation in the analysis processes. To get a higher trustable estimation, the sample size for each group is better to be larger than 10. 

### <font color=#8A2BE2F>Output of PM2RA@XYL_Lab</font>
PM2RA@XYL_Lab helps researchers to conduct analysis on the taxa relationship alternation between experimental groups. 

## <font color=#8A2BE2F>Getting started</font>

PM2RA@XYL_Lab is open source software. One can find source code of PM2RA@XYL_Lab at [github](https://github.com/Xingyinliu-Lab/PM2RA). Both GUI and CUI version of the software are provided. Users familiar with python or R lanugage can start with the CUI version. One can also try this graphical user interface version. 

The PM2RA@XYL_Lab is a cpu-intensive software. It takes about 10 minutes on a laptop computer with I7-8750H using 8 threads to conduct an PM 1D and 2D analysis for an approximately 60 taxa composed relative abundance table. If one have to analysis more taxas, it is recommended to run the scripts on a server PC. It is suggested one conduct a proper prevalence fiter strategy to remove taxa with low abundance or low detection rate before PM2RA analysis. PM2RA@XYL_Lab also provides the built-in prevalence filter.


## <font color=#8A2BE2F>Citing PM2RA@XYL_Lab</font>
If you use PM2RA@XYL_Lab for any published research, please include the following citation:
```
PM2RA: A Framework for Detecting and Quantifying Relationship Alterations in Microbial Community
```
