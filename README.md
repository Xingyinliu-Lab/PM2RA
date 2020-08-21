# <center><font color=#8A2BE2F>PM2RA@XYL_Lab</font></center>
PM2RA@XYL_Lab: Detecting and Quantifying Relationship Alterations in Microbial Community

The human microbiome is a complex bacterial community that forms sub-communities based on shared niche specializations and specific interactions between individual microbes. The dysbiosis of gut microbiota relationship is associated with the pathogenesis of human disease. This software, *PM2RA@XYL_Lab* examines the relationship alteration (RA) in the microbiome between health status  and can provide additional hints about the pathogenesis of human disease. For more details, please refer to:.

## <font color=#8A2BE2F>Input of PM2RA@XYL_Lab</font>

>1. One can feed either microbiome absoulte abundance table or relative abundance table to the software. The PM2RA uses relative abundance data as its input. If raw counts from sequencing is fed, it will convert the table to relative abundance automatically.
>2. PM2RA use kernel distribution estimation in the analysis processes. To get a higher trustable estimation, the sample size for each group is better to be larger than 10. 

## <font color=#8A2BE2F>About PM2RA@XYL_Lab</font>
PM2RA@XYL_Lab helps researchers to conduct analysis on the taxa relationship alternation between experimental groups. We provided both R and python implementation of this software. PM2RA@XYL_Lab is open source software. One can find source code of PM2RA@XYL_Lab at [github](https://github.com/Xingyinliu-Lab/PM2RA). Both GUI and CUI version of the software are provided. 

## <font color=#8A2BE2F>Folder Structure of this repository</font>
```
PM2RA # root of the repository
│
├─PM2RA_R_implementation # R implementation of PM2RA
│          
├─PM2RA_python_implementation # Python implementation of PM2RA
│  
├─PM2RA_software_distribution # GUI version of PM2RA@XYL_Lab
│  
└─README.md # this file
```
## <font color=#8A2BE2F>Getting started</font>

Both character user interface and graphical user interface version PM2RA@XYL_Lab are provided in this repo. Users familiar with python or R lanugage can start with the CUI version. This version provides more flexible functions. One can also try the graphical user interface version. 

The PM2RA@XYL_Lab is a cpu-intensive software. It takes about 10 minutes on a laptop computer with I7-8750H using 8 threads to conduct an PM 1D and 2D analysis for an approximately 60 taxa composed relative abundance table. If one have to analysis more taxas, it is recommended to run the scripts on a server PC. It is suggested one conduct a proper prevalence fiter strategy to remove taxa with low abundance or low detection rate before PM2RA analysis. PM2RA@XYL_Lab also provides the built-in prevalence filter.


## <font color=#8A2BE2F>Citing PM2RA@XYL_Lab</font>
If you use PM2RA@XYL_Lab for any published research, please include the following citation:
```
PM2RA: A Framework for Detecting and Quantifying Relationship Alterations in Microbial Community
```

