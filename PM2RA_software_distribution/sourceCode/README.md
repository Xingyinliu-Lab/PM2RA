# Graphic User Interface PM2RA@XYL_Lab

This is the source code for Graphic User Interface version (Python implemented) of PM2RA@XYL_Lab.

## Installation recommendations
It is recommended to install GUI-PM2RA@XYL_Lab as follows
### using conda
```python
conda create -n PM2RA python=3.6
source activate PM2RA
conda install numpy==1.19.1
conda install pandas==1.1.0
conda install scipy==1.5.2
conda install scikit-learn==0.23.2
conda install networkx==2.4
conda install matplotlib==3.1.1
conda install PyQt5==5.15.0
```
### using pipenv
```python
pipenv install --python 3.6.7
pipenv shell
pip install numpy==1.19.1
pip install pandas==1.1.0
pip install scipy==1.5.2
pip install scikit-learn==0.23.2
pip install networkx==2.4
pip install matplotlib==3.1.1
pip install PyQt5==5.15.0
```
Note: only python 3.6.7 was tested.


## Getting started

```python
pipenv shell#source activate PM2RA
python PM.py
```
