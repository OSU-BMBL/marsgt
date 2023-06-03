<h1 align="center">MarsGT: Multi-omics data analysis for rare population inference using single-cell graph transformer</h1>

## Description

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;DeepMARS, for rare cell identification from matched scRNA-seq (snRNA-seq) and scATAC-seq (snATAC-seq),includes genes, enhancers, and cells in a heterogeneous graph to simultaneously identify major cell clusters and rare cell clusters based on eRegulon.
<p align="center">
  <img src="./images/MarsGT%20Flowchart%201.jpg" alt="MarsGT Flowchart" width="900">
</p>

## Installation

### System Requirements

* Python Version ≥ 3.8.0
* Hardware Architecture: x86_64
* Operating System: GNU/Linux

### Dependencies, MarsGT has the following dependencies:
* anndata==0.8.0
* dill==0.3.4
* matplotlib==3.5.1
* numpy==1.22.3
* pandas==1.4.2
* scipy==1.9.1
* seaborn==0.11.2
* scikit-learn==1.1.2
* torch==1.12.0
* torch-geometric==2.1.0.post1
* torchmetrics==0.9.3
* xlwt==1.3.0
* tqdm==4.64.0
* scanpy==1.9.1
* leidenalg==0.8.10
* ipywidgets==8.0.6

### Installation Steps
* Install the required dependencies using pip:
```bash
pip install -r requirements.txt
```
* use pip to install MarsGT:
```bash
pip install MarsGT
```
## GPU Acceleration

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If you are interested in accelerating your software's runtime using a GPU, I kindly recommend visiting the official PyTorch website at [https://pytorch.org/get-started/previous-versions/](https://pytorch.org/get-started/previous-versions/). There, you will have the opportunity to select and download the CUDA-enabled version of PyTorch that best matches your system configuration.
