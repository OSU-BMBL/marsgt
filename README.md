<h1 align="center">MarsGT: Multi-omics data analysis for rare population inference using single-cell graph transformer</h1>

## Description

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MarsGT, for rare cell identification from matched scRNA-seq (snRNA-seq) and scATAC-seq (snATAC-seq),includes genes, enhancers, and cells in a heterogeneous graph to simultaneously identify major cell clusters and rare cell clusters based on eRegulon.
<p align="center">
  <img src="./images/MarsGT%20Flowchart%201.jpg" alt="MarsGT Flowchart" width="900">
</p>

## Installation

### System Requirements

* Python Version == 3.8.0
* Hardware Architecture: x86_64
* Operating System: GNU/Linux(windows)

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
* You can directly install MarsGT using the pip command:
```bash
pip install --upgrade MarsGT
```

### Then, you can choose to use either the CPU or GPU version of PyTorch as per your preference
#### CPU Version
##### For linux system(torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.0.9+ torch_sparse-0.6.14):
    pip install https://download.pytorch.org/whl/cpu/torch-1.12.0%2Bcpu-cp38-cp38-linux_x86_64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_cluster-1.6.0%2Bpt112cpu-cp38-cp38-linux_x86_64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_scatter-2.0.9-cp38-cp38-linux_x86_64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_sparse-0.6.14-cp38-cp38-linux_x86_64.whl
##### For windows system(torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.0.9+ torch_sparse-0.6.14):
    pip install https://download.pytorch.org/whl/cpu/torch-1.12.0%2Bcpu-cp38-cp38-win_amd64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_scatter-2.0.9-cp38-cp38-win_amd64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_sparse-0.6.14-cp38-cp38-win_amd64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_cluster-1.6.0%2Bpt112cpu-cp38-cp38-win_amd64.whl
#### GPU Version
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If you are interested in accelerating your software's runtime using a GPU, I kindly recommend visiting the official PyTorch website at [https://pytorch.org/]. There, you will have the opportunity to select and download the CUDA-enabled version of PyTorch that best matches your system configuration. For example:
##### For linux system(torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.1.0+ torch_sparse-0.6.16)：
    pip install https://download.pytorch.org/whl/cu102/torch-1.12.0%2Bcu102-cp38-cp38-linux_x86_64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu102/torch_scatter-2.1.0%2Bpt112cu102-cp38-cp38-linux_x86_64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu102/torch_sparse-0.6.16%2Bpt112cu102-cp38-cp38-linux_x86_64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu102/torch_cluster-1.6.0%2Bpt112cu102-cp38-cp38-linux_x86_64.whl
##### For Windows system(torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.1.0+ torch_sparse-0.6.16)：
    pip install https://download.pytorch.org/whl/cu116/torch-1.12.0%2Bcu116-cp38-cp38-win_amd64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu116/torch_scatter-2.1.0%2Bpt112cu116-cp38-cp38-win_amd64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu116/torch_sparse-0.6.15%2Bpt112cu116-cp38-cp38-win_amd64.whl
    pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu116/torch_cluster-1.6.0%2Bpt112cu116-cp38-cp38-win_amd64.whl

## Data Download
### You can choose to download them directly from your browser by visiting the following URL: 
    https://zenodo.org/record/8160180
### You also can download them from cmd command：
    curl -o /path/to/save/location/B_lymphoma.zip https://zenodo.org/api/files/cf6453c0-0853-4633-a9d3-e571fb8ba47e/B_lymphoma.zip
    curl -o /path/to/save/location/Mouse_retina.zip https://zenodo.org/api/files/d749ff9e-ff3a-41a2-a922-c862cf962e66/Mouse_retina.zip
    curl -o /path/to/save/location/Tutorial_example.zip https://zenodo.org/api/files/d749ff9e-ff3a-41a2-a922-c862cf962e66/Tutorial_example.zip

