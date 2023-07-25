<h1 align="center">MarsGT: Multi-omics data analysis for rare population inference using single-cell graph transformer</h1>

## Description

MarsGT, for rare cell identification from matched scRNA-seq (snRNA-seq) and scATAC-seq (snATAC-seq),includes genes, enhancers, and cells in a heterogeneous graph to simultaneously identify major cell clusters and rare cell clusters based on eRegulon.
<p align="center">
  <img src="./images/MarsGT%20Flowchart%201.jpg" alt="MarsGT Flowchart" width="900">
</p>

## Installation

### System Requirements

* Python Version >=3.8.0
* Hardware Architecture: x86_64
* Operating System: GNU/Linux or Windows or MacOS

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

The installation process involves some optional and necessary steps. Here's the detailed breakdown:

1. **Recommended Step:** Create a new environment, you should use python 3.8.0.

    ```bash
    conda create --name marsgt python=3.8
    conda activate marsgt
    ```

2. **Necessary Step:** You need to install either the CPU or GPU version of PyTorch as per your preference, We recommend using the GPU version, which has a faster running speed compared to the CPU version:

    - **CPU Version**
        - For Linux system (torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.0.9+ torch_sparse-0.6.14):
        
            ```bash
            pip install https://download.pytorch.org/whl/cpu/torch-1.12.0%2Bcpu-cp38-cp38-linux_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_cluster-1.6.0%2Bpt112cpu-cp38-cp38-linux_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_scatter-2.0.9-cp38-cp38-linux_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_sparse-0.6.14-cp38-cp38-linux_x86_64.whl
            ```

        - For Windows system (torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.0.9+ torch_sparse-0.6.14):

            ```bash
            pip install https://download.pytorch.org/whl/cpu/torch-1.12.0%2Bcpu-cp38-cp38-win_amd64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_scatter-2.0.9-cp38-cp38-win_amd64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_sparse-0.6.14-cp38-cp38-win_amd64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_cluster-1.6.0%2Bpt112cpu-cp38-cp38-win_amd64.whl
            ```
       - For MacOS system (torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.0.9+ torch_sparse-0.6.14):

            ```bash
            conda install pytorch==1.12.0
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_scatter-2.0.9-cp38-cp38-macosx_10_15_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_cluster-1.6.0-cp38-cp38-macosx_10_15_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcpu/torch_sparse-0.6.14-cp38-cp38-macosx_10_15_x86_64.whl
            ```

    - **GPU Version**
        - Please visit the official PyTorch website at [PyTorch](https://pytorch.org/) to select and download the CUDA-enabled version of PyTorch that best matches your system configuration.
        - For linux system(You need to select the version that is compatible with your system's graphics card. For example: torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.1.0+ torch_sparse-0.6.16):
          
             ```bash
            pip install https://download.pytorch.org/whl/cu102/torch-1.12.0%2Bcu102-cp38-cp38-linux_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu102/torch_scatter-2.1.0%2Bpt112cu102-cp38-cp38-linux_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu102/torch_sparse-0.6.16%2Bpt112cu102-cp38-cp38-linux_x86_64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu102/torch_cluster-1.6.0%2Bpt112cu102-cp38-cp38-linux_x86_64.whl
             ```
        - For Windows system(You need to select the version that is compatible with your system's graphics card. For example: torch-1.12.0+ torch_cluster-1.6.0+ torch_scatter-2.1.0+ torch_sparse-0.6.16):

             ```bash
            pip install https://download.pytorch.org/whl/cu116/torch-1.12.0%2Bcu116-cp38-cp38-win_amd64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu116/torch_scatter-2.1.0%2Bpt112cu116-cp38-cp38-win_amd64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu116/torch_sparse-0.6.15%2Bpt112cu116-cp38-cp38-win_amd64.whl
            pip install https://data.pyg.org/whl/torch-1.12.0%2Bcu116/torch_cluster-1.6.0%2Bpt112cu116-cp38-cp38-win_amd64.whl
            ```
             
        - For MacOS system(According to the official PyTorch documentation, CUDA is not available on MacOS, please use the default package):

3. **Necessary Step:** You can directly install MarsGT using the pip command:

    ```bash
    pip install --upgrade MarsGT
    ```
             
## Data Download
### You can choose to download them directly from your browser by visiting the following URL: 
- [Tutorial_example (.zip) (5.3 MB)](https://zenodo.org/record/8163160/files/Tutorial_example.zip?download=1)
- [Case1: Mouse_retina Dataset (.zip) (123.7 MB)](https://zenodo.org/record/8163160/files/Mouse_retina.zip?download=1)
- [Case2: B_lymphoma Dataset (.zip) (288.0 MB)](https://zenodo.org/record/8163160/files/B_lymphoma.zip?download=1)
- [Case3: PBMCs_Dataset (.zip) (2.5 GB)](https://zenodo.org/record/8163160/files/PBMCs.zip?download=1)

### You also can download them from cmd command：
    curl -o Tutorial_example.zip https://zenodo.org/api/files/7ca78984-0e31-48cf-8b48-9544099d57bb/Tutorial_example.zip
    curl -o Mouse_retina.zip https://zenodo.org/api/files/7ca78984-0e31-48cf-8b48-9544099d57bb/Mouse_retina.zip
    curl -o B_lymphoma.zip https://zenodo.org/api/files/7ca78984-0e31-48cf-8b48-9544099d57bb/B_lymphoma.zip
    curl -o PBMCs.zip https://zenodo.org/api/files/7ca78984-0e31-48cf-8b48-9544099d57bb/PBMCs.zip

## MarsGT Model Tutorial

We have curated tutorials to assist you in operating the MarsGT model. You can locate these tutorials in the marsgt/Tutorial directory of the project. Additionally, we have provided an Example Dataset to aid you in testing and acquainting yourself with the MarsGT functionality:  
**Note:The execution of this code approximately takes 2 hours and requires about 90GB of memory.**     
&nbsp;&nbsp;&nbsp;&nbsp;[**Example Dataset** ↗](https://github.com/mtduan/marsgt/blob/main/Tutorial/Tutorial_for_example_data.ipynb)   
 

Beyond the aforementioned resources, we offer two versions of the tutorial designed to reproduce the results of the paper:

1. [**Local Version** ↗](https://github.com/mtduan/marsgt/tree/main/Tutorial/Tutorial_local_version): This version is specifically created for running the model locally on your computer. Utilize this version if you are aiming to reproduce the results documented in the paper.

2. [**Server Version** ↗](https://github.com/mtduan/marsgt/tree/main/Tutorial/Turtorial_server_version): This version is crafted for running the model on a server. Opt for this version if your goal is to reproduce both the model running process and the subsequent results.

To begin the tutorial, select the version that suits your needs or choose the example dataset. Subsequently, follow the directives provided in the README file within the respective directory or refer to the Jupyter notebook for the example dataset.
