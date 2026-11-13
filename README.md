<h1 align="center">MarsGT: Multi-omics data analysis for rare population inference using single-cell graph transformer</h1>

> **Fork note (branch `modern-deps`).** This fork runs on current Python/PyTorch (tested: Python 3.13, torch 2.11 + CUDA 13, torch_geometric 2.8, pandas 2.2, anndata 0.13) and fixes the stage-2 NaN embeddings reported in OSU-BMBL/marsgt#1, #3, #4. Install with:
>
> ```bash
> pip install "git+https://github.com/guarelin/marsgt.git@modern-deps"
> ```
>
> Changes (one commit each): NumPy ≥1.23 label indexing; AnnData `dtype=` removal; `torch.sparse_coo_tensor`; `ReduceLROnPlateau(verbose=)` removal; `DataFrame.append` → `pd.concat` (previously made eGRN silently return empty); stage-2 KL computed with `log_target=True` (the NaN fix). `pyproject.toml` declares version ranges; the pinned `requirements.txt` below is the original upstream environment. The PyG extensions (`torch-sparse`, `torch-scatter`) are not imported by MarsGT and are not needed.
>
> **Sub-sampling modes** (`batch_select_whole(..., sampling=...)`): `'upstream'` (default) is the original code, unchanged. Its per-feature probabilities are misaligned with the features, because `list(set(...))` reorders them after the probabilities are computed, so the weights are effectively shuffled. `'aligned'` is the same weighting with that bug fixed. `'paper'` implements the specificity-based scheme of the MarsGT paper (Wang et al. 2024, *Nat Commun* 15:338, doi:10.1038/s41467-023-44570-8, eqs. 2-5): keep genes above the cell's first-quartile count, then sample with probability proportional to the share of each feature's total counts found in the cell. See `tests/test_sampling.py`.
>
> **Batch features built on the GPU** (0.2.1.post3): each batch's dense gene and peak inputs (rows of `RNA_matrix`/`ATAC_matrix` across *all* cells) are built on the GPU from one resident sparse copy of each matrix, so the original per-batch CPU densification and copy are gone. Values, memory layout and therefore training are identical (`tests/test_device_rows.py`). The resident copy costs about 1-3 GB of GPU memory at 180k cells. Turn it off with `MarsGT.utils.GPU_FEATURES = False`, and free it with `MarsGT.utils.clear_feature_cache()`.
>
> **Prediction with dropout off, and per-epoch loss history** (0.2.1.post4):
> - The released `MarsGT_pred` predicted with dropout (0.3) still active, so repeated predictions of the same trained model disagreed: ARI about 0.3, with about 40% of cells changing cluster. It now predicts in eval mode by default; `eval_mode=False` restores the original behavior.
> - Both training stages record mean losses per epoch in `.history`.
> - Training one epoch at a time on the same objects is bitwise-identical to one multi-epoch call, so wrappers can add early stopping without changing the math (`tests/test_training_controls.py`).
>
> **Paper-faithful options** (0.2.1.post5; all off by default, so the released behavior is unchanged):
> - `NodeDimensionReduction(..., num_relations=4)`: all four relations pass messages, including peak↔cell, as the paper describes. The released `num_relations=2` silently drops relations 2 and 3, so ATAC values have no effect on stage 1.
> - `MarsGT_pred(..., pred_graph='training')`: prediction uses the same four-relation, both-direction graph as training. The released prediction graph labels peak→cell edges as relation 1, which in training is cell→gene.
> - `MarsGT(..., cosine_loss=True)`: adds the paper's stage-2 cosine term (eq. 12), computed on the current batch. The released version used a stale stage-1 tensor and was commented out.
>
> - `cosine_reduce='mean'` and `cosine_self_pairs=False` (0.2.1.post6; on both `NodeDimensionReduction` and `MarsGT`): average the within-cluster cosine term over clusters, rather than summing it, and drop each cell's similarity with itself. That bounds the term to [−1, 1]. The released sum grows with the number of clusters per batch and outweighed the other losses in our runs. The defaults reproduce the released code exactly.
>
> See `tests/test_paper_options.py`.

## Description

MarsGT: Multi-omics analysis for rare population inference using single-cell Graph Transformer, which is an end-to-end deep learning model for rare cell population identification from scMulti-omics data.

<p align="center">
  <img src="./images/MarsGT%20Flowchart%203.svg" alt="MarsGT Flowchart" width="900">
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

1. **Recommended Step:** Create a new environment, you should use python 3.8.

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
### You can choose to download them directly from your browser by the following URL: 
- [Tutorial_example (.zip) (5.3 MB)](https://zenodo.org/record/8163160/files/Tutorial_example.zip?download=1)
- [Case1: Mouse_retina Dataset (.zip) (123.7 MB)](https://zenodo.org/record/8163160/files/Mouse_retina.zip?download=1)
- [Case2: B_lymphoma Dataset (.zip) (288.0 MB)](https://zenodo.org/record/8163160/files/B_lymphoma.zip?download=1)
- [Case3: PBMCs_Dataset (.zip) (2.5 GB)](https://zenodo.org/record/8163160/files/PBMCs.zip?download=1)

### You also can download them from cmd command：
    curl -L -o Tutorial_example.zip https://zenodo.org/record/8163160/files/Tutorial_example.zip?download=1
    curl -L -o Mouse_retina.zip https://zenodo.org/record/8163160/files/Mouse_retina.zip?download=1
    curl -L -o B_lymphoma.zip https://zenodo.org/record/8163160/files/B_lymphoma.zip?download=1
    curl -L -o PBMCs.zip https://zenodo.org/record/8163160/files/PBMCs.zip?download=1

### TF Binding Sites Data

We retrieved the genome browser track file from JASPAR, which stores all known TF binding sites for each TF. A p-value score was provided in JASPAR. Download Links：
- **Human (hg38 genome)**: [hg38_lisa_500.qsave](https://zenodo.org/records/14207730/files/hg38_lisa_500.qsave?download=1)
- **Mouse (mm10 genome)**: [mm10.qsave](https://zenodo.org/records/14207730/files/mm10.qsave?download=1)

## MarsGT Model Tutorial

We have curated tutorials to assist you in operating the MarsGT model. You can locate these tutorials in the marsgt/Tutorial directory of the project. Additionally, we have provided an Example Dataset to aid you in testing and acquainting yourself with the MarsGT functionality: 

&nbsp;&nbsp;&nbsp;&nbsp;[**Example Dataset** ↗](https://github.com/mtduan/marsgt/blob/main/Tutorial/Tutorial_for_example_data.ipynb)  

### Note-1
The execution of this code approximately takes 2 hours and requires about 250GB of memory.

### Note-2
More epochs will lead to better performance, but will also require a longer duration. Please balance this according to your needs to set parameters 'epochs' and 'num_epochs'.

### Note-3
For the usage of MAESTRO in this context, please refer to the process in `Tutorial/Tutorial_server_version/Case1/MAESTRO_Gene_Peak_Calculation.R`.

Beyond the aforementioned resources, we offer two versions of the tutorial designed to reproduce the results of the paper:

1. [**Local Version** ↗](https://github.com/mtduan/marsgt/tree/main/Tutorial/Tutorial_local_version): This version is specifically created for running the model locally on your computer. Utilize this version if you are aiming to reproduce the results documented in the paper.

2. [**Server Version** ↗](https://github.com/mtduan/marsgt/tree/main/Tutorial/Turtorial_server_version): This version is crafted for running the model on a server. Opt for this version if your goal is to reproduce both the model running process and the subsequent results.

To begin the tutorial, select the version that suits your needs or choose the example dataset. Subsequently, follow the directives provided in the README file within the respective directory or refer to the Jupyter notebook for the example dataset.
