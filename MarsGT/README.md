# MarsGT Project

## Introduction

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Welcome to the MarsGT project. This repository contains all the source code for the MarsGT project, which is a sophisticated and groundbreaking initiative aimed at rare cell identification from matched scRNA-seq (snRNA-seq) and scATAC-seq (snATAC-seq), includes genes, enhancers, and cells in a heterogeneous graph to simultaneously identify major cell clusters and rare cell clusters based on eRegulon.

## Source Code

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The source code located in this repository is organized into various files, each serving a specific purpose in the overall functioning of the MarsGT project.Here is an overview of the main files and their content:

### MarsGT  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── __init__.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── conv.py : Perform convolutional operations on graph data to extract valuable information.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── egrn.py : Compute and extract the EGRN between RNA sequencing data and ATAC sequencing data.   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── marsgt_model.py : Node dimensionality reduction as well as predicting cell type or gene-peak relationships.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── utils.py : Perform data preprocessing, initialize clusters, calculate performance metrics, and more.  
