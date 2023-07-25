# Tutorial Server Version Case3

## Input Data Download
### You can choose to download them directly from your browser by visiting the following URL: 
- [PBMCs Dataset (.zip) (2.5 GB)](https://zenodo.org/record/8163160/files/PBMCs.zip?download=1)
### You also can download them from cmd command：
    curl -o PBMCs.zip https://zenodo.org/api/files/7ca78984-0e31-48cf-8b48-9544099d57bb/PBMCs.zip
### Dataset Structure:
#### **PBMCS_Dataset.zip**
<small>

- Cell_emb10.csv (80.6 MB)
- Cell_names.tsv (2.0 MB)
- Gene_Cell.mtx (2.5 GB)
- Gene_names.tsv (107.7 kB)
- Gene_Peak.mtx (158.7 MB)
- Peak_Cell.mtx (3.6 GB)
- Peak_names.tsv (303.4 kB)
- ycpu.h5ad (1.6 GB)

- Sample
    - HD1 
    - HD2 
    - P1 ...    
- Dataframe
    - 1_total_egrn.csv (15.4 MB)
    - 9_total_egrn.csv (14.5 MB)
    - 12_total_egrn.csv (3.4 MB) ...  

</small>

This project includes a tutorial that demonstrates the use and understanding of different parts of the project, using Server Version Case3 as an example. Before starting with this tutorial, make sure you have downloaded all the files of the project and placed them at the correct paths.

## Intermediate Result Data Download

All intermediate results files are in the `marsgt/Tutorial/Turtorial_server_version/Case3` directory. Before starting with this tutorial, you need to download these files.

After downloading the files, you need to alter the file paths in the code to match your own system environment. These file paths are typically in the file reading statements, such as in `open()` or `read()` functions. Ensure that you've correctly configured all file paths.

For example, if you see the following code:

```python
with open('path/to/your/file.txt', 'r') as file:
    data = file.read()
```
You need to change 'path/to/your/file.txt' to the actual path of your file.




