# Tutorial Server Version Case2

## Input Data Download
### You can choose to download them directly from your browser by visiting the following URL: 
- [B_lymphoma Dataset (.zip) (288.0 MB)](https://zenodo.org/record/8160180/files/B_lymphoma.zip?download=1)
### You also can download them from cmd command：
    curl -o /path/to/save/location/B_lymphoma.zip https://zenodo.org/api/files/cf6453c0-0853-4633-a9d3-e571fb8ba47e/B_lymphoma.zip
### Dataset Structure:
#### **B_lymphoma_Dataset.zip**

<small>

- &nbsp;&nbsp;&nbsp;B_lymphoma_B.h5ad (66.4 MB)
- &nbsp;&nbsp;&nbsp;B_lymphoma_all.h5ad (182.6 MB)
- &nbsp;&nbsp;&nbsp;B_socre_df.csv (60.4 MB)
- &nbsp;&nbsp;&nbsp;Cell_names.tsv (268.8 kB)
- &nbsp;&nbsp;&nbsp;Gene_Cell.mtx (271.5 MB)
- &nbsp;&nbsp;&nbsp;Gene_Peak.mtx (8.9 MB)
- &nbsp;&nbsp;&nbsp;Gene_names.tsv (129.9 kB)
- &nbsp;&nbsp;&nbsp;Peak_Cell.mtx (859.7 MB)
- &nbsp;&nbsp;&nbsp;Peak_names.tsv (2.6 MB)

</small>

This project includes a tutorial that demonstrates the use and understanding of different parts of the project, using Server Version Case1 as an example. Before starting with this tutorial, make sure you have downloaded all the files of the project and placed them at the correct paths.

## Intermediate Result Data Download

All intermediate results files are in the `marsgt/Tutorial/Turtorial_server_version/Case2` directory. Before starting with this tutorial, you need to download these files.

After downloading the files, you need to alter the file paths in the code to match your own system environment. These file paths are typically in the file reading statements, such as in `open()` or `read()` functions. Ensure that you've correctly configured all file paths.

For example, if you see the following code:

```python
with open('path/to/your/file.txt', 'r') as file:
    data = file.read()
```
You need to change 'path/to/your/file.txt' to the actual path of your file.


