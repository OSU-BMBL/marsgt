# Tutorial Server Version Case1

## Data Download
### You can choose to download them directly from your browser by visiting the following URL: 
- [Mouse_retina Dataset (.zip) (123.7 MB)](https://zenodo.org/record/8160180/files/Mouse_retina.zip?download=1)
### You also can download them from cmd command：
    curl -o /path/to/save/location/Mouse_retina.zip https://zenodo.org/api/files/d749ff9e-ff3a-41a2-a922-c862cf962e66/Mouse_retina.zip
### Dataset Structure:
#### **Mouse_retina_Dataset.zip**

<small>

- &nbsp;&nbsp;&nbsp;Cell_names.tsv (178.3 kB)
- &nbsp;&nbsp;&nbsp;Cell_type.tsv (37.4 kB)
- &nbsp;&nbsp;&nbsp;Gene_Cell.mtx (41.7 MB)
- &nbsp;&nbsp;&nbsp;Gene_Peak.mtx (2.1 MB)
- &nbsp;&nbsp;&nbsp;Gene_names.tsv (39.3 kB)
- &nbsp;&nbsp;&nbsp;Mouse_retina_BC.h5ad (39.7 MB)
- &nbsp;&nbsp;&nbsp;Mouse_retina_all.h5ad (64.9 MB)
- &nbsp;&nbsp;&nbsp;Peak_Cell.mtx (372.8 MB)
- &nbsp;&nbsp;&nbsp;Peak_names.tsv (1.4 MB)
- &nbsp;&nbsp;&nbsp;Total_socre_df.csv (46.6 MB)

</small>

This project includes a tutorial that demonstrates the use and understanding of different parts of the project, using Server Version Case1 as an example. Before starting with this tutorial, make sure you have downloaded all the files of the project and placed them at the correct paths.

## Download and Configure Files

All intermediate results files are in the `marsgt/Tutorial/Turtorial_server_version/Case1` directory. Before starting with this tutorial, you need to download these files.

After downloading the files, you need to alter the file paths in the code to match your own system environment. These file paths are typically in the file reading statements, such as in `open()` or `read()` functions. Ensure that you've correctly configured all file paths.

For example, if you see the following code:

```python
with open('path/to/your/file.txt', 'r') as file:
    data = file.read()
You need to change 'path/to/your/file.txt' to the actual path of your file.

