# multiGRNtools
### This Repository has Scripts for Multiomics based GRN inference methods

### For Tools installtion

## CellOracle installation using conda and pip

We recommend installing CellOracle in an independent conda environment to avoid dependent software conflicts. Please make a new python environment for celloracle and install dependent libraries in it.

```

conda create -n celloracle_env python=3.10

conda activate celloracle_env



```



You can install CellOracle and all dependencies with the following command.



```

pip install celloracle

```



### Check installation

You can check the installed library version as follows. Please make sure that all dependencies are appropriately installed.

In python console,



```

import celloracle as co

co.check_python_requirements()

```

## Pando Installation (https://github.com/quadbio/Pando)

```
devtools::install_github('quadbio/Pando')

```


## FigR Installation (https://buenrostrolab.github.io/FigR/)

```
devtools::install_github("buenrostrolab/FigR")

```

## DIRECT-NET Installation (https://github.com/zhanglhbioinfor/DIRECT-NET)
To make it easy to run DIRECT-NET in most common scRNA-seq and scATAC-seq data analysis pipelines, DIRECT-NET is now implemented within Seurat V4/Signac workflow. Please first install Seurat R pacakge (>= 4.0) via install.packages('Seurat').

DIRECT-NET R package can then be easily installed from Github using devtools:

```
devtools::install_github("zhanglhbioinfor/DIRECT-NET")

```
### Installation of other dependencies
### Install Signac pacakge : 

```
devtools::install_github("timoast/signac", ref = "develop")
```

### Install Cicero package: 
```
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
```
