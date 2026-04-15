# multiGRNtools
##This Repository has Scripts for Multiomics based GRN inference methods

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
