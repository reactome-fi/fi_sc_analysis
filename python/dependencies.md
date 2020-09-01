## The following python packages are required:

- python version: 3.7; conda install -c anaconda python=3.7

- jsonrpclib: ~/miniconda3/envs/scanpy/bin/pip install jsonrpclib-pelix, to support json rpc srever. This is required for python 3. The original package jsonrpclib (https://pypi.org/project/jsonrpclib/) cannot work anymore.

- scanpy: https://anaconda.org/bioconda/scanpy, version 1.5.1

- scVelo: installed via pip, version 0.2.2. This is used for RNA velocity analysis and visualization. For details, see https://scvelo.readthedocs.io/installation.html.  

- magic: https://github.com/KrishnaswamyLab/MAGIC#installation-with-pip for imputation (note: dca cannot work for some reason)

- leidenalg: conda install -c conda-forge leidenalg

- ipython: conda install -c conda-forge ipython, version 7.16.1, for better python development

- Scribe-py: Used to calculate restricted directed mutual information to infer causality between TFs and targets. Install by /Users/wug/miniconda3/envs/scanpy/bin/pip install git+https://github.com/aristoteleo/Scribe-py.git (August 26, 2020 snapshot).