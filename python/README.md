### Note

This Python package is used to help **building regulatory networks** composed of both transcription factors and pathways in Reactome. The actual Python code used to build regulatory networks is hosted at [https://github.com/reactome-fi/sc_regulatory_network](https://github.com/reactome-fi/sc_regulatory_network "sc_regulatory_network"). 

### Installation 

- PyPI install using pip from the source directly to get the latest version. **It is strongly recommended to create a conda environment for python-based single cell data analysis using this or other packages.**

   ``` 
   git clone https://github.com/reactome-fi/fi_sc_analysis
   cd python
   pip install .
   ```

**If you have a difficult to install it, following the following steps after creating an conda env:**

1. pip install scanpy
2. pip install pyscenic
3. pip install scvelo
4. pip install .

**Note: The setup.py has been updated for the above procedures**

- You may install our released version by simply running the following:

```
pip install scpy4reactome
```

### Build the standalone Python applications

- To build the standalone Python application using shiv under Mac, use build.sh. Some configuration in the build.sh may 
  be adjusted for the actual local system.

- To build the standalone Python application using shiv for Windows, use build.bat. Some configuration in this file may 
  be adjusted based on the actual local system.

### Notes to be deleted:

- For fast cluster, install fastcluster via pip: pip install fastcluster (https://pypi.org/project/fastcluster/). 
  This is used for hierarchical clustering view.
