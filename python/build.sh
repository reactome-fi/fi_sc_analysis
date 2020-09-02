#!/bin/bash

# To run this, make sure the scanpy env is activated
# Remove all old filess
rm -r dist scpy4reactome.pyz

# Install all dependencies, including scpy4reactome
/Users/wug/miniconda3/envs/scanpy/bin/pip install . --target=dist

# Build the bundles
shiv --site-packages dist --compressed -p '/Users/wug/miniconda3/envs/scanpy python' -o scpy4reactome.pyz -e scpy4reactome.ScJSONServer

# To run the package use. Make sure python version is 3.7.0.
# python scpy4reactome.pyz