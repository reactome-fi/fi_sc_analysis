#!/bin/bash

# To run this, make sure the scanpy env is activated
# Remove all old files
rm -r dist_mac scpy4reactome_mac.pyz

# Install all dependencies, including scpy4reactome
pip install . --target=dist_mac

# Build the bundles. Make sure the main method should be entered like this.
shiv --site-packages dist_mac --compressed -p '/usr/bin/env python3' -o scpy4reactome_mac.pyz -e scpy4reactome.sc_json_server:main

# To run the package use. Make sure python version is 3.7.0.
# python scpy4reactome_mac.pyz
# Edit the version file and then copy the above pyz file to the server in their own OS folder and remove _mac in the file name.
