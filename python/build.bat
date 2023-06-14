@echo off
REM This batch script is modified directly from build.sh
REM To run this batch file, make sure python3.7, pip and shiv are installed correctly and set up as in PATH var

REM Remove all old files
ECHO Remove the old file and directory...
DEL dist_win scpy4reactome_win.pyz

REM Install all dependencies, including scpy4reactome
ECHO Install all dependencies...
pip install . --target=dist_win

REM Build the bundles. Make sure the main method should be entered like this.
ECHO Build the bundle...
REM This command is a little bit different from the .sh version: no -o  is set.
shiv --site-packages dist_win --compressed -o scpy4reactome_win.pyz -e scpy4reactome.sc_json_server:main

REM To run the package use. Make sure python version is 3.7.0.
REM python scpy4reactome_win.pyz
REM Edit the version file and then copy the above pyz file to the server in their own OS folder and remove _win in the file name.