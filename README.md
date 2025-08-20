# data_unibremen
Scripts for downloading, opening, read data files from the Universitat Bremen

### To build matching environment:

(1) Navigate to folder in terminal and run:
conda env create --file=environment.yml

(2) Optional: if running in Jupyter Notebook environment
After creating the environment, you will next need to add it to the list of available kernels. 
This does not happen automatically. In the command line, run:
python -m ipykernel install --user --name <NAME>
* Make sure to replace <NAME> with the name of the environment

If conda command not recognized by default with terminal launch:
source /opt/anaconda3/bin/activate

(3) Optional, to update environment file after modifying packages:
conda env export > environment.yml


# coordinates.py
functions to name and open files with coordinates data for different grids
https://data.seaice.uni-bremen.de/grid_coordinates/
*** currently only support remote file access ***

# SIC.py
functions to open AMSR-based SIC files at various resolutions
https://data.seaice.uni-bremen.de/amsr2/asi_daygrid_swath/
*** currently only support remote file access ***

# MultiYearIce.py
functions to open MYI concentration files
https://data.seaice.uni-bremen.de/MultiYearIce/ 
*** currently only support remote file access ***
