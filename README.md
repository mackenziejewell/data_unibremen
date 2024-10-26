# data_UniBremen
Scripts for downloading, opening, read data files from the Universitat Bremen

# To build matching environment:
Navigate to folder in terminal and run:
conda env create --file=environment.yml

To update environment file after modifying packages:
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