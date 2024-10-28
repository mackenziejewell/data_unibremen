# Description: Functions to generate url path to coordinate files for Uni Bremen data

from pyhdf.SD  import *
from metpy.units import units

def name_area_file(res: str, hem: str, include_url: bool = True):

    """Generate the name of the cell area file for a given resolution and hemisphere.
    INPUT:
    - res: str, resolution of desired file ('3125, '6250', '12500', '25000')
    - hem: str, hemisphere of desired file ('n', 's')
    - include_url: whether to include remote url in filename (True) or file name only (False)
    OUTPUT:
    - file: str, url path to coordinate file"""

    area_path = 'https://data.seaice.uni-bremen.de/grid_coordinates/GridCellArea/'

    # convert to km
    resstr = str(float(res)*units('m').to('km').magnitude)

    if include_url:
        file = area_path + f'PolStereo_GridCellArea_{hem}{resstr}km_Arctic.nc'
    else:
        file = f'PolStereo_GridCellArea_{hem}{resstr}km_Arctic.nc'
    
    return file
    
    
def name_coordinate_file(res: str, hem: str, include_url: bool = True):

    """Generate the name of the coordinate file for a given resolution and hemisphere.
    INPUT:
    - res: str, resolution of desired file ('3125, '6250', '12500', '25000')
    - hem: str, hemisphere of desired file ('n', 's')
    - include_url: whether to include remote url in filename (True) or file name only (False)
    OUTPUT:
    - file: str, url path to coordinate file"""

    coordinates_path = 'https://data.seaice.uni-bremen.de/grid_coordinates/'

    group = hem+res

    if hem == 'n':
        region = 'Arctic'
    elif hem == 's':
        region = 'Antarctic'

    # cell coordinate (lat/lon) file
    filepath = coordinates_path + f'{group}/'

    if res == '3125':
        coord_file = f'LongitudeLatitudeGrid-{group}-{region}3125.hdf'
    elif res == '25000':
        coord_file = f'PolStereo_LonLat-{hem}25km-{region}.hdf'
    else:
        coord_file = f'LongitudeLatitudeGrid-{group}-{region}.hdf'

    if include_url:
        file = filepath + coord_file
    else:
        file = coord_file

    return file


def open_coord_file(file):
    """Open coordinate file from locally downlaoded file."""
    # open file
    try:
        f = SD(file, mode=1) 
    # raise an error if file can't be opened 
    except Exception as e:
        print(f"{e}, error opening {file}")
    lon = f.select('Longitudes')[:]
    lat = f.select('Latitudes')[:]
    return lon, lat