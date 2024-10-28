from datetime import datetime, timedelta
import xarray as xr
import numpy as np
import cartopy
import cartopy.crs as ccrs
from metpy.units import units
import urllib
import os

from .coordinates import open_coord_file, name_area_file, name_coordinate_file



def extract_variables(ds):
    """Extract variables from dataset."""
    data = {}
    data['proj'] = grab_projection(ds)
    data['ds'] = ds

    data['x'] = ds.x.values * units(ds.x.units)
    data['y'] = ds.y.values * units(ds.y.units)
    data['xx'], data['yy'], = np.meshgrid(data['x'], data['y'])
    data['sic'] = ds.z.values * units('%')
    
    return data

def grab_projection(ds, quiet = True):

    """Grab projection info from Uni B AMSR2-AMSRE sea ice concentration data (doi: 10.1029/2005JC003384)
    https://seaice.uni-bremen.de/sea-ice-concentration/amsre-amsr2/

INPUT: 
- ds: data opened with xarray
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- projection: cartopy projection from data projection info

Latest recorded update:
10-25-2024
    """
    # grab projection info
    #---------------------
    CRS = ds.polar_stereographic.attrs

    # grab parameters from crs spatial attributes
    central_meridian = int(CRS['straight_vertical_longitude_from_pole'])
    semimajor = CRS['semi_major_axis']
    inv_flat = CRS['inverse_flattening']
    standard_parallel = int(CRS['standard_parallel'])
    
    if quiet != True:
        print(f'>>> data provided in polar_stereographic projection')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - inverse_flattening: {inv_flat}')
        print(f'  - straight_vertical_longitude_from_pole: {central_meridian}')
        print(f'  - standard_parallel: {standard_parallel}')
        print(f'  - proj4text: {ds.attrs["proj4string"]}')
        
    # create projection from info
    projection = ccrs.NorthPolarStereo(central_longitude=central_meridian, 
                                           globe=ccrs.Globe(semimajor_axis = semimajor,inverse_flattening=inv_flat),
                                           true_scale_latitude=standard_parallel)

    return projection


def name_SIC_file(date, res: str, hem: str, include_url: bool = True):

    """Construct filename for Uni Bremen AMSR2-AMSRE sea ice concentration data (doi: 10.1029/2005JC003384)
    INPUT: 
    - date: datetime object for desired file
    - res: str, resolution of desired file ('3125, '6250', '12500', '25000')
    - hem: str, hemisphere of desired file ('n', 's')
    - include_url: whether to include remote url in filename (True) or file name only (False)
    """
    group = hem+res

    dataset_path = 'https://data.seaice.uni-bremen.de/'

    # create strftime for filename based on date
    if date >= datetime(2012, 7, 2):
        file_strftime = f'amsr2/asi_daygrid_swath/{group}/netcdf/%Y/asi-AMSR2-{group}-%Y%m%d-v5.4.nc'
        
    elif (date >= datetime(2002, 6, 1)) & (date <= datetime(2011, 10, 4)):  
        file_strftime = f'amsre/asi_daygrid_swath/{group}/netcdf/%Y/asi-n6250-%Y%m%d-v5.4.nc'

    # construct filename
    if include_url:
        file = dataset_path + date.strftime(file_strftime)
    else:
        file = date.strftime(file_strftime).split('/')[-1]

    return file


def open_remote_file(date, res : str='6250', hem : str='n', 
                     method : str = 'xarray', 
                     coordinates: bool = False,
                     area: bool = False,
                     include_units = False, quiet = True):
    
    """Use xarray to open remote files from HTTPS server:
    https://data.seaice.uni-bremen.de/

INPUT: 
- date: datetime object for desired file
- res: str, resolution of desired file ('6250' or '3125')
- hem: str, hemisphere of desired file ('n' or 's')
- method: str, method to attempt to laod data 
    ('xarray' for reading into memory only, sometimes has trouble
    'urllib' to temporarily download file to disk and then read with xarray)
- coordinates: bool, whether or not to download lat/lon coordinate data
- area: bool, whether or not to download cell area data
- include_units: bool, whether or not to return data with units
- quiet: bool, whether or not to supress print statements

OUTPUT:
- data: dictionary with nc data

Latest recorded update:
10-25-2024
    """
    # construct filename
    file = name_SIC_file(date, res, hem)

    if not quiet:
        print(f'>>> opening {file}')

    # use xarray to open remote file
    if method == 'xarray':

        # for some reason need to add this for now
        # see https://github.com/Unidata/netcdf4-python/issues/1043
        extrabit = '#mode=bytes'
        ds = xr.open_dataset(file+extrabit, decode_cf=False)
        
        # Extract variables from dataset
        data = extract_variables(ds)

    elif method == 'urllib':
        urllib.request.urlretrieve(file, "temp.nc")
        ds = xr.open_dataset("temp.nc")
        
        # Extract variables from dataset
        data = extract_variables(ds)

        # remove temporary file
        os.remove("temp.nc")

    # download cell area file
    #------------------------------
    if area:

        # cell area file
        area_file = name_area_file(res, hem)

        # download area file
        urllib.request.urlretrieve(area_file, 'tmp_area.nc')

        # open area file
        ds_area = xr.open_dataset('tmp_area.nc')
        data['area'] = ds_area['data'].values * units('km^2')
        os.remove("tmp_area.nc")

    # download cell coordinate file
    #------------------------------
    if coordinates:

        # cell coordinate (lat/lon) file
        coord_file = name_coordinate_file(res, hem)

        # download coordinate file
        urllib.request.urlretrieve(coord_file, 'tmp_coord.hdf')

        # extract lon/lat from coordinate file
        lon, lat = open_coord_file('tmp_coord.hdf')
        data['lon'] = lon * units('degreeE')
        data['lat'] = lat * units('degreeN')
        os.remove("tmp_coord.hdf")



    # remove units if desired
    if not include_units:
        for key in data.keys():
            if key not in ['proj', 'ds']:
                data[key] = data[key].magnitude



    return data




def open_local_file(date, res : str='6250', hem : str='n', 
                    main_path : str = '/Volumes/Seagate_Jewell/KenzieStuff/',
                     coordinates: bool = False,
                     area: bool = False,
                     include_units = False, quiet = True,):
    
    """Use xarray to open locally stored SIC files from UniBremen (https://data.seaice.uni-bremen.de/)

    Assumes following directory structure:

    main_path/
    ├── UniB-ASI-SIC-{res}{hem}/
    |   ├── LongitudeLatitudeGrid-{hem}{res}-Arctic.nc
    |   ├── PolStereo_GridCellArea_{hem}{res}_Arctic.nc
    |   |
    │   ├── asi-AMSR2-{hem}{res}-{date.year}/
    │       ├── asi-AMSR2-{hem}{res}-{date}-v5.4.nc


INPUT: 
- date: datetime object for desired file
- res: str, resolution of desired file ('6250' or '3125')
- hem: str, hemisphere of desired file ('n' or 's')
- main_path: str, directory where files are locally stored 
    - looks for subfolder named UniB-ASI-SIC-{}{} where {} will be replaced with res and hem
- coordinates: bool, whether or not to download lat/lon coordinate data
- area: bool, whether or not to download cell area data
- include_units: bool, whether or not to return data with units
- quiet: bool, whether or not to supress print statements

OUTPUT:
- data: dictionary with nc data

Latest recorded update:
10-25-2024
    """
    # construct filename
    filename = name_SIC_file(date, res, hem, include_url = False)

    # construct path
    path = main_path + f'UniB-ASI-SIC-{hem}{res}/'

    # determine platform for annual subfolder
    Y = date.year
    if Y <= 2011:
        annual_folder = f'asi-AMSRE-{hem}{res}-{Y}/'
    else:
        annual_folder = f'asi-AMSR2-{hem}{res}-{Y}/'

    # search for file
    file = path+annual_folder+filename

    if os.path.isfile(path+annual_folder+filename):
        file = path+annual_folder+filename
    else:
        if not quiet:
            print(f'>>> {filename} not found in {path+annual_folder}')

    if not quiet:
        print(f'>>> opening {file}')

    # open data set
    ds = xr.open_dataset(file)
        
    # Extract variables from dataset
    data = extract_variables(ds)

    # download cell area file
    #------------------------------
    if area:

        # cell area file
        area_file = name_area_file(res, hem, include_url=False)

        # check if local area file exists
        if not os.path.isfile(path+area_file):
            print(f'>>> {area_file} not found in {path}')

        else:     
            # open local area file
            ds_area = xr.open_dataset(path + area_file)
            data['area'] = ds_area['data'].values * units('km^2')

    # download cell coordinate file
    #------------------------------
    if coordinates:

        # cell coordinate (lat/lon) file
        coord_file = name_coordinate_file(res, hem, include_url=False)

        # check if local coord file exists
        if not os.path.isfile(path+coord_file):
            print(f'>>> {coord_file} not found in {path}')

        else:
            # extract lon/lat from coordinate file
            lon, lat = open_coord_file(path+coord_file)
            data['lon'] = lon * units('degreeE')
            data['lat'] = lat * units('degreeN')

    # remove units if desired
    if not include_units:
        for key in data.keys():
            if key not in ['proj', 'ds']:
                data[key] = data[key].magnitude



    return data