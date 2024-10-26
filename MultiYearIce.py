from datetime import datetime, timedelta
import xarray as xr
import numpy as np
import cartopy
import cartopy.crs as ccrs
from metpy.units import units
import urllib
import os

# import SIC
from .coordinates import open_coord_file, name_area_file, name_coordinate_file

def grab_projection(ds, quiet = True):

    """Grab projection info from Uni B MYI concentration data (doi: 10.1029/2005JC003384)
    https://seaice.uni-bremen.de/sea-ice-concentration/amsre-amsr2/

INPUT: 
- ds: data opened with xarray
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- projection: cartopy projection from data projection info

Latest recorded update:
10-25-2024
    """
    
    # functino to grab var from string
    def grab_var(text, var: str):
        value = (text[text.find(var):].split(' ')[0]).split('=')[1]
        return value

    # grab parameters from crs spatial attributes
    central_meridian = float(grab_var(ds.comment, 'lon_0'))
    standard_parallel = float(grab_var(ds.comment, 'lat_ts'))
    semimajor = float(grab_var(ds.comment, '+a'))
    semiminor = float(grab_var(ds.comment, '+b'))
    
    if not quiet:
        print(f'>>> data provided in polar_stereographic projection')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - semi_minor_axis: {semiminor}')
        print(f'  - straight_vertical_longitude_from_pole: {central_meridian}')
        print(f'  - standard_parallel: {standard_parallel}')
        print(f'  - proj4text: {ds.attrs["proj4string"]}')
        
    # create projection from info
    projection = ccrs.NorthPolarStereo(central_longitude=central_meridian, 
                                           globe=ccrs.Globe(semimajor_axis = semimajor,
                                                            semiminor_axis = semiminor),
                                           true_scale_latitude=standard_parallel)

    return projection

def extract_variables(ds):
    """Extract variables from dataset."""
    data = {}
    data['ds'] = ds
    data['proj'] = grab_projection(ds)
    data['MYI'] = ds['MYI'].values * units('%')
    
    return data


def name_MYI_file(date, res, hem):
    """Generate the name of the MYI file for a given date and hemisphere."""

    dataset_path = 'https://data.seaice.uni-bremen.de/MultiYearIce/'

    group = hem+res

    # create strftime for filename based on date
    #--------------------------------------
    # 2012/09/22 - present
    if date >= datetime(2012, 9, 22):
        file_strftime = 'ascat-amsr2/final/Arctic/netcdf/%Y/MultiYearIce-Arctic-%Y%m%dv1.1.nc'

    # 2009/09/26 - 2011/05/09 
    elif (date >= datetime(2009, 9, 26)) & (date <= datetime(2011, 5, 9)):  
        file_strftime = 'ascat-amsre/final/netcdf/%Y/MultiYearIce-Arctic-%Y%m%dv1.0.nc'
        
    # I THINK THESE FILES ARE NOT YET PREPARED?
    # 2002/10/01 - 2009/04/30 
    elif date <= datetime(2009, 4, 30):
        file_strftime = 'qscat-amsre/final/Arctic/netcdf/%Y/MYI-NMYI-CORRECTION-%Y%m%d.nc'

    # construct filename
    file = dataset_path + date.strftime(file_strftime) 

    return file


def open_remote_file(date, hem : str='n', 
                     method : str = 'xarray', 
                     coordinates: bool = False,
                     area: bool = False,
                     include_units = False, quiet=True):
    
    """Use xarray to open remote files from HTTPS server:
    https://data.seaice.uni-bremen.de/

INPUT: 
- date: datetime object for desired file
- hem: str, hemisphere of desired file ('n' or 's')
    *** FOR NOW ONLY NORTHERN HEMISPHERE IS SUPPORTED ***
- method: str, method to attempt to laod data 
    ('xarray' for reading into memory only, sometimes has trouble
    'urllib' to temporarily download file to disk and then read with xarray)
- coordinates: bool, whether or not to download lat/lon coordinate data
- area: bool, whether or not to download cell area data
- include_units: bool, whether or not to return data with units
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- data: dictionary with nc data

Latest recorded update:
10-25-2024
    """

    res = '12500'

    # construct filename
    file = name_MYI_file(date, res, hem)

    if not quiet:
        print(f'>>> opening {file}')

    # use xarray to open remote file
    if method == 'xarray':
        # for some reason need to add this for now
        # see https://github.com/Unidata/netcdf4-python/issues/1043
        extrabit = '#mode=bytes'
        ds = xr.open_dataset(file + extrabit, decode_cf=False)

    elif method == 'urllib':
        urllib.request.urlretrieve(file, "temp.nc")
        ds = xr.open_dataset("temp.nc")
        # remove temporary file
        os.remove("temp.nc")

    # Extract variables from dataset
    data = extract_variables(ds)

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

    # open landmask file to grab corresponding x,y coordinates
    ref_file = 'https://data.seaice.uni-bremen.de/landmasks/landmask_Arctic_12.500km.nc'
    urllib.request.urlretrieve(ref_file, 'tmp_xy.hdf')
    ds_xy = xr.open_dataset('tmp_xy.hdf')
    x = ds_xy['x'].values * units('km').to('m')
    # for some reason y values are reversed relative to other projections
    y = ds_xy['y'].values[::-1] * units('km').to('m')
    data['xx'], data['yy'] = np.meshgrid(x,y)
    os.remove("tmp_xy.hdf")

    
    # remove units if desired
    if not include_units:
        for key in data.keys():
            if key not in ['proj', 'ds']:
                data[key] = data[key].magnitude



    return data

