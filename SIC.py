from datetime import datetime, timedelta
import xarray as xr
import numpy as np
import cartopy
import cartopy.crs as ccrs
from metpy.units import units
import urllib
import os

from .coordinates import open_coord_file, name_area_file, name_coordinate_file



def extract_variables(ds, res = '6250'):
    """Extract variables from dataset. Works for 3125, 6250, 1000ma2

Latest recorded update:
02-13-2025
    """

    data = {}
    data['proj'] = grab_projection(ds)
    data['ds'] = ds

    data['x'] = ds.x.values * units(ds.x.units)
    data['y'] = ds.y.values * units(ds.y.units)
    data['xx'], data['yy'], = np.meshgrid(data['x'], data['y'])

    if str(res) == '1000ma2':
        data['sic_merged'] = ds.sic_merged.values * units('%')
        data['sic_modis'] = ds.sic_modis.values * units('%')
        data['sic_amsr2'] = ds.sic_amsr2.values * units('%')
        data['unc_sic_merged'] = ds.unc_sic_merged.values * units('%')

    else:
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
    - res: str, resolution of desired file ('3125, '6250', '12500', '25000', '1000MA2)
    - hem: str, hemisphere of desired file ('n', 's')
    - include_url: whether to include remote url in filename (True) or file name only (False)
    
    Latest recorded update:
    02-13-2025
    """

    group = hem+res

    dataset_path = 'https://data.seaice.uni-bremen.de/'

    # create strftime for filename based on date
    if date >= datetime(2012, 7, 2):
        file_strftime = f'amsr2/asi_daygrid_swath/{group}/netcdf/%Y/asi-AMSR2-{group}-%Y%m%d-v5.4.nc'
        
    elif (date >= datetime(2002, 6, 1)) & (date <= datetime(2011, 10, 4)):  
        file_strftime = f'amsre/asi_daygrid_swath/{group}/netcdf/%Y/asi-n6250-%Y%m%d-v5.4.nc'

    if str(res) == '1000ma2':
        # NEED TO DOUBLE CHECK ONLINE PATH - WEBSITE CURRENTLY DOWN!
        file_strftime = f'modis-amsr2/asi_daygrid_swath/nh/1000m/%Y/sic_modis-aqua_amsr2-gcom-w1_merged_{hem}h_1000m_%Y%m%d.nc'

    # construct filename
    if include_url:
        file = dataset_path + date.strftime(file_strftime)
    else:
        file = date.strftime(file_strftime).split('/')[-1]

    return file


def open_remote_file(date, res ='6250', hem ='n', 
                     method = 'xarray', 
                     crop = [0, None, 0, None],
                     coordinates = False,
                     area = False,
                     include_units = False, 
                     quiet = True):
    
    """Use xarray to open remote files from HTTPS server:
    https://data.seaice.uni-bremen.de/

INPUT: 
- date: datetime object for desired file
- res: str, resolution of desired file ('6250' or '3125')
- hem: str, hemisphere of desired file ('n' or 's')
- method: str, method to attempt to laod data 
    ('xarray' for reading into memory only, sometimes has trouble
    'urllib' to temporarily download file to disk and then read with xarray)
- crop: indices along dim1 ("i"), dim2 ("j") to crop to area of interest [ai, bi, aj, bj]
- coordinates: bool, whether or not to download lat/lon coordinate data
- area: bool, whether or not to download cell area data
- include_units: bool, whether or not to return data with units
- quiet: bool, whether or not to supress print statements

OUTPUT:
- data: dictionary with nc data

Latest recorded update:
01-31-2025
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

    # crop data values here (and others below)
    ai, bi, aj, bj = crop[0], crop[1], crop[2], crop[3]
    data['y'] = data['y'][ai:bi]
    data['x'] = data['x'][aj:bj]
    data['xx'] = data['xx'][ai:bi, aj:bj]
    data['yy'] = data['yy'][ai:bi, aj:bj]
    data['sic'] = data['sic'][ai:bi, aj:bj]

    # download cell area file
    #------------------------------
    if area:

        # cell area file
        area_file = name_area_file(res, hem)

        # download area file
        urllib.request.urlretrieve(area_file, 'tmp_area.nc')

        # open area file
        ds_area = xr.open_dataset('tmp_area.nc')
        data['area'] = ds_area['data'].values[ai:bi, aj:bj] * units('km^2')
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
        data['lon'] = lon[ai:bi, aj:bj] * units('degreeE')
        data['lat'] = lat[ai:bi, aj:bj] * units('degreeN')
        os.remove("tmp_coord.hdf")



    # remove units if desired
    if not include_units:
        for key in data.keys():
            if key not in ['proj', 'ds']:
                data[key] = data[key].magnitude



    return data


def open_local_file(date, res = '6250', hem ='n', 
                    main_path = '/Volumes/Seagate_Jewell/KenzieStuff/',
                    crop = [0, None, 0, None],
                     coordinates = False, area = False,
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
- res: str, resolution of desired file ('6250' or '3125' or '1000ma2' (1km MODIS-AMSR2 product))
- hem: str, hemisphere of desired file ('n' or 's')
- main_path: str, directory where files are locally stored 
    - looks for subfolder named UniB-ASI-SIC-{}{} where {} will be replaced with res and hem
- crop: indices along dim1 ("i"), dim2 ("j") to crop to area of interest [ai, bi, aj, bj]
- coordinates: bool, whether or not to download lat/lon coordinate data
- area: bool, whether or not to download cell area data
- include_units: bool, whether or not to return data with units
- quiet: bool, whether or not to supress print statements

OUTPUT:
- data: dictionary with nc data

Latest recorded update:
02-13-2025
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

    # 1 km product
    if str(res) == '1000ma2':
        path = main_path + f'UniB-ASI-modis-amsr2-SIC/' # 1km MODIS-AMSR2 product
        annual_folder = f'{Y}/'

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
    data = extract_variables(ds, res=res)

    # crop data values here (and others below)
    ai, bi, aj, bj = crop[0], crop[1], crop[2], crop[3]
    data['y'] = data['y'][ai:bi]
    data['x'] = data['x'][aj:bj]
    data['xx'] = data['xx'][ai:bi, aj:bj]
    data['yy'] = data['yy'][ai:bi, aj:bj]

    if str(res) == '1000ma2':
        # (for some reason dimension 1 must be reversed to match x/y and lat/lon)
        data['sic_merged'] = data['sic_merged'][::-1,:][ai:bi, aj:bj]
        data['sic_modis'] = data['sic_modis'][::-1,:][ai:bi, aj:bj]
        data['sic_amsr2'] = data['sic_amsr2'][::-1,:][ai:bi, aj:bj]
        data['unc_sic_merged'] = data['unc_sic_merged'][::-1,:][ai:bi, aj:bj]
    else:
        data['sic'] = data['sic'][ai:bi, aj:bj]


    # download cell area file
    #------------------------------

    #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # LOOK FOR THIS WITH 1KM PRODUCT
    #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    if area:

        # cell area file
        area_file = name_area_file(res, hem, include_url=False)

        # check if local area file exists
        if not os.path.isfile(path+area_file):
            print(f'>>> {area_file} not found in {path}')

        else:     
            # open local area file
            ds_area = xr.open_dataset(path + area_file)
            data['area'] = ds_area['data'].values[ai:bi, aj:bj] * units('km^2')

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
            if str(res) == '1000ma2':
                lon[lon<0]+=360
            data['lon'] = lon[ai:bi, aj:bj] * units('degreeE')
            data['lat'] = lat[ai:bi, aj:bj] * units('degreeN')

    # remove units if desired
    if not include_units:
        for key in data.keys():
            if key not in ['proj', 'ds']:
                data[key] = data[key].magnitude



    return data



def calc_meansic_openfreq(dates, crop = [0,None,0,None], open_thresh = 70,
                  res = '6250', hem ='n', sic_key = 'sic',
                  main_path = '/Volumes/Seagate_Jewell/KenzieStuff/',
                  coordinates = False, area = False, quiet=True):
    
    """Use xarray to calculate mean open frequency and sic across dates.
     So far only applies to locally stored SIC files from UniBremen (https://data.seaice.uni-bremen.de/)

INPUT: 
- date: datetime object for desired file
- crop: indices along dim1 ("a"), dim2 ("b") to crop to area of interest [ai, aj, bi, bj]
    default is no cropping
- open_thresh: threshold to indicate open water (default: 70%)
- res: str, resolution of desired file ('6250' or '3125')
- hem: str, hemisphere of desired file ('n' or 's')
- sic_key: str, key for sic data in dictionary (default: 'sic', but needs to be adjusted for 1km product)
- main_path: str, directory where files are locally stored 
    - looks for subfolder named UniB-ASI-SIC-{}{} where {} will be replaced with res and hem
- coordinates: bool, whether or not to download lat/lon coordinate data
- area: bool, whether or not to download cell area data
- quiet: bool, whether or not to supress print statements

OUTPUT:
- data: dictionary with coordinates (if desired), mean sic and open frequency across dates, and missing dates

Latest recorded update:
03-30-2025
    """

    # extract cropping info
    ai, aj = crop[0], crop[1]
    bi, bj = crop[2], crop[3]


    # dictionary to store data
    data = {}

    counter = 0 # counts any files that have opened
    missing_dates = np.array([])
    
    for date in dates:

        # open local sic file
        try:
            sic = open_local_file(date, res = res, hem = hem,
                                         main_path = main_path,
                                         coordinates = coordinates, area = area, 
                                         include_units=False)
            exists = True
            counter+=1

        except:
            exists = False

            # record missing dates
            missing_dates = np.append(missing_dates, date)

        if counter == 1:
            
            # mean
            #-------------------------------------------------
            sic_sum = np.zeros(sic[sic_key][ai:aj, bi:bj].shape)
            open_sum = np.zeros(sic[sic_key][ai:aj, bi:bj].shape)
            non_nan = np.zeros(sic[sic_key][ai:aj, bi:bj].shape)

            # add coordinates or cell areas if desired
            data['xx'] = sic['xx'][ai:aj, bi:bj]
            data['yy'] = sic['yy'][ai:aj, bi:bj]

            # add coordinates or cell areas if desired
            if coordinates:
                data['lon'] = sic['lon'][ai:aj, bi:bj]
                data['lat'] = sic['lat'][ai:aj, bi:bj]
            if area:
                data['area'] = sic['area'][ai:aj, bi:bj]
            
        if exists:

            sic_data = np.copy(sic[sic_key][ai:aj, bi:bj])
            
            # if sic_key != 'sic':
            #     sic_data[sic_data>100] = np.nan
                # print(np.sum(np.isnan(sic_data)))

            # record non-nan values, then replace nans with zeros
            if sic_key == 'sic':
                non_nan += np.isfinite(sic_data).astype(int)
                sic_data[np.isnan(sic_data)] = 0
            else:
                non_nan += (sic_data <= 100).astype(int)
                sic_data[sic_data>100] = 0
                
            sic_sum += sic_data
            open_sum += (sic_data < open_thresh).astype(int)

        # record progress
        if not quiet:
            if counter%100 == 0:

                print(f'[{counter}/{len(dates)}] {date}')

    non_nan[non_nan == 0] = 9999

    sic_mean = sic_sum/non_nan
    open_freq = open_sum/non_nan

    sic_mean[non_nan == 9999] = np.nan
    open_freq[non_nan == 9999] = np.nan
    
    data['sic_mean'] = sic_mean
    data['open_freq'] = open_freq
    data['missing_dates'] = missing_dates

    return data
