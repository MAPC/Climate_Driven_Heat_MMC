# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:24:31 2021

Clean version for github.

@author: cspence
"""

import os
import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
import geopandas as gpd
import matplotlib.pyplot as plt
import gc
import pandas as pd
from scipy import ndimage as nd

from shapely.geometry.multipolygon import MultiPolygon

''' Getting set up: Set parameters, directories, and custom functions.'''

# Set the minimum percent "clear" we want in a single date RPA image to use it
# in the RPA's LST index; and minimum daily high air temperature to use a Landsat date.
minclearp = 10  # Starting with 10% clear minimum.
minTmax = 70    # Starting with 70 degrees Fahrenheit.

# Where are all the downloaded and unzipped Landsat data stored? 
madir = ""
# Where to find the tables we need for the analysis (and write tables we create)
tabledir = ""
# Where to find the vector data we need for the analysis
vectordir = ""
# Where to store intermediate raster files
dirtosave = ""
# Where to save the LST index rasters when created?
resdir = ""


# Three custom functions we need for the analysis
def maxminnorm(array):
    '''
    This function normalizes an indicator value by where it lies between the
    maximum and minimum of all non-NaN values in the array.
    For example, if the max value of the array is 8 and the min value is 0,
    the max-min normed value of 4 would be 0.5.
    
    Inputs: array (array to be normalized)
    Outputs: normedarray: Normalized array, same shape as input array, with 
              NaN values in same places.
    '''
    
    maxval = np.nanmax(array)
    minval = np.nanmin(array)
    normedarray = (array - minval)/(maxval - minval)
    
    return(normedarray)

def fill(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    """
    #import numpy as np
    #import scipy.ndimage as nd

    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]

def mosaicdate(files, date, varname, shapes, rpaname, tiles, dst_crs, dirname):
    '''
    this function opens rasters from a list of file pathsthe list of file names "files"
    whose filenames indicate they were captured on "date" and represent 
    "varname." 
    The function mosaics multiple rasters from the tiles covering the RPA
    region together, then uses the "shapes" parameter (the geometry of the 
    shapefile polygon used to mask the raster) to mask the mosaiced raster to
    a rectangle bounded by the northernmost, southernmost, easternmost, and 
    westernmost extents of the shapefile.
    "dst_crs" is the coordinate reference system of the shapefile. 
    Saves results in "dirname."
    
    Pick files that match "varname" and the date of interest. Read them in and
    add them to a list.
    
    Deletes intermediate rasters.
    
    returns: masked/cropped raster, transform info, number of tiles combined
    to form masked raster, and masked raster's metadata.
    
    Inputs:
      files: List of string-type file paths to all the files downloaded and
              from EarthExplorer and hen unzipped.
      date: The particular image date we are looking for in YYYYMMDD format.
      varname: The default-format name of the variable we are looking for.
          For example, "_ST" for surface temperature and "_PIXELQA" for surface 
          temperature quality assessment code. 
      shapes: SHapely geometry object of the RPA outline. Because of issues 
          with the mask function, needs to be a MultiPolygon even if it's just 
          one polygon.
      rpaname: String name of RPA we are working with for file naming
      tiles: List of LandSat 8 tiles we need to encompass the RPA geometry.
          Format for each tile name shouold be 'hhhvvv' (h = horizontal, v=vertical).
          Left pad with zeroes.
      dst_crs: String name of Coordinate Reference System we want to use in the 
          files we create. 
      dirname: String, name of directory in which to save. 
    
    Returns:
      
    '''
    
    # 1. Create list to store the rasters that match the "date" and "varname"
    #    we're looking for.
    rasdate = list()
    tilelist = list()
    for k in range(len(files)):
        if (files[k][15:23] == date) and (files[k][40:] == (varname + '.tif')):
            tiletemp = files[k][8:14] 
            tilelist.append(tiletemp)
            if tiletemp in tiles:
                ras = rasterio.open(files[k])
                rasdate.append(ras)
            else:
                pass
        else:
            pass
            
    # 2. Check to make sure we found sufficient matching files. If not, skip.
    #    If so, merge them with the function from rasterio. Save a copy of the 
    #    original raster's metadata and update the height, width, and transform 
    #    to suit the new mosaiced raster.
    if len(rasdate) == len(tiles): 
        # Mosaic the rasters
        mosaic, out_trans = merge(rasdate)
        out_meta = ras.meta.copy()
        out_meta.update({"driver": "GTiff",
                         "height": mosaic.shape[1],
                         "width": mosaic.shape[2],
                         "transform": out_trans})
        
        # Write mosaiced raster to file
        out_fp1 = dirname + "\\tiledt_" + date + varname + ".tif"
        with rasterio.open(out_fp1, "w", **out_meta) as dest:
            dest.write(mosaic)
        
        # 3. Read the file we just wrote, and reproject to MA state plane
        #    coordinate system (from the shapefile)
        with rasterio.open(out_fp1) as src:
            transform, width, height = calculate_default_transform(
                src.crs, dst_crs, src.width, src.height, *src.bounds)
            kwargs = src.meta.copy()
            kwargs.update({
                            'crs': dst_crs,
                            'transform': transform,
                            'width': width,
                            'height': height
                            })
            
            # Still with src, create new file path for reprojected version of outfp1.
            proj_fp = out_fp1[:-4] + str('stpln.tif')
            with rasterio.open(proj_fp, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                              source=rasterio.band(src, i),
                              destination=rasterio.band(dst, i),
                              src_transform=src.transform,
                              src_crs=src.crs,
                              dst_transform=transform,
                              dst_crs=dst_crs,
                              resampling=Resampling.nearest)
        
        # 4. Read the file we just wrote, and mask w/ RPA shape. 
        # 4a. We need to change the "nodata" in the original (unmasked) tile 
        #     so that we do not automatically crop it out later (instead we will
        #     deal with it separately).
        rastiledmat = rasterio.open(proj_fp).read(1)
        rastiledmeta = rasterio.open(proj_fp)
        rastiledmat[rastiledmat == -9999] = -999
        rastiledmat = np.expand_dims(rastiledmat, axis=0)
        # now write to new file
        proj_fprevised = out_fp1[:-4] + str('stpln_rev.tif')
        out_meta = rastiledmeta.meta.copy()
        out_meta.update({"driver": "GTiff",
                         "height": rastiledmat.shape[1],
                         "width": rastiledmat.shape[2]})
        
        # Write revised raster to file
        with rasterio.open(proj_fprevised, "w", **out_meta) as dest:
            dest.write(rastiledmat)
        
        rastiled = rasterio.open(proj_fprevised)
        regionmeta = rastiled.meta.copy()
        # Mask with the RPA polygon. set nodata to -9998 (original is -9999) 
        # so we can later distinguish outside-RPA-nodata from missing/blank-nodata
        region_image, region_transform = mask(rastiled, shapes, crop=True, filled=True, nodata=-9998)
        
        out_trans = region_transform
        # Update with latest metadata. 
        regionmeta.update({"driver": "GTiff",
                           "height": region_image.shape[1],
                           "width": region_image.shape[2],
                           "transform": region_transform})
        
        # Write it to a file. We will read it later.
        out_fp = dirname + "\\" + rpaname + "_" + date + varname + "_stF.tif"
        with rasterio.open(out_fp, "w", **regionmeta) as dest:
            dest.write(region_image)
            
    else:
        print('Skipping ' + date + ' ' + varname + ' because not enough tiles')
        out_trans = np.empty([1,2,3])
        out_fp = ''
        regionmeta = ''
    ntiles = len(rasdate)
    
    return(out_fp, out_trans, ntiles, regionmeta)

''' Functions done, set up files '''

os.getcwd()
   
archives = os.listdir(madir)
os.chdir(madir)
    

# Establish wihch Landsat tiles and which weather station we need for each RPA. 
RPAfp = tabledir + "\\RPA_landsat_noaa_lookup.xlsx"
RPAs = pd.read_excel(RPAfp, sheet_name='Edit')
RPAs['RPA'] =  RPAs['RPA'].replace(to_replace = 'Southeast Regional Planning_Economic Development District', value="Southeast Regional Planning & Economic Development District")
RPAs['RPA'] =  RPAs['RPA'].replace(to_replace = 'Nantucket Planning_Economic Development Commission', value="Nantucket Planning & Economic Development Commission")

RPAloc = vectordir + "\\RPA_bounds.shp"
RPApolys = gpd.read_file(RPAloc)
RPApolys['rpa_name'] =  RPApolys['rpa_name'].replace(to_replace = "Belongs to both MAPC & OCPC", value="Metropolitan Area Planning Council")

RPApolys = RPApolys.dissolve(by='rpa_name')

# Get the coordinate reference system (CRS) and store as "shapecrs" for future reference.
shapecrs = RPApolys.crs

# Join the RPA info to RPApolys.
RPApolys = RPApolys.merge(RPAs, left_on = 'rpa_name', right_on = 'RPA')


''' Quick side quest to get air temperature data on the dates we used (from NOAA's GHCDN). 
    Methods with help from https://towardsdatascience.com/getting-weather-data-in-3-easy-steps-8dc10cc5c859.'''

print("Beginning weather station air temperature retrieval")

import requests
import json
from datetime import datetime

airport_ids = RPApolys['NOAA Station ID'].unique()
stationids = [airpt[6:] for airpt in airport_ids]

# read in station point shapefile; adjust to MA stateplane CRS
allstationsfp = tabledir + "\\NOAA_GHCDN_metadata.csv"
allstationmeta = pd.read_csv(allstationsfp)
stationmeta = allstationmeta.loc[allstationmeta['STATIONID'].isin(stationids)]
station_coords = gpd.GeoDataFrame(stationmeta, geometry=gpd.points_from_xy(stationmeta.LONGITUDE, stationmeta.LATITUDE))
# The coordinates have no CRS, so give then WGS84 before resetting to the same CRS as the RPA polygons.
station_coords.crs = 'epsg:4326'
station_coords = station_coords.to_crs(shapecrs)

# Start NOAA download
Token = 'cUYQjnHYxVRfgCCBKegiLWHbEeggLgaO'

datesT = list()
TsF = list()
SID = list()
years = ['2018', '2019', '2020']


#Grab TMax data for each station and year
for idn in airport_ids:
    for year in years:
        year = str(year)
        
        #make the api call
        req = "https://www.ncdc.noaa.gov/cdo-web/api/v2/data?datasetid=GHCND&datatypeid=TMAX&LATITUDE&LONGITUDE&limit=1000&stationid=" + idn + "&startdate="+year+"-04-01&enddate="+year+"-10-31"
        r = requests.get(req, headers={'token':Token})
        
        try:
            #load the api response as a json
            d = json.loads(r.text)
            #get all items in the response which are average temperature readings
            avg_temps = [item for item in d['results'] if item['datatype']=='TMAX']
            #get the date field from all average temperature readings
            datesT += [item['date'] for item in avg_temps]
            #get the actual average temperature from all average temperature readings
            TsF += [item['value'] for item in avg_temps]
            SID += [item['station'] for item in avg_temps]
        except:
            pass

# Create a Pandas (pd) DataFrame (table) to store the information from each date. 
dfTF = pd.DataFrame()

dfTF['station'] = [v for v in SID]
dfTF['date'] = [datetime.strptime(d, "%Y-%m-%dT%H:%M:%S") for d in datesT]
dfTF['maxTemp'] = [float(v)/10.0*1.8 + 32 for v in TsF]
dfTF['dtstr'] = [d[0:4] + d[5:7] + d[8:10] for d in datesT]
dfTF['statshrt']  = [v[6:] for v in SID]

# Check how many stations available on each date. Rule out dates with less than the max.
ndates = len(dfTF['dtstr'].unique())
statdates = dfTF['dtstr'].unique()
datelist = list()
statlist = list()
for k in range(ndates):
    day = dfTF.loc[dfTF['dtstr'] == statdates[k]]
    nstats = len(day['statshrt'].unique())
    datelist.append(statdates[k])
    statlist.append(nstats)
    
lookupdict = {'Date': datelist,
              'Number of stations with data': statlist}
stationdatedf = pd.DataFrame(data=lookupdict)

''' End air temperature data side quest '''


''' Begin raster processing! Set up tile processing: Read .tif files by type '''
# one for each with the same 0:40
# one for [-9:] = 'PIXELQA.tif' (pixel quality assessment- fill, clear vs. snow/water/cloud/cloud shadow/terrain occlusion; cloud conf. + cirrus confidence. We pick just the "clear" pixels.)
# one for [-9:] = '_STQA.tif' (surface temperature quality assessment; Kelvin*100, stdv. of emissivity + dist to cloud. We remove stqa > 10F.)
# one for [-7:] = '_ST.tif' (surface temperature) We remove fill, uncertainty > 10F, non-clear pixels. Filter to RP region. 

print('Beginning Landsat Raster processing')

# Setup for all RPAs: reate a unique list of tile/images
names = list()
for k in range(len(archives)):
    name = archives[k][:40]
    if name in names: 
        pass
    else:
        names.append(name)
        
ntiles = len(names)
  
# Create a list of unique dates 
dates = list()
for k in range(len(names)):
    date = names[k][15:23]  # Acquisition date
    if date in dates:
        pass
    else: 
        dates.append(date)
  
# Create a dictionary with the Landsat dates, then convert to pandas dataframe.
datedict = {'Landsat Date': dates}
lsdf = pd.DataFrame(data=datedict)
# Write the dataframe  to a csv file.
fp = tabledir + "\\landsatdates_ma.csv"
lsdf.to_csv(fp)
# Expand the dataframe to include some other attributes.
lsdf = lsdf.join(dfTF.set_index('dtstr'), on='Landsat Date')


''' Mosaic the tiles, clip to RPA region, adjust values, then write to new file.'''
# # Initiate lists in which we will store info about each mosaiced and masked raster.
# # We'll do this for all the dates, then filter down to those with sufficiently high-quality data for our purposes.

# For each date
for k in range(len(RPApolys)):
    # Set up the things we need to merge rasters from multiple tiles and mask the
    # merged rasters to the RPA shape.
    shapes = RPApolys["geometry"].values[k]
    RPAname = RPApolys['RPAshrt'].values[k]
    tiles = eval(RPApolys['Tiles'].values[k])
    station = RPApolys['NOAA Station ID'].values[k].strip()
    print('Starting ' + RPApolys['RPA'].values[k])
    
    # Select just the dates where the weather station data says it was warm enough
    # for our purposes.
    lsdf_rpa = lsdf.loc[(lsdf['maxTemp'] >= minTmax) & (lsdf['station'] == station)]
    # Set list to examine to only dates with logan air temp 70F or higher. 
    dates = lsdf_rpa['Landsat Date'].values

    # Initiate lists that we'll fill in as we cycle through the landsat dates 
    # that were sufficiently warm at our station.
    pctclear = list()   # Percent of pixels "clear" on date
    tilenst = list()    # Number of tiles found for date
    minst = list()      # Min surface temp (F) on valid pixels
    maxst = list()      # Max surface temp (F) on valid pixels
    meanst = list()     # Mean surface temp (F) on valid pixels
    medst = list()      # Median surface temp (F) on valid pixels
    year = list()       # Year on which image captured
    month = list()      # Month on which image captured
    day = list()        # Day of month on which image captured
    datesfull = list()  # List of dates (original YYYYMMDD string format)
    clearbool = list()  # Boolean array: Is pixel clear?
    transf = list()     # List of transformation information for each synthesized raster
    fps = list()        # File paths of successfully synthesized ST rasters
    emptyp = list()     # Stores how much of clipped-to-rpa tile is empty
    
    # Make sure the RPA geometry is a Shapely Multipolygon because rasterio.mask needs this.
    if type(shapes) is MultiPolygon: 
        pass
    else:
        shapes = MultiPolygon([shapes])
    
    ngood= 0    # This is a counter we will use to make sure we are initiating
                # things rather than adding to them if this is the first date
                # we've found that meets quality standards for the RPA region.
    for date in dates:
        
        # mosaic rasters (if they exist, outputs empty if not) and mask to RPA region only
        (stmosfp, stinf, sttilen, rasmeta) = mosaicdate(archives, str(date), '_ST', shapes, RPAname, tiles, shapecrs, dirtosave)
        (stqamosfp, stqainf, stqatilen, stqameta) = mosaicdate(archives, str(date), '_STQA', shapes, RPAname, tiles, shapecrs, dirtosave)
        (stpixmosfp, stpixinf, stpixtilen, pixmeta) = mosaicdate(archives, str(date), '_PIXELQA', shapes, RPAname, tiles, shapecrs, dirtosave)
        
        # def mosaicdate(files, date, varname, shapes, rpaname, tiles, dst_crs, dirname):
        # Next: Create boolean array matching dimensinos of stpixmos indicating where image is clear.
        # Get binary value as string: '{0:b}'.format(900) for example
        # pixelQA = 322: Clear, not fill, cloud conf low. cirrus conf low
        # If rasters exist, read in, remove obscured & high-uncertainty pixels, and convert to Fahrenheit
        if len(stmosfp) > 0:
            
            stmos = (rasterio.open(stmosfp).read(1))
            stqamos = (rasterio.open(stqamosfp).read(1))
            stpixmos = (rasterio.open(stpixmosfp).read(1))
            
            # Sad :/ but this is how I did it. "clear" bit is second from right; 1 if clear, 0 if not. 
            # ("Clear": not cloud, water, snow, cloud shadow, terrain occlusion)
            stpixmos[stpixmos > 2047] = 1   # Set quality band to "fill" if invalid
            clear = np.logical_and(((((((((((stpixmos % 1024) % 512) % 256) % 128) % 64) % 32) % 16) % 8) % 4) % 2) < (((((((((stpixmos % 1024) % 512) % 256) % 128) % 64) % 32) % 16) % 8) % 4), stpixmos > 1)
            water = np.logical_and((((((((((stpixmos % 1024) % 512) % 256) % 128) % 64) % 32) % 16) % 8) % 4) < ((((((((stpixmos % 1024) % 512) % 256) % 128) % 64) % 32) % 16) % 8), stpixmos > 1)
            terrocc = np.logical_and((stpixmos % 1024), stpixmos > 1)
            # make an exception for water
            clear[water] = True
            del stpixmos    # Deleting this to free space because we are done with it at this point.
            
            # If the resulting raster is an expected size, covered by the right 
            # number of tiles, and has at least some clear areas, process further.
            if (np.sum(np.shape(stmos)) > 6 and sttilen == len(tiles)) and np.sum(clear) > 0:
                ngood = ngood + 1
                
                # Translate surface temp quality assessment (Kelvin * 100) to Fahrenheit, remove fill
                    # First, calculate pct ocean. Ocean (and outside region shape)= -9998, all other original "nodata" = -9999.
                    # Turn this into an RPA mask. 
                dims = np.shape(stqamos)
                nonrpamask = np.ones(dims)
                nonrpamask[stqamos == -9998] = 0
                
                pctempty = np.sum(stqamos < 0)/(np.shape(stqamos)[0]*np.shape(stqamos)[1])
                emptyp.append(pctempty)
                print(str(date) + ": " + str(pctempty*100) + "% empty")
                if pctempty <= np.min(emptyp): maskrpa = nonrpamask    # Store the version with least "empty" pixels
                
                stqamos = stqamos.astype(float)
                stqamos[stqamos == -9998] = np.nan
                stqamos[stqamos < 0] = np.nan
                stqamosF = np.squeeze(1.8*(stqamos/100.0))
                del stqamos
                
                # Process surface temperature itself:
                #   Remove any pixels that are fill and that aren't "clear" 
                stmos = stmos.astype(float)
                stmos[stmos < 0] = np.nan
                stmos[np.isnan(stqamosF)] = np.nan
                stmos[~clear] = np.nan
                # Convert surface temp to degrees Fahrenheit
                stmosF = np.squeeze(1.8*((stmos/10.0) - 273.15) + 32)
                del stmos
                
                # Remove any surface temp pixels where uncertainty > 10 degrees Fahrenheit
                stmosF[stqamosF > 10] = np.nan
                del stqamosF
                # Add another dimension (band 1)
                stmosF = np.expand_dims(stmosF, axis=0)
                if ngood == 1: stmosF_norm = np.empty(np.shape(stmosF) + (1,))
                
                if np.sum(~np.isnan(stmosF.flatten())) > 0:
                    
                    # write to a new file.
                    # Update metadata
                    stmetasrc = rasterio.open(stmosfp)
                    out_meta = stmetasrc.meta.copy()
                    out_meta.update({"driver": "GTiff",
                                     "dtype": "float64",
                                     "height": stmosF.shape[0],
                                     "width": stmosF.shape[1],
                                     "transform": stinf})
                    out_fp = dirtosave + "//flt_" + str(date) + "_stF.tif"
                    with rasterio.open(out_fp, "w", **out_meta) as dest:
                        dest.write(stmosF)
                        
                    # # Calculate percent of mosaiced tiles that are clear, + stats of remaining cell values.
                    pclear = (np.sum(~np.isnan(stmosF[0,maskrpa.astype(bool)]))/np.product(np.shape(stmosF)))*100
                    
                    gc.collect()
                    
                    # Create maxmin normed array: [1, nhpixels, nvpixels, ndays]
                    if (pclear > minclearp) & (np.max(stmosF_norm) < 0.01): 
                        print("Haven't yet filled out normed one-date image (this is the first good date)")
                        stmosF_norm = np.empty(np.shape(stmosF) + (1,))
                        stnorm = maxminnorm(stmosF)
                        
                        stmosF_norm[:,:,:,0] = stnorm
                        
                        # Save the file path for later
                        fps.append(out_fp)
                        
                        # Document the image capture dates
                        datesfull.append(str(date))
                        year.append(str(date)[0:4])
                        month.append(str(date)[4:6])
                        day.append(str(date)[6:])
                        
                        pctclear.append(pclear)
                        
                        # Other dist. stats
                        maxst.append(np.nanmax(stmosF))
                        minst.append(np.nanmin(stmosF))
                        meanst.append(np.nanmean(stmosF))
                        medst.append(np.median(stmosF))
                        
                        print('Max T: ' + str(np.nanmax(stmosF)))
                        print('Min T: ' + str(np.nanmin(stmosF)))
                        
                        tilenst.append(sttilen) 
                    elif pclear > minclearp:
                        stnorm = maxminnorm(stmosF)
                        stmosF_norm = np.append(stmosF_norm, np.expand_dims(stnorm, axis=3), axis = 3)
                        # Save the file path for later
                        fps.append(out_fp)
                        
                        # Document the image capture dates
                        datesfull.append(str(date))
                        year.append(str(date)[0:4])
                        month.append(str(date)[4:6])
                        day.append(str(date)[6:])
                        
                        pctclear.append(pclear)
                        
                        # Other dist. stats
                        maxst.append(np.nanmax(stmosF))
                        minst.append(np.nanmin(stmosF))
                        meanst.append(np.nanmean(stmosF))
                        medst.append(np.median(stmosF))
                        
                        print('Max T: ' + str(np.nanmax(stmosF)))
                        print('Min T: ' + str(np.nanmin(stmosF)))
                        
                        tilenst.append(sttilen) 
                    else:
                        pass
                    del stmosF
                else:
                    print('Skipping ' + str(date) + ' because 0% clear')
            else:
                print('Skipping ' + str(date))
        else:
            pass
    
    # With all dates processed, create some metadata from the lists we collected along the way. 
    # Store as a dictionary, then convert to pandas Dataframe. Add the Logan air temperature info 
    # we collected earlier by joining the two dataframes by the "date" attribute.   
    dateinfo = {"Date": datesfull,
                "Year": year,
                "Month": month,
                "Day": day,
                "Number of Tiles": tilenst,
                "Percent Clear": pctclear,
                "Max ST (F)": maxst,
                "Min ST (F)": minst,
                "Mean ST (F)": meanst,
                "Median ST (F)": medst}
    
    # Visualize T in Fahrenheit.
    datedf = pd.DataFrame(data=dateinfo)
    datedf = datedf.join(dfTF.set_index('dtstr'), on='Date')
    datedf = datedf.loc[datedf['station'] == station]
    
    # Send to csv
    metaloc = tabledir + "\\results_info_10F_clear_" + RPAname + ".csv"
    datedf.to_csv(metaloc)
    
    ''' At last: Create the LST index from the max-min normed single-day surface temperature images.'''
    # First, make sure we use the correct water mask.
    date = datesfull[0]
    (stpixmosfp, stpixinf, stpixtilen, pixmeta) = mosaicdate(archives, str(date), '_PIXELQA', shapes, RPAname, tiles, shapecrs, dirtosave)
    stpixmos = (rasterio.open(stpixmosfp).read(1))
    stpixmos[stpixmos > 2047] = 1   # Set quality band to "fill" if invalid
    water = np.logical_and((((((((((stpixmos % 1024) % 512) % 256) % 128) % 64) % 32) % 16) % 8) % 4) < ((((((((stpixmos % 1024) % 512) % 256) % 128) % 64) % 32) % 16) % 8), stpixmos > 1)
    
    maskrpa = np.array(maskrpa, dtype = bool)
    st_norm_mn = np.nanmean(stmosF_norm, axis=3)
    st_norm_mn = fill(st_norm_mn)   # Fill in missing data with nearest values
    st_norm_mn[:,~maskrpa] = np.nan # Then make sure the parts outside the RPA are nodata
    #st_norm_mn[:,water] = np.nan
    st_norm_std = np.nanstd(stmosF_norm, axis=3)
    st_norm_std[:, ~maskrpa] = np.nan
    #st_norm_std[:,water] = np.nan
    
    indexdest = resdir + "\\LST_index_normwater_" + RPAname + ".tif"
    # Copy metadata from last raster processed
    stmetasrc = rasterio.open(fps[-1])
    out_meta = stmetasrc.meta.copy()
    out_meta.update({"driver": "GTiff",
                     "height": st_norm_mn.shape[1],
                     "width": st_norm_mn.shape[2]})
              
    with rasterio.open(indexdest, "w", **out_meta) as dest:
        dest.write(st_norm_mn)
        
        
    vardest = resdir + "\\LST_std_normwater_" + RPAname + ".tif"
    stmetasrc = rasterio.open(fps[-2])
    out_meta = stmetasrc.meta.copy()
    out_meta.update({"driver": "GTiff",
                     "height": st_norm_std.shape[1],
                     "width": st_norm_std.shape[2]})
               
    with rasterio.open(vardest, "w", **out_meta) as dest:
        dest.write(st_norm_std)
    
    # Clear leftover, now invalid data before moving on to next RPA
    del water
    del nonrpamask
    del stmosF_norm
    del stpixmos
    del maskrpa
    # End of RPA loop, move on to next RPA.

# End RPA loop; all done. 
print('Completed LST index for all RPAs.')

''' Merge individual RPA rasters into a single dataset ''' 

print('Merging RPA rasters into MA raster.')
# Get RPA LST index and single-date index variability rasters for each RPA, 
# append to a list of file paths for each. 
outras = list()
outfps = list()
outvar = list()
for k in range(len(RPApolys)):
    RPAname = RPApolys['RPAshrt'].values[k]
    
    indexfp = resdir + "\\LST_index_normwater_" + RPAname + ".tif"
    outfps.append(rasterio.open(indexfp))
    
    vardest = resdir + "\\LST_std_normwater_" + RPAname + ".tif"
    outvar.append(rasterio.open(vardest))
    
# Merge all into one re-scaled raster
merge_indexfp = resdir + "\\Statewide_LST_RPA.tif"
statewide_index_output = merge(outfps, dst_path = merge_indexfp, nodata=np.nan)

merge_varfp = resdir + "\\Statewide_LST_var_RPA.tif"
statewide_var_output = merge(outvar, dst_path = merge_varfp, nodata=np.nan)

# Now open the (LST index) result and plot it!
statewide_index = rasterio.open(merge_indexfp).read(1)

fig, ax = plt.subplots(figsize=(10,10))
lstplot = ax.imshow(np.squeeze(statewide_index), cmap = 'gist_rainbow', vmin=0, vmax = 1)
ax.set_title("LST index", fontsize=14)
cbar = fig.colorbar(lstplot, fraction=0.035, pad=0.01)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel("Normalized ST", rotation=270)
ax.set_axis_off()
plt.show()
    

