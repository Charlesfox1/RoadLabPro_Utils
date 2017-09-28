# -*- coding: utf-8 -*-
###################################################################################################
# Calculate Risk for roads
# Charles Fox and Benjmain Stewart, September 2017
# Purpose: Intersect the risk rasters with a network dataset
###################################################################################################
import os, sys, inspect

import fiona
import rasterio
import rasterstats as rs
import pandas as pd
import geopandas as gpd
import shapely.wkt
import numpy as np

verbose = 1
checkcols = 0

def main(district="YD", curVerbose=0, curCheckCols=0):
    #district = str(raw_input('\nDistrict Code: (YD | TT) '))
    verbose = curVerbose
    checkcols = curCheckCols
    
    path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    hazardPath = os.path.join(path, "PCS/Hazard/HazardData") #Stores the hazard raster data

    outpath = os.path.join(path, 'Outputs','%s' % district)
    runtime = os.path.join(path, 'PCS', 'Hazard', 'Runtime')
    for d in [outpath, runtime]:
        if not os.path.isdir(d):
            os.mkdir(d)

    dash = os.path.join(path, 'PCS','dashboard.xlsx')
    Network = os.path.join(path, 'runtime', '%s' % district,'Network.csv')

    #read CMS output in as GDF
    df = pd.read_csv(Network)
    geometry = df['Line_Geometry'].map(shapely.wkt.loads)
    crs = {'init': 'epsg:4326'}
    gdf = gpd.GeoDataFrame(df, crs = crs, geometry = geometry)

    #import risk layers
    Riskdf = pd.read_excel(dash, sheetname = "HAZARD", index_col = 0)
    
    for RASTER in Riskdf.index:
        print '\nComputing exposure for risk layer: %s' % RASTER
        a = '%s_mean' % RASTER
        try:
            curRasterPath = os.path.join(hazardPath,str(Riskdf['PATH'][RASTER]))
            if Riskdf['REL_ABS'][RASTER] == "Relative":
                zstat = Relative(gdf, a, RASTER,curRasterPath, Riskdf['MIN'][RASTER], Riskdf['MAX'][RASTER])
            elif Riskdf['REL_ABS'][RASTER] == "Absolute":
                zstat = Absolute(gdf, a, RASTER,curRasterPath, Riskdf['ABS_VAL'][RASTER], Riskdf['MIN'][RASTER], Riskdf['MAX'][RASTER])
            else:
                print "\n** ERROR: relative / absolute decision not made for risk layer %s **" % RASTER
            df = df.merge(zstat, left_index=True, right_index=True)
        except Exception as e:
            print "\n** ERROR: Risk Layer %s Did not compute correctly **" % RASTER
            print e.message
    
    #Calulcate Risk Score
    b = []
    for RASTER in Riskdf.index:
        a = df['%s_mean' % RASTER] * Riskdf['WEIGHT'][RASTER]
        a[np.isnan(a)] = 0
        b.append(a)
    df['RISK_SCORE'] = sum(b)
    df['RISK_SCORE'] = ((df['RISK_SCORE'] - df['RISK_SCORE'].min()) / (df['RISK_SCORE'].max() - df['RISK_SCORE'].min()))

    if checkcols == 1:
        pass
    elif checkcols == 0:
        for RASTER in Riskdf.index:
            df = df.drop('%s_mean' % RASTER, axis = 1)
    else:
        raise ERROR          

    #Output
    df = df.drop('geometry',axis = 1)
    df = pd.DataFrame(df)
    df.to_csv(os.path.join(outpath,'risk_output.csv'), index=False)  
    
        
#Define Functions
def make_cmean(minval, maxval):
    def _func(x, minval = minval, maxval = maxval):
        x = x[x <= maxval]
        x = x[x >= minval]
        return (np.mean(x))
    return _func

def Vprint(a):
    if verbose == 1:
        print a
    else:
        pass

def Value_Mask(rastpath, minval, maxval):
    Vprint("Loading: %s, \nfloor value: %f, headlimit value : %f" % (rastpath, minval, maxval))
    targ_ras = rasterio.open(rastpath).read(1)
    x,y = targ_ras.shape[0],targ_ras.shape[1]
    if x < 1000 and y < 1000:
        targ_ras_temp = np.ma.masked_where(targ_ras > maxval, targ_ras)
        targ_ras_masked = np.ma.masked_where(targ_ras_temp < minval, targ_ras_temp)
        mini = targ_ras_masked.min()
        maxi = targ_ras_masked.max()
    else:
        x2, y2 = int(x / 10), int(y / 10)
        minlist, maxlist = [], []
        for i in range (1,11):
            for j in range (1,11):
                temp =  targ_ras[((i-1)*x2):(i*x2),((j-1)*y2):(j*y2)]
                temp = np.ma.masked_where(temp > maxval, temp)
                temp = np.ma.masked_where(temp < minval, temp)
                minlist.append(temp.min())
                maxlist.append(temp.max())
        mini = min(minlist)
        maxi = max(maxlist)
    return mini, maxi

def Relative(gdf, a, raster_name, rastpath, minval, maxval):
    mini, maxi = Value_Mask(rastpath, minval, maxval)
    Vprint("Raster %s will be measured relative to a masked max of %f and a masked min of %f" % (raster_name, maxi, mini))
    cmean = make_cmean(minval, maxval)
    zs = rs.zonal_stats(gdf,rastpath,all_touched=True,stats="mean",add_stats={"cmean":cmean})
    stat = pd.DataFrame(zs)
    stat.columns = [a,'rawmean']
    stat['normed'] = ((stat[a]-mini)/(maxi-mini))
    if checkcols == 1:
        stat.to_csv(os.path.join(runtime,'temp_%s.csv' % raster_name))
    else:
        pass
    stat[a] = stat['normed'].astype(float)
    stat = stat.drop(['rawmean','normed'], axis = 1)
    return stat

def Absolute(gdf, a, raster_name, rastpath,thresh, minval, maxval):
    Vprint("Raster %s will be measured against an absolute threshold of %f" % (raster_name, thresh))
    cmean = make_cmean(minval, maxval)
    zs = rs.zonal_stats(gdf,rastpath,all_touched=True,stats="mean",add_stats={"cmean":cmean})
    stat = pd.DataFrame(zs)
    stat.columns = [a,'rawmean']
    stat['normed'] = (stat[a] > thresh).astype(int)
    if checkcols == 1:
        stat.to_csv(os.path.join(runtime,'temp_%s.csv' % raster_name))
    else:
        pass
    stat[a] = stat['normed'].astype(float)
    stat = stat.drop(['rawmean','normed'], axis = 1)
    return stat

main()