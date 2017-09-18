import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point, LineString, MultiPoint
import shapely.wkt
import glob2
import re
import sys, os
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns
import pyproj
import time
from datetime import datetime

#Inputs
valid_roads = raw_input('\nRemove all non-VPROMMS_ID roads? y = yes, n = no\n')
district = str(raw_input('\nDistrict Code: (YD | TT) '))
InPath = r'C:\Users\charl\Documents\Vietnam\Fieldwork\\roadlab_bin_%s' % district
OutPath = r'C:\Users\charl\Documents\GitHub\RoadLabPro_Utils\\Runtime\\%s' % district
runtime = r'C:\Users\charl\Documents\GitHub\RoadLabPro_Utils\\Runtime\\%s' % district
sdthresh = 0.5
verbose = 1
crs_in = {'init': 'epsg:4326'}   #WGS84
crs_measure = {'init': 'epsg:32648'} #UTM zone 48N
bufwidth = 0.0003

# Utility Functions
def Filedump(df, name):
    if verbose == 1:
        try:
            df.to_csv(os.path.join(runtime,'%s.csv' % name))
        except:
            pass

def FileOut(df, name):
    df.to_csv(os.path.join(OutPath, '%s.csv' % name), index = False)

def Vprint(s):
    if verbose == 1:
        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print '\n%s -- %s' % (st, s)

# Step 1: Aggregate original intervals
intervals_list = glob2.glob(InPath+'\\'+'**'+'\\*intervals*.csv')
dataframes = []
for file_ in intervals_list:
    df = pd.DataFrame(pd.read_csv(file_))
    df['Input_file'] = str(file_)
    dataframes.append(df)
X = pd.concat(dataframes, ignore_index=True)
X['Line_Geometry'] = 'LINESTRING ('+X['start_lon'].map(str)+' '+X['start_lat'].map(str)+', '+X['end_lon'].map(str)+' '+X['end_lat'].map(str)+')'
X['Part1'] = X['Input_file'].map(str).str.extract(('(bin.*\.csv)'), expand=False)
X['Part1'] = X['Part1'].map(str).str.replace('bin', '').str.replace('\.csv', '')
X['Part1'] = X['Part1'].str.split('\\').str.get(-3)
X['Part2'] = X['Input_file'].map(str).str.extract(('(bin.*\.csv)'),expand=False)
X['Part2'] = X['Part2'].map(str).str.replace('bin', '').str.replace('\.csv', '')
X['Part2'] = X['Part2'].str.split('\\').str.get(-2)
for y in X.Part1.unique():
    miniframe = X.loc[X['Part1'] == y]
    n = 1
    for y2 in miniframe.Part2.unique():
        X.loc[(X['Part2'] == y2)&(X['Part1']==y),'Part1'] = X.loc[(X['Part2'] == y2)&(X['Part1']==y),'Part1']+'_'+str(n)
        n += 1
X['VPROMMS_ID'] = X['Part1']
X = X.drop(['Part2','Part1','Input_file'], axis = 1)
if valid_roads == 'y':
    X = X[X.VPROMMS_ID.str.contains('OTHER*') == False]
else:
   pass
Original_Intervals = X
FileOut(Original_Intervals,'Original_Intervals')

# Step 2: Aggregate original Point files
pointframe_list = glob2.glob(InPath+'\\'+'**'+'\\*RoadPath'+'*.csv')
dataframes = []
for file_ in pointframe_list:
    df = pd.DataFrame(pd.read_csv(file_))
    df['Input_file'] = str(file_)
    dataframes.append(df)
Y = pd.concat(dataframes, ignore_index=True)
Y['Point_Geometry'] = 'POINT ('+Y['longitude'].map(str)+' '+Y['latitude'].map(str)+')'
Y['Part1'] = Y['Input_file'].map(str).str.extract(('(bin.*\.csv)'),expand=False)
Y['Part1'] = Y['Part1'].map(str).str.replace('bin', '').str.replace('\.csv', '')
Y['Part1'] = Y['Part1'].str.split('\\').str.get(-3)
Y['Part2'] = Y['Input_file'].map(str).str.extract(('(bin.*\.csv)'),expand=False)
Y['Part2'] = Y['Part2'].map(str).str.replace('bin', '').str.replace('\.csv', '')
Y['Part2'] = Y['Part2'].str.split('\\').str.get(-2)
for y in Y.Part1.unique():
    miniframe = Y.loc[Y['Part1'] == y]
    n = 1
    for y2 in miniframe.Part2.unique():
        Y.loc[(Y['Part2'] == y2)&(Y['Part1']==y),'Part1'] = Y.loc[(Y['Part2'] == y2)&(Y['Part1']==y),'Part1']+'_'+str(n)
        n += 1
Y['VPROMMS_ID'] = Y['Part1']
Y = Y.drop(['Part2','Part1','Input_file'], axis = 1)
if valid_roads == 'y':
    Y = Y[Y.VPROMMS_ID.str.contains('OTHER*') == False]
else:
    pass

#Step 3: clean points based on outlier distance / time deltas between point pairs
def main(X):
    timeconvert(X)
    deltas(X)
    X['disterror'] = outliers(X['distdiff'].as_matrix(), 'distance')
    X['timeerror'] = outliers(X['timediff'].as_matrix(), 'time')
    X['error'] = ((X['timeerror'] == True) | (X['disterror'] == True))
    X = X.drop(['latdiff','longdiff'], axis = 1)
    return X

def timeconvert(X):
    # Convert RLP timestamp to pandas time object
    X['time'] = pd.to_datetime(X['time'],infer_datetime_format = True)
    try:
        X = X.drop(['Unnamed: 0'], axis =1)
    except:
        pass
    return X

def deltas(X):
    for y in X.VPROMMS_ID.unique():
        MF = X.loc[X['VPROMMS_ID'] == y]
        X.loc[X['VPROMMS_ID'] == y, 'timediff'] = MF.time.diff() / np.timedelta64(1, 's')
        X.loc[X['VPROMMS_ID'] == y, 'latdiff'] = MF.latitude.diff()
        X.loc[X['VPROMMS_ID'] == y, 'longdiff'] = MF.longitude.diff()
    X['distdiff'] = np.sqrt(np.square(X['latdiff'])+np.square(X['longdiff']))

def outliers(points,name,thresh = sdthresh):
    Vprint('---Summary of outlier detection: %s---' % name)
    mean = np.nanmean(points, axis=0)
    Vprint('mean = %s' % mean)
    squarediff = (points - mean)**2
    Vprint('squared differences = %s' % squarediff)
    sumsquarediff = np.nansum(squarediff)
    Vprint('sumsquarediff = %s' % sumsquarediff)
    N = len(points)
    Vprint('N = %s' % N)
    sd_deviation = np.sqrt((sumsquarediff / N))
    Vprint('sd_deviation = %s' % sd_deviation)
    z_score = (points - mean) / sd_deviation
    Vprint(z_score)
    return z_score > thresh

Original_Points = main(Y)
FileOut(Original_Points,'Original_Points')

#Step 4: chain points to lines, attach IRI info
## USER INPUTS ##
geod = pyproj.Geod(ellps='WGS84')  #for calculating distance measurements

## FUNCTIONS ##

def LINEGROUPER(x2):
    y = pd.DataFrame()
    try:
        y['Line_Geometry'] = [LineString(x2.geometry.tolist())]
    except:
        y['Line_Geometry'] = 'null'
    y['iri_mean'] = [x2.iri.mean()]
    y['iri_med'] = [x2.iri.median()]
    y['iri_min'] = [x2.iri.min()]
    y['iri_max'] = [x2.iri.max()]
    y['iri_StDev'] = [x2.iri.std()]
    y['speed_mean'] = [x2.speed.mean()]
    y['speed_med'] = [x2.speed.median()]
    y['speed_min'] = [x2.speed.min()]
    y['speed_max'] = [x2.speed.max()]
    y['VPROMMS_ID'] = [x2.VPROMMS_ID.iloc[0]]
    y['npoints'] = len(x2)
    return y

def POINTGROUPER(x2):
    if x2['cerror'].iloc[0] == 1:
        pass
    else:
        x2['VPROMMS_ID'] = "%s_seg%s" % (x2['VPROMMS_ID'].iloc[0], x2['cerror'].iloc[0])
    return x2

def GROUPER(x,b):
    x = x.sort_values(by='time')
    if b == 'line':
        x = x.groupby(['cerror']).apply(lambda x: LINEGROUPER(x))
    elif b =='point':
        x['cerror'] = x['error'].cumsum(axis=0)+1
        x = x.groupby(['cerror']).apply(lambda x: POINTGROUPER(x))
    else:
        pass
    return x

## IMPORT DATA AS GEODATAFRAMES ##
gdfp = gpd.GeoDataFrame(Original_Points, crs = crs_in, geometry = Original_Points['Point_Geometry'].map(shapely.wkt.loads))
gdfl = gpd.GeoDataFrame(Original_Intervals, crs = crs_in, geometry = Original_Intervals['Line_Geometry'].map(shapely.wkt.loads))

## MANIPULATIONS ##
gdfbuffer = gdfp
gdfbuffer['geometry'] = gdfbuffer['geometry'].buffer(bufwidth)
gdfbuffer = gdfbuffer.set_geometry('geometry')

#Execute spatial join on buffers to points
#gdfbuffer.to_file(os.path.join(Path, 'InputBUFFER.shp'), driver = 'ESRI Shapefile') drop bool fields first
points2 = gpd.sjoin(gdfbuffer,gdfl,how="inner",op='intersects')
points2 = points2[points2['VPROMMS_ID_left'] == points2['VPROMMS_ID_right']]

#rationalise multiple buffers per point
points2 = points2[['latitude','longitude','iri','speed']]
points2['Point_Geometry'] = [Point(xy) for xy in zip(points2.longitude, points2.latitude)]
points2['Point_Geometry'] = points2['Point_Geometry'].astype(str)
points2 = points2.groupby('Point_Geometry').mean()
points2['Point_Geometry'] = [Point(xy) for xy in zip(points2.longitude, points2.latitude)]
points2['Point_Geometry'] = points2['Point_Geometry'].astype(str)
points2 = points2.drop(['latitude','longitude'], axis =1)
IRIpoints = Original_Points.merge(points2, how = 'left', on='Point_Geometry')
IRIpoints = IRIpoints.groupby(['VPROMMS_ID']).apply(lambda x: GROUPER(x, 'point')).drop('geometry', axis = 1)

#Make those points a line
geop = [Point(xy) for xy in zip(IRIpoints.longitude, IRIpoints.latitude)]
IRIpoints2 = gpd.GeoDataFrame(IRIpoints, geometry = geop, crs = crs_in)
IRIlines = IRIpoints2.groupby(['VPROMMS_ID']).apply(lambda x: GROUPER(x, 'line'))
IRIlines = IRIlines.loc[IRIlines['npoints'] > 4]
gIRIlines = gpd.GeoDataFrame(IRIlines, crs = crs_in, geometry = IRIlines['Line_Geometry'])
gIRIlines = gIRIlines.to_crs(crs_measure)
gIRIlines['length'] = gIRIlines['geometry'].apply(lambda x: (x.length)).map(float)
gIRIlines = gIRIlines.drop(['geometry'], axis = 1)
IRIpoints = IRIpoints.loc[IRIpoints['VPROMMS_ID'].isin(gIRIlines['VPROMMS_ID']) == True]

## OUTPUTS ##
FileOut(gIRIlines,'Adj_lines')
FileOut(IRIpoints, 'Adj_points')
print 'fraction with IRI: %f' % float((IRIpoints['iri'].count() / float(IRIpoints['time'].count())))
