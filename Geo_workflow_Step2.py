## LIBRARY IMPORTS ##
import pandas as pd
import shapely
import gdal
import fiona
import sys, os
from shapely.geometry import Point, LineString, MultiPoint
from shapely.wkt import loads
import geopandas as gpd
import rtree
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns
import pyproj
import numpy as np

## USER INPUTS ##
graphs = 'y'
crs_in = {'init': 'epsg:4326'}   #WGS84
crs_out = {'init': 'epsg:4326'} #{'init':'epsg:32648'}   #UTM48N
bufwidth = 0.0003 #CONVERT TO METRES
RoadlabPoints = "FOR_JOIN_PATHS.csv"
RoadlabIntervals = "FOR_JOIN_INTS.csv"
geod = pyproj.Geod(ellps='WGS84')  #for calculating distance measurements
print '\n## Welcome to the Data Manipulation Script for Roadlab Pro ##'
print '\nInitiating Data import. Data will be imported in projection %s and exported in %s projection' % (crs_in['init'], crs_out['init'])
district = str(raw_input('\nDistrict Code: (YD | TT) '))
Path = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_%s' % district
print '\n\n## Commence manipulations ##'

## FUNCTIONS ##
def GDFdescriber(dataframe):
    print 'Data frame sample:'
    print 'Geometry = ', dataframe.geometry.name,' crs =',dataframe.crs, dataframe.head(5), dataframe.dtypes

def plotterhist(Z,z,R,r):
    sns.distplot(Z, label = z)
    sns.distplot(R, label = r)
    plt.legend()
    plt.savefig("%s\\%s_%s_hist_%s.png" % (Path, z, r, bufwidth))

def converter(x):
    if pd.isnull(x['prevgeo']) == True:
        x['Distance_meters'] = 0
        return x
    else:
        x['prevgeo'] = loads(x['prevgeo'])
        x['Point_Geometry'] = loads(x['Point_Geometry'])
        angle1,angle2,distance = geod.inv(x['prevgeo'].x, x['prevgeo'].y, x['Point_Geometry'].x, x['Point_Geometry'].y)
        x['Distance_meters'] = distance
        return x

def LINEGROUPER(x2):
    y = pd.DataFrame()
    try:
        y['LINER'] = [LineString(x2.geometry.tolist())]
    except:
        y['LINER'] = 'null'
    y['iri_mean'] = [x2.iri.mean()]
    y['iri_med'] = [x2.iri.median()]
    y['iri_min'] = [x2.iri.min()]
    y['iri_max'] = [x2.iri.max()]
    y['iri_StDev'] = [x2.iri.std()]
    y['speed_mean'] = [x2.speed.mean()]
    y['speed_med'] = [x2.speed.median()]
    y['speed_min'] = [x2.speed.min()]
    y['speed_max'] = [x2.speed.max()]
    y['length'] = [x2.Distance_meters.sum(skipna=True)]
    y['VPROMMS_ID'] = [x2.VPROMMS_ID.iloc[0]]
    y['npoints'] = len(x2)
    return y

def POINTGROUPER(x2):
    if x2['cerror'].iloc[0] == 1:
        pass
    else:
        x2['VPROMMS_ID'] = "%s_seg%s" % (x2['VPROMMS_ID'].iloc[0], x2['cerror'].iloc[0])
    x2['prevgeo'] = x2.Point_Geometry.shift(1)
    x2['RelativeIRI'] = x2['iri']*(10/x2.iri.max())
    return x2

def DISTCALC(x3):
    x3 = x3.apply(lambda x: converter(x), axis = 1)
    x3['Cm_Distance_metres'] = x3.Distance_meters.cumsum(axis =0)
    return x3

def GROUPER(x,b):
    x = x.sort_values(by='time')
    if b == 'line':
        x = x.groupby(['cerror']).apply(lambda x: LINEGROUPER(x))
    elif b =='point':
        x['cerror'] = x['error'].cumsum(axis=0)+1
        x = x.groupby(['cerror']).apply(lambda x: POINTGROUPER(x))
        x = x.groupby(['VPROMMS_ID']).apply(lambda x: DISTCALC(x))
    else:
        pass
    return x

## IMPORT DATA AS GEODATAFRAMES ##
#point data set first
print '\nStep 1 : Importing original point dataframe'
dfp = pd.DataFrame(pd.read_csv(os.path.join(Path, RoadlabPoints)))
geometryp = [Point(xy) for xy in zip(dfp.longitude, dfp.latitude)]
gdfp = gpd.GeoDataFrame(dfp, crs=crs_in, geometry=geometryp)
#gdfp = gdfp.rename(columns={'geometry': 'points'}).set_geometry('points')
pcount = gdfp['geometry'].count()
#print 'identified %d points in supplied point frame' % pcount
#GDFdescriber(gdfp)

#now line data set
print '\nStep 2 : Importing original line dataframe'
dfl = pd.DataFrame(pd.read_csv(os.path.join(Path, RoadlabIntervals))).drop('Line_Geometry',1)
geometryl = dfl.apply(lambda row: LineString([(row["start_lon"],row["start_lat"]),(row["end_lon"],row["end_lat"])]), axis=1)
gdfl = gpd.GeoDataFrame(dfl, crs=crs_in, geometry=geometryl).drop('is_fixed',1)
#gdfl = gdfl.rename(columns={'geometry': 'lines'}).set_geometry('lines')
lcount = gdfl['geometry'].count()
#print 'identified %d lines in supplied point frame' % lcount
#GDFdescriber(gdfl)

## MANIPULATIONS ##
#Add buffer to all points
print "\nStep 3: Applying buffer to all points, buffer width of %f" % bufwidth
gdfbuffer = gdfp
gdfbuffer['geometry'] = gdfbuffer['geometry'].buffer(bufwidth)
gdfbuffer = gdfbuffer.set_geometry('geometry')
#print "\nSample of 5 buffers:\n"
#GDFdescriber(gdfbuffer)

#Execute spatial join on buffers to points
print '\nStep 4: executing spatial join - lines intersecting buffers\n'
gdfl = gdfl[gdfl.is_valid]
#gdfbuffer.to_file(os.path.join(Path, 'InputBUFFER.shp'), driver = 'ESRI Shapefile') drop bool fields first
points2 = gpd.sjoin(gdfbuffer,gdfl,how="inner",op='intersects')
points2 = points2[points2['VPROMMS_ID_left'] == points2['VPROMMS_ID_right']]
#GDFdescriber(points2)

#rationalise multiple buffers per point
print '\nStep 5a: rationalise multiple buffers on same point; average IRIs for each buffer'
points2 = points2.filter(['latitude','longitude','iri','speed'],axis =1)
points2['Point_Geometry'] = 'POINT ('+points2['longitude'].map(str)+' '+points2['latitude'].map(str)+')'
points2 = points2.groupby('Point_Geometry').mean()
points2['Point_Geometry'] = 'POINT ('+points2['longitude'].map(str)+' '+points2['latitude'].map(str)+')'
points2 = points2.drop(['latitude','longitude'], axis =1)
IRIpoints = pd.DataFrame(pd.read_csv(os.path.join(Path, RoadlabPoints))).merge(points2, how = 'left', on='Point_Geometry')

#Add distance between each point in metres
print '\nStep 5b: Calculate real-world distance between points'
IRIpoints = IRIpoints.groupby(['VPROMMS_ID']).apply(lambda x: GROUPER(x, 'point'))
IRIpoints = IRIpoints.drop(['prevgeo', 'Point_Geometry', 'Unnamed: 0'], axis =1)

#Make those points a line
print '\nStep 6: Join VPROMMS_ID segments, breaking for errors:'
geop = [Point(xy) for xy in zip(IRIpoints.longitude, IRIpoints.latitude)]
IRIpoints2 = gpd.GeoDataFrame(IRIpoints, geometry=geop, crs = crs_in)
IRIlines = IRIpoints2.groupby(['VPROMMS_ID']).apply(lambda x: GROUPER(x, 'line'))
IRIlines = IRIlines.loc[IRIlines['npoints'] > 4]
IRIpoints = IRIpoints.loc[IRIpoints['VPROMMS_ID'].isin(IRIlines['VPROMMS_ID']) == True]

## OUTPUTS ##
#Input RLPlines as shapefile
print 'Step 7: Begin output of files:'
print "\nPre-write file format: RLP input Line frame:\n"
#GDFdescriber(gdfl)
gdfl.to_file(os.path.join(Path, 'InputLinesRLP.shp'), driver = 'ESRI Shapefile')

#Output new lines as CSV
print '\nSending lines to CSV...'
IRIlines.to_csv(os.path.join(Path, 'IRIlines_csv_format.csv'), index = False)

#Output new lines as shapefile
print "\nPre-write file format: Point frame to SHP:\n"
IRIlines2 = gpd.GeoDataFrame(IRIlines, geometry=IRIlines['LINER'], crs=crs_in)
IRIlines2 = IRIlines2.to_crs(crs_out)
IRIlines2 = IRIlines2.drop('LINER', axis = 1)
GDFdescriber(IRIlines2)
IRIlines2.to_file(os.path.join(Path, 'OutputLines.shp'), driver = 'ESRI Shapefile')

#Output  points file now with IRI as CSV
print "\nSending point frame to CSV...\n"
IRIpoints.to_csv(os.path.join(Path, 'IRIpoints_csv_format.csv'), index = False)

#Output original points file now with IRI as shapefile
print "\nPre-write file format: Point frame to SHP:\n"
IRIpoints2 = IRIpoints2.drop(['error','disterror','timeerror'], axis =1)
IRIpoints2 = IRIpoints2.to_crs(crs_out)
#GDFdescriber(IRIpoints2)
IRIpoints2.to_file(os.path.join(Path, 'OutputPoints.shp'), driver = 'ESRI Shapefile')

#Graphs
if graphs == 'y':
    OriginalplotIRI = dfl['iri'].dropna()
    Original_IRI = OriginalplotIRI.astype('int')
    Post_IRI = IRIlines['iri_med'].dropna()
    Post_IRI = Post_IRI.astype('int')
    plotterhist(Original_IRI,'Original IRI', Post_IRI, 'Post S-Join IRI')
else:
    pass

print 'fraction with IRI: %f' % float((IRIpoints['iri'].count() / float(IRIpoints['time'].count())))
print '\nAll exports successful :) Have a great day! \n\n ------END SCRIPT------- \n'
