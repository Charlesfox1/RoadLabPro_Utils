import geopandas as gpd
from centerline.shp2centerline import *
import time
import datetime
from datetime import datetime
import pandas as pd
import shapely.wkt
import os
import re
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString
from osgeo import ogr
import shapely.ops

#Settings
crs_in = {'init': 'epsg:4326'}   #WGS 84
crs_measure = {'init': 'epsg:32648'} #UTM zone 48N
bufwidth = 0.00025 # width to buffer initial roads
road_too_small_m = 30 # length of roads to delete (metres, crs_measure)
accuracy_of_centerlines = bufwidth / 4 # simplification of polygons for calculating centerline. Seems to work well.
verbose = 0 # Print mid-script actions or not
acceptable_suspension_settings = ['MEDIUM', 'HARD-MEDIUM'] #Filter for Point cloud on suspension
low_speed = 20 #Filter for Point cloud on speed: lowest speed for IRI to be considered valid
high_speed = 40 #Filter for Point cloud on speed: highest speed for IRI to be considered valid
timethresh = r'01/01/2017' #Filter for Point cloud on time: earliest IRI reading to be considered valid
timethresh_format = '%d/%m/%Y'
timethresh2 = datetime.strptime(timethresh, timethresh_format)
VPROMMS_ID_certainty_thresh = 0.75 #When re-attaching VPROMMS_IDs, assign certainty after this percentage of points have the same ID

#Input files
linefile = r'IRIlines_csv_format_n.csv'
path = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_YD\experimental\\'
runtime = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_YD\experimental\\temp\\'
RawPoints = r'FOR_JOIN_PATHS.csv'
RawLines = r'FOR_JOIN_INTS.csv'

#User Defined Functions
def Filedump(df, name):
    if verbose == 1:
        try:
            df.to_csv(runtime+r'%s.csv' % name)
        except:
            pass
    else:
        pass

def FileOut(df, name):
    df.to_csv(path+r'%s.csv' % name, index = False)

def Vprint(s):
    if verbose == 1:
        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print '\n%s -- %s' % (st, s)
    else:
        pass

def opIRI(x):
    y = pd.DataFrame()
    y['ID'] = [x.ID.iloc[0]]
    y['iri_mean'] = [x.iri.mean()]
    y['iri_med'] = [x.iri.median()]
    y['iri_min'] = [x.iri.min()]
    y['iri_max'] = [x.iri.max()]
	y['iri_stdev'] = [x.iri.std()]
    y['speed_mean'] = [x.speed.mean()]
    y['speed_med'] = [x.speed.median()]
    return y

def opVID(x):
    z = pd.DataFrame()
    z['ID'] = [x.ID.iloc[0]]
    y = pd.DataFrame(x['VPROMMS_ID'].str[:10].value_counts(normalize = True).reset_index())
    y.columns = ['VID','Fraction']
    a, b = y.VID.tolist(), y.Fraction.tolist()
    if b[0] > VPROMMS_ID_certainty_thresh:
        d = a[0]
    else:
        d = 'Input VPROMMS: '
        for x in range(0,len(a)):
            d = d + '%s (%s pct); ' % (a[x], round((b[x]*100),2))
    z['VPROMMS_ID'] = d
    return z

################################################################################
Vprint('**Phase 1: Create snapped road network**')

Vprint('Read in Roads file, buffer all lines')
roadlines = pd.read_csv(path+linefile)
gdf_lines = gpd.GeoDataFrame(roadlines, crs=crs_in, geometry = roadlines['LINER'].map(shapely.wkt.loads))
gdf_lines['Buffer'] = gdf_lines['geometry'].buffer(bufwidth)
gdf_lines = gdf_lines.set_geometry('Buffer')
gdf_lines = gdf_lines.drop('geometry', axis = 1)
union_obj = str(gdf_lines.unary_union)

Vprint('join roads into single multipolygon, then find centerline')
union_obj = union_obj.replace('MULTIPOLYGON (','').replace(')))','))')
pat = r'(?<=\(\().+?(?=\)\))'
match = re.findall(pat, union_obj)
sausages = pd.DataFrame(match, columns = ['geometry'])
sausages['geometry'] = 'POLYGON ((' + sausages['geometry'].map(str) + '))'
sausages['id'] = sausages.index
Filedump(sausages,'temp')
sausages = gpd.GeoDataFrame(sausages, crs = crs_in, geometry = sausages['geometry'].map(shapely.wkt.loads))
sausages.to_file(os.path.join(runtime, 'temp.shp'), driver = 'ESRI Shapefile')
tempfile = 'temp.shp'
targ = runtime + tempfile

Shp2centerline(targ, runtime+r'\\output_%s.shp' % accuracy_of_centerlines, accuracy_of_centerlines)

a = gpd.read_file(runtime + r'\\output_%s.shp' % accuracy_of_centerlines)
a['Shapelyobj'] = a['geometry']
a['WKTgeometry'] = a['geometry'].map(str)

Phase1_final_df = pd.DataFrame(columns = ['WKT'])
roads = []
junctions = pd.DataFrame(columns = ['Junction'])

#split up the geometry collections, add useful roads to df
for row in a.index:
    Vprint('splittling collection number %s' % row)
    obj = a['WKTgeometry'].iloc[row]
    obj = obj.replace('MULTILINESTRING (','').replace('))',')')
    pat = r'(?<=\().+?(?=\))'
    match = re.findall(pat, obj)
    b = pd.DataFrame(match, columns = ['geometry'])
    b['Line'] = 'LINESTRING (' + b['geometry'].map(str) + ')'
    b['Point1'] = b['geometry'].str.split(', ').str.get(0)
    b['Point2'] = b['geometry'].str.split(', ').str.get(-1)
    c = b['Point1'].tolist() + b['Point2'].tolist()
    e = pd.DataFrame(c,columns = ['points'])
    f = pd.DataFrame(e.apply(pd.value_counts),columns = ['points'])
    f = f.loc[f['points'] >= 3]
    f['Geom'] = f.index
    f['Junction'] = 'POINT (' + f['Geom'].map(str) + ')'
    junctions = pd.concat([junctions, f], axis = 0, ignore_index=True)
    f['Junction'] = f['Junction'].map(shapely.wkt.loads)
    points = MultiPoint(f['Junction'])
    line = a['Shapelyobj'].iloc[row]
    line = shapely.ops.linemerge(line)
    try:
        splitted = shapely.ops.split(line, points)
        for x in splitted:
            roads.append(str(x))
    except:
        Vprint('        ** This collection cannot be split **')
        roads.append(str(line))

Phase_1_roads = pd.DataFrame(roads)
Phase_1_junctions = junctions.drop(['Geom'], axis = 1)
Filedump(Phase_1_junctions, 'Phase_1_junctions')
Filedump(Phase_1_roads, 'Phase_1_roads')
################################################################################
Vprint('**Phase 2: Deleteing small road artifacts**')
n = 1

def RoadCleanup(inputroads, inputjunctions, n):
    Vprint('Road cleanup: iteration %d' % n)
    bad_roads, good_roads = inputroads, inputroads
    bad_roads.columns, good_roads.columns = ['lines'],['lines']

    Vprint('    Drop Small Roads that do not glue junctions. Threshold length: %dm' % road_too_small_m)
    bad_roads['shapely'] = bad_roads['lines'].map(shapely.wkt.loads)
    gbad_roads = gpd.GeoDataFrame(bad_roads, crs = crs_in, geometry = bad_roads['shapely'])
    gbad_roads = gbad_roads.to_crs(crs_measure)
    gbad_roads['length'] = gbad_roads['geometry'].apply(lambda x: (x.length)).map(float)
    bad_roads = gbad_roads.drop('geometry',axis = 1)
    bad_roads['start'] = bad_roads['shapely'].apply(lambda x: Point(x.coords[0])).map(str)
    bad_roads['end'] = bad_roads['shapely'].apply(lambda x: Point(x.coords[-1])).map(str)
    Vprint('Input Junctions length = %d' % len(inputjunctions['Junction']))
    Vprint('length of bad roads before junction check: %d' % len(bad_roads))
    Filedump(inputjunctions,'Input_Junctions_%d' % n)
    Filedump(bad_roads,'bad_roads_%d' % n)
    inputjunctions['Junction'] = inputjunctions['Junction'].map(str)
    bad_roads = bad_roads.drop(bad_roads[(bad_roads['start'].isin(inputjunctions['Junction'].map(str)) == True) & (bad_roads['end'].isin(inputjunctions['Junction'].map(str)) == True)].index)
    Vprint('Length of bad roads after junction check: %d' % len(bad_roads))
    bad_roads = bad_roads.drop(bad_roads[(bad_roads['length'] > road_too_small_m)].index)

    bad_roads = bad_roads[['lines']]
    good_roads = good_roads.drop(good_roads[(good_roads['lines'].isin(bad_roads['lines']) == True)].index)
    Vprint('    Define new true junctions')
    good_roads['shapely'] = good_roads['lines'].map(shapely.wkt.loads)
    good_roads['length'] = good_roads['shapely'].apply(lambda x: x.length).map(float)
    good_roads['start'] = good_roads['shapely'].apply(lambda x: Point(x.coords[0])).map(str)
    good_roads['end'] = good_roads['shapely'].apply(lambda x: Point(x.coords[-1])).map(str)
    c = good_roads['start'].tolist() + good_roads['end'].tolist()
    e = pd.DataFrame(c,columns = ['points'])
    OutJunctions = pd.DataFrame(e.apply(pd.value_counts))
    OutJunctions['Junction'] = OutJunctions.index
    OutJunctions = OutJunctions.loc[OutJunctions['points'] >= 3]
    OutJunctions['Junction'] = OutJunctions['Junction'].map(shapely.wkt.loads)

    Vprint('    Resplit roads on new junctions')
    Phase2_final_df = pd.DataFrame(columns = ['WKT'])
    points = MultiPoint(OutJunctions['Junction'])
    good_roads['shapely'] = good_roads['lines'].map(shapely.wkt.loads)
    line = MultiLineString(good_roads['shapely'].tolist())
    line = shapely.ops.linemerge(line)
    splitted = shapely.ops.split(line, points)
    roads = []
    for x in splitted:
        roads.append(str(x))
    OutRoads = pd.DataFrame(roads)
    Filedump(OutRoads, 'Outroads_%d'%n)
    Filedump(OutJunctions, 'Outjunctions_%d'%n)
    return OutRoads, OutJunctions

iteration_roads, iteration_junctions = RoadCleanup(Phase_1_roads, Phase_1_junctions, 1)
phase_2_final_roads, phase_2_final_junctions = RoadCleanup(iteration_roads, iteration_junctions, 2)

################################################################################
Vprint('**Phase 3: Re-assign attribute data**')
#Prep Attribute data as centroids of line features
IRI = pd.read_csv(path+RawLines)
gIRI = gpd.GeoDataFrame(IRI, crs = crs_in, geometry = IRI['Line_Geometry'].map(shapely.wkt.loads))
gIRI['centre'] = gIRI['geometry'].centroid
gIRI = gIRI.set_geometry('centre')
gIRI['time'] = pd.to_datetime(gIRI['time'],infer_datetime_format = True)
gIRI = gIRI[['time','speed','suspension','iri','centre']]
gIRI = gIRI.loc[
    (gIRI['time'] > timethresh) &
    (gIRI['speed'] > low_speed) &
    (gIRI['speed'] < high_speed) &
    (gIRI['suspension'].isin(acceptable_suspension_settings))
    ]
gIRI['geometry'] = gIRI['centre']

Vprint('    Prepare the Road network')
Network = phase_2_final_roads
Network['ID'] = Network.index
Network.columns = ['Line_Geometry','ID']
gNetwork = gpd.GeoDataFrame(Network, crs = crs_in, geometry = Network['Line_Geometry'].map(shapely.wkt.loads))
gNetwork['Buffer'] = gNetwork.geometry.apply(lambda g: g.buffer(bufwidth, cap_style=2))
gNetwork['geometry'] = gNetwork['Buffer']

Vprint('    Stitch on to new network the original IRI readings')
JoinIRI = gpd.sjoin(gNetwork, gIRI, how="inner",op='intersects')
JoinIRI = JoinIRI[['ID','Line_Geometry','speed','iri']]
JoinIRI = JoinIRI.groupby(['ID']).apply(lambda x: opIRI(x))
Network = Network.drop(['Buffer','geometry'], axis = 1).merge(JoinIRI, how = 'left', on= 'ID')

Vprint('    Stitch on to new network the original VPROMMS_ID tags')
VID = pd.read_csv(path+RawPoints)
VID = VID[['Point_Geometry','VPROMMS_ID']]
gVID = gpd.GeoDataFrame(VID, crs = crs_in, geometry = VID['Point_Geometry'].map(shapely.wkt.loads))
JoinVID = gpd.sjoin(gNetwork, gVID, how="inner",op='intersects')
JoinVID = JoinVID[['ID','VPROMMS_ID']]
JoinVID = JoinVID.groupby(['ID']).apply(lambda x: opVID(x))
Network = Network.merge(JoinVID, how = 'left', on= 'ID')
gNetwork = gpd.GeoDataFrame(Network, crs = crs_in, geometry = Network['Line_Geometry'].map(shapely.wkt.loads))
gNetwork = gNetwork.to_crs(crs_measure)
gNetwork['length'] = gNetwork['geometry'].apply(lambda x: (x.length)/1000).map(float)

Vprint('    Output the new files')
gNetwork = gNetwork.drop(['geometry'], axis = 1)
FileOut(gNetwork, 'Network')
FileOut(phase_2_final_junctions, 'Junctions')
