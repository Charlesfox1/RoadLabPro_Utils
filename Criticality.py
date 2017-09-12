#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 15:54:18 2017
@author: yan
Re-written by by C Fox for use in Vietnam
"""
# criticality analysis
#
# get_ipython().magic(u'matplotlib inline')
import networkx as nx
import pandas as pd
import geopandas as gpd
import shapely.geometry.base
import shapely.wkt
import copy
import time
import datetime
from datetime import datetime
import os
import sys
import warnings
warnings.filterwarnings("ignore")

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

# Modules developed by TU Delft team
from network_lib import network_prep as net_p
from network_lib import od_prep as od_p

import geopandas as gp
import numpy as np

verbose = 0
dump = 1
# Path Settings
outpath = r'C:\Users\charl\Documents\Vietnam\PCS\\'
runtime = outpath+r'Criticality\Stack\runtime\\'
inpath = r'C:\Users\charl\Documents\Vietnam\PCS\Criticality\Stack\input\Vietnam'
crs_in = {'init': 'epsg:4326'}   #WGS 84

#Utility Functions
def Filedump(df, name):
    # Dumps file to runtime folder
    if dump == 1:
        df.to_csv(runtime+r'%s.csv' % name, encoding = 'utf-8')

def FileOut(df, name):
    # Dumps file to outpath folder
    df.to_csv(outpath+r'%s.csv' % name, index = False, encoding = 'utf-8')

def Vprint(s):
    if verbose == 1:
        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print '\n%s -- %s' % (st, s)

#Network Prep
inNetwork = pd.read_csv(os.path.join(inpath, 'Network.csv'))
fillvalue = inNetwork['iri_med'].mean()
inNetwork['TC_iri_med'] = inNetwork['iri_med'].fillna(fillvalue)
inNetwork['total_cost'] = inNetwork['TC_iri_med'] * inNetwork['length']
ginNetwork = gpd.GeoDataFrame(inNetwork,crs = crs_in, geometry = inNetwork['Line_Geometry'].map(shapely.wkt.loads))
ginNetwork.to_file(os.path.join(runtime, 'Network.shp'), driver = 'ESRI Shapefile')

#Centroid Prep
SingleNetworkObj = pd.DataFrame([str(ginNetwork.unary_union)], columns = ['geom'])
gSingleNetworkObj = gpd.GeoDataFrame(SingleNetworkObj,crs = crs_in, geometry = SingleNetworkObj['geom'].map(shapely.wkt.loads))
inAdmin = os.path.join(outpath+'Poverty\Poverty_Communes_2009.shp')
ginAdmin = gpd.read_file(os.path.join(outpath+'Poverty\Poverty_Communes_2009.shp'))
ginAdmin = ginAdmin.to_crs(crs_in)
ginAdmin = ginAdmin[['P_EName','D_EName','EN_name','Admin_EN','COMPOPULA','PROCODE04','DISTCODE04','CCode04New','geometry']]
ginAdmin = ginAdmin.rename(columns={'COMPOPULA': 'Pop'})
ginAdmin['ID'] = ginAdmin.index
SelectedAdmins = gpd.sjoin(gSingleNetworkObj, ginAdmin, how="inner",op='intersects')
SelectedAdmins = SelectedAdmins.drop(['geom','geometry'], axis = 1)
ginAdmin = ginAdmin.loc[ginAdmin['ID'].isin(SelectedAdmins['ID']) == True]
ginAdmin['geometry'] = ginAdmin['geometry'].centroid
ginAdmin.to_file(os.path.join(runtime, 'adm_centroids.shp'), driver = 'ESRI Shapefile')

##########Criticality Script #########
Outputs = []   # List of final outputs, each object is a dictionary corresponding to an origin:destination pair
#Inputs
network = os.path.join(runtime,'Network.shp')
origin_points = {
    'name': 'adm_centroids',  # must be string
    'file': os.path.join(runtime, 'adm_centroids.shp'),  # Must be shapefile
    'penalty': 1000,   # USD cost of isolated trip
    'annual trips':50,  # scalar for number of trips from / to this destination type
    'scalar_column': 'Pop' # relative attractiveness scalar for each node
    }
hospitals = {
    'name': 'hospitals',
    'file': os.path.join(runtime, 'Hospitals.shp'),
    'penalty': 10000,
    'annual trips':3,
    'scalar_column': 'Beds'
    }
origin_list = [origin_points]
destination_list = [origin_points,hospitals]

# Prepation of network
gdf_points, gdf_node_pos, gdf = net_p.prepare_centroids_network(origin_points['file'], network)
# Create Networkx MultiGraph object from the GeoDataFrame
G = net_p.gdf_to_simplified_multidigraph(gdf_node_pos, gdf, simplify=False)
# Change the MultiGraph object to Graph object to reduce computation cost
G_tograph = net_p.multigraph_to_graph(G)
# Observe the properties of the Graph object
Vprint('number of disconnected components is: %d' % nx.number_connected_components(G_tograph))
nx.info(G_tograph)
# Take only the largest subgraph with all connected links
len_old = 0
for g in nx.connected_component_subgraphs(G_tograph):
    if len(list(g.edges())) > len_old:
        G1 = g
        len_old = len(list(g.edges()))
G_sub = G1.copy()

Vprint('number of disconnected components is: %d' % nx.number_connected_components(G_sub))
nx.info(G_sub)

# Save the simplified transport network into a GeoDataFrame
gdf_sub = net_p.graph_to_df(G_sub)
blank, gdf_node_pos2, gdf_new = net_p.prepare_newOD(origin_list[0]['file'], gdf_sub)

#Road Network Graph prep
G2_multi = net_p.gdf_to_simplified_multidigraph(gdf_node_pos2, gdf_new, simplify=False)
Filedump(gdf_new, 'Road_Lines'), Filedump(gdf_node_pos2,'Road_Nodes')
G2 = net_p.multigraph_to_graph(G2_multi)
gdf2 = net_p.graph_to_df(G2)
nLink = len(G2.edges())

def PrepSet(point_set):
    """
    Prepares a small df of a given origin / destination set, expressed as 'item : Nearest Node ID'
    """
    Prepared_point_set, gdf_node_pos2, gdf_new = net_p.prepare_newOD(point_set['file'], gdf_sub)
    Prepared_point_set = Prepared_point_set['Node']
    return Prepared_point_set

def Pathfinder(origin, destination):
    """
    find the shortest path for each od to minimize the total travel cost;
    output:
    1) baseCost ($): total travel cost between all OD pairs;
    2) basePath the shortest path between all OD pairs;
    3) baseLength: total travel distance between all OD pairs.
    """
    basePath = [[ [] for d in range(len(destination))] for d in range(len(origin))]
    baseCost = np.zeros((len(origin),len(destination)))
    baseLength = np.zeros((len(origin),len(destination)))
    for o in range(len(origin)):
        for d in range(len(destination)):
            basePath[o][d] = nx.dijkstra_path(G2,origin[o],destination[d], weight = 'total_cost')
            baseCost[o][d] = nx.dijkstra_path_length(G2,origin[o],destination[d], weight = 'total_cost')
            baseLength[o][d] = nx.dijkstra_path_length(G2,origin[o],destination[d], weight = 'length')
    return basePath, baseCost, baseLength

def BreakEdge(origin, destination, penalty, baseCost, name):
    """
    fundamental criticaility analysis: remove each link by changing it to a extremely large value,
    calculate the shortest path, then re-set total cost of traversing link.
    output:
    1.) diff: a 3D matrix describing the change in cost for each 2D O-Dmatrix when a given link is broken
    2.) iso: a 3D matrix describing the cost of isolated journeys for each 2D O-D matrix when a given link is broken
    """
    to_allNode = []
    G = copy.deepcopy(G2)
    cost_disrupt = np.zeros((nLink, len(origin),len(destination)))
    for road in range(len(gdf2)):
        start = time.clock()
        w = G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost']
        G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost'] = 1e10
        for o in range(len(origin)):
            to_allNode.append(nx.single_source_dijkstra_path_length(G,origin[o],weight = 'total_cost'))
        G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost'] = w # Resetting traverse cost to original (w)
    for item_number in range(len(to_allNode)):
        road_link_id = int(item_number / len(origin))
        o = item_number % len(origin)
        for d in range(len(destination)):
            cost_disrupt[road_link_id][o][d] = to_allNode[item_number].get(destination[d]) # shortest path cost for each OD; # of link * OD *OD
    diff = cost_disrupt - baseCost
    Filedump(pd.DataFrame(cost_disrupt[1]),'cost_disrupt_%s' % name)
    # change cost of the isolated OD to penalty value
    iso = np.zeros((nLink, len(origin),len(destination)))
    for index,item in np.ndenumerate(cost_disrupt):
        if item >= 1e9:
            diff[index] = penalty
            iso[index] = 1
    return diff, iso

def GenerateDemandFunction(origin,destination,baseLength):
    """
    gravity models for travel demand.
    output:
    1.) output: a demand function of form  Trips[o,d] = k * Pop[o] * Pop[d] / Distance[o,d]
    """
    popu={}
    demand = np.zeros((len(origin['P']),len(destination['P'])))
    O_DF = gpd.read_file(origin['file'])
    D_DF = gpd.read_file(destination['file'])
    for o in O_DF.index:
        for d in D_DF.index:
            if (o == d) & (origin['name'] == destination['name']):
                pass
            else:
                if baseLength[o][d] < 1:
                    demand[o][d] = (O_DF[origin['scalar_column']][o] * D_DF[destination['scalar_column']][d])
                else:
                    demand[o][d] = ((O_DF[origin['scalar_column']][o] * D_DF[destination['scalar_column']][d]) / (baseLength[o][d]))
    Filedump(pd.DataFrame(demand),'demand_%s_%s' % (origin['name'], destination['name']))
    return demand

def summarise(diff, iso, demand, origin, destination):
    """
    For each link disrupted, we find the number of trips that cannot be completed due to link disruption: isolate_sumTrip
    Total social cost is defined as a link disrupted under traffic: disrupt_sumCost
    Total number of trips being disrupted: disrupt_sumTrip
    """
    disrupt_sumCost, disrupt_sumTrip, isolate_sumTrip = np.zeros(nLink),np.zeros(nLink),np.zeros(nLink)
    for i in range((nLink)):
        disrupt_sumCost[i] = np.nansum(np.multiply(diff[i,:,:],demand))
        disrupt_sumTrip[i] = np.nansum((diff[i,:,:]>0)*demand)
        isolate_sumTrip[i] = np.nansum(np.multiply(iso[i,:,:],demand))
    # criticality analysis with traffic for total user cost
    col_1 = gdf2['ID']
    col_2 = disrupt_sumCost
    col_3 = disrupt_sumTrip
    col_4 = isolate_sumTrip
    col_5 = np.array(col_2)/np.array(col_3)  # this is the average disuption cost for a link being disrupted
    criticality = np.column_stack((col_1,col_2,col_3, col_4,col_5))
    criticality[:,0] = criticality[:,0].astype(int)
    out = pd.DataFrame({'ID':criticality[:,0],'Social_Cost':criticality[:,1],'Disrupted_Trips':criticality[:,2],'Isolated_Trips':criticality[:,3],'averageCost':criticality[:,4]})
    a = out['Social_Cost']
    b = out['Disrupted_Trips']
    c = out['Isolated_Trips']
    out['Social_Cost_score'] = ((a - a.min()) / (a.max()- a.min()))
    out['Disrupted_score'] = ((b - b.min()) / (b.max()- b.min()))
    out['Isolated_score'] = ((c - c.min()) / (c.max()- c.min()))
    out['crit_score_%s_%s' % (origin['name'], destination['name'])] = (
        0.5 * out['Social_Cost_score'] +
        0.5 * out['Disrupted_score'] +
        0 * out['Isolated_score'])
    out = out[['ID','crit_score_%s_%s' % (origin['name'], destination['name'])]]
    return out

def Main(origin, destination):
    Vprint('\ncomputing for origin = %s and destination = %s\n' % (origin['name'], destination['name']))
    origin['P'], destination['P'] = PrepSet(origin), PrepSet(destination)
    basePath, baseCost, baseLength = Pathfinder(origin['P'], destination['P'])
    diff, iso = BreakEdge(origin['P'], destination['P'], destination['penalty'], baseCost, destination['name'])
    demand = GenerateDemandFunction(origin, destination, baseLength)
    summary = summarise(diff, iso, demand, origin, destination)
    Outputs.append({
        'origin': origin['name'],
        'destination': destination['name'],
        'basePath': basePath,
        'baseCost': baseCost,
        'baseLength': baseLength,
        'diff': diff,
        'iso': iso,
        'summary': summary})

Main(origin_points, origin_points)
Main(origin_points, hospitals)

Output = inNetwork.drop(['geometry','TC_iri_med','total_cost'],axis =1)
for o_d_calc in range(0,len(Outputs)):
    Output = Output.merge(Outputs[o_d_calc]['summary'],how = 'left', on = 'ID')
FileOut(Output,'criticality_output')
