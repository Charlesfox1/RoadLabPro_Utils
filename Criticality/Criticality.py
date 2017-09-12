#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 15:54:18 2017

@author: yan
"""
# criticality analysis
#
#get_ipython().magic(u'matplotlib inline')
from matplotlib import pyplot as plt
import networkx as nx
import pandas as pd
import copy
import time
import csv
import os
import sys

from collections import Counter

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

#Modules developed by TU Delft team
from network_lib import network_prep as net_p
from network_lib import network_visualization as net_v
from network_lib import od_prep as od_p
from network_lib import weighted_betweenness as betw_w

import geopandas as gp
import numpy as np
from simpledbf import Dbf5
def Vprint(a):
    if verbose == 1:
        print a
    else:
        pass
verbose = 1
# Network Preparation
outpath = r'C:\Users\charl\Documents\Vietnam\PCS'
network = r'./input/Vietnam/corrected.shp'
centroid = r'./input/Vietnam/YD_centroids.shp'
#dbf = Dbf5(r'./input/Vietnam/fj_bridges_south_FRA_v1.dbf')
#df_structure = dbf.to_dataframe()
gdf_points, gdf_node_pos, gdf = net_p.prepare_centroids_network(centroid, network)

# Create Networkx MultiGraph object from the GeoDataFrame
G = net_p.gdf_to_simplified_multidigraph(gdf_node_pos, gdf, simplify=False)

# Change the MultiGraph object to Graph object to reduce computation cost
G_tograph = net_p.multigraph_to_graph(G)

# Observe the properties of the Graph object
Vprint('number of disconnected compoents is: %d' % nx.number_connected_components(G_tograph))
nx.info(G_tograph)

# Take only the largest subgraph which all connected links
len_old = 0
for g in nx.connected_component_subgraphs(G_tograph):
    if len(list(g.edges())) > len_old:
        G1 = g
        len_old = len(list(g.edges()))
G_sub = G1.copy()

Vprint('number of disconnected compoents is: %d' % nx.number_connected_components(G_sub))
nx.info(G_sub)

# Save the simplified transport network back into GeoDataFrame
gdf_sub = net_p.graph_to_df(G_sub)
gdf_sub.to_csv(r'./runtime_working_files/RoadGraph.csv')

# assign the OD to the closest node of the biggest subgraph:
gdf_points2, gdf_node_pos2, gdf_new=net_p.prepare_newOD(centroid, gdf_sub)
G2_multi = net_p.gdf_to_simplified_multidigraph(gdf_node_pos2, gdf_new, simplify=False)
G2 = net_p.multigraph_to_graph(G2_multi)
gdf2 = net_p.graph_to_df(G2)
allNode = G2.nodes()
allEdge = G2.edges()
od = gdf_points2['Node']
outlist = [od, gdf2]
gdf_points2.to_csv(r'./runtime_working_files/gdf_points2.csv')
od.to_csv(r'./runtime_working_files/od.csv')
gdf2.to_csv(r'./runtime_working_files/gdf2.csv')
gdf_node_pos2.to_csv(r'./runtime_working_files/gdf_node_pos2.csv')
# the output of this section is gdf_points2, gdf_node_pos2, gdf2, G2

# find the shortest path for each OD to minimize the total travel cost;
# output:
# 1) baseCost ($): total travel cost between all OD pairs;
# 2) basePath the shortest path between all OD pairs

n=0
basePath = [[[]for i in range(len(od))] for j in range(len(od))]
baseCost=np.zeros((len(od),len(od)))
for i in range(len(od)):
    for j in range(i+1,len(od)):
        Vprint(n)
        basePath[i][j]=nx.dijkstra_path(G2,od[i],od[j],weight = 'total_cost')
        baseCost[i][j]=nx.dijkstra_path_length(G2,od[i],od[j],weight = 'total_cost')
        n=n+1

Vprint( "basePath: ")
Vprint( basePath)
Vprint( "baseCost: ")
Vprint( baseCost)
# find the shortest path for each od to minimize the total travel cost;
# output: baseLength: total travel distance between all OD pairs;
n=0
baseLength=np.zeros((len(od),len(od)))
for i in range(len(od)):
    for j in range(i+1,len(od)):
        Vprint(n)
        baseLength[i][j]=nx.dijkstra_path_length(G2,od[i],od[j],weight = 'length')
        n=n+1
Vprint('Baselengths:' )
Vprint( baseLength)
# criticaility analysis: remove each link by changing it to a extremely large value, calculate the shortest path
# Here, if the criticaility is for only a subset, instead of all link: change the for loop in range(nLink) to the
# subset you would like to change, and "G[gdf2['FNODE_'][n]][gdf2['TNODE_'][n]]['total_cost']=1e10  " is used
# to change the weight (travel time or travel cost) of the graph.
nLink = len(allEdge)
Vprint('allEdge variable:')
Vprint(allEdge)
to_allNode = []
G = copy.deepcopy(G2)
for n in range(len(gdf2)):
    start = time.clock()
    Vprint('\nRunning calculations for road n (%d) of %d (range len(gdf2)):' % (n+1, len(gdf2)))
    w = G[gdf2['FNODE_'][n]][gdf2['TNODE_'][n]]['total_cost']
    Vprint("\n'w', or original cost of traversing, is %d" % w)
    Vprint('\nMaking this segment, defined by a starting node of "%s", and an end node of "%s", too expensive to traverse' % (gdf2['FNODE_'][n], gdf2['TNODE_'][n]))
    G[gdf2['FNODE_'][n]][gdf2['TNODE_'][n]]['total_cost']=1e10
    Vprint('Iterating through all OD destinations')
    for i in range(len(od)):
        to_allNode.append(nx.single_source_dijkstra_path_length(G,od[i],weight = 'total_cost'))
        Vprint('\nCalculate cost of moving from node number (%d) [node code "%d"], of %d total nodes, to all other nodes:' % (i+1, od[i], len(od)))
        Vprint('cost from %d to other nodes in network is: %s ' % (od[i], nx.single_source_dijkstra_path_length(G,od[i],weight = 'total_cost')))
    G[gdf2['FNODE_'][n]][gdf2['TNODE_'][n]]['total_cost']= w # Resetting traverse cost to original (w)
    Vprint('%s seconds' % (time.clock() - start))

# the disruption cost
penalty = 1000 # the penalty is the dollar value assigned as a penalty when a trip is isolated
cost_disrupt= np.zeros((nLink, len(od),len(od)))
for i in range(len(to_allNode)):
    id_link = int(i/len(od))
    Vprint('id_link, or ID of the road segment: %s' % id_link)
    id_o = i%len(od)
    Vprint("id_o, or ID of the origin: %s" % id_o)
    for j in range(len(od)):
        if j>id_o:
            cost_disrupt[id_link ][id_o][j] = to_allNode[i].get(od[j]) # shortest path cost for each OD; # of link * OD *OD
diff = np.zeros((nLink, len(od),len(od)))
base = np.zeros((nLink, len(od),len(od)))
iso = np.zeros((nLink, len(od),len(od)))
for i in range(len(gdf2)):
    base[i,:,:] = baseCost
diff = cost_disrupt - base

# change cost of the isolated OD to penalty value
for index,item in np.ndenumerate(cost_disrupt):
    if item>=1e9:
        diff[index]=penalty
        iso[index]=1

# gravity model for the travel demand information: Tij=k*pi*pj/dij; travel demand per day (AADT)
k=1.3
popu={}
for i in range(len(od)): # puts the population of the cities into a dict with node code as key
    idx = gdf_points2['Node'].tolist().index(od[i])
    popu[od[i]]=gdf_points2['Pop'].tolist()[idx]

T=np.zeros((len(od),len(od)))
for i in range(len(od)):
    for j in range(i+1,len(od)):
        T[i][j]= k*popu[od[i]]*popu[od[j]]/baseLength[i][j]/30*1e3
#constant X population of i X population of J / length i to j / constant
#
# For each link disrupted, we find the number of trips that cannot be completed due to link disruption: isolate_sumTrip
# Total social cost is defined as a link disrupted under traffic: disrupt_sumCost
# Total number of trips being disrupted: disrupt_sumTrip
disrupt_sumCost=np.zeros(nLink)
disrupt_sumTrip=np.zeros((nLink))
isolate_sumTrip = np.zeros((nLink))
for i in range((nLink)):
    disrupt_sumCost[i] = np.sum(np.multiply(diff[i,:,:],T))
    disrupt_sumTrip[i] = np.sum((diff[i,:,:]>0)*T)
    isolate_sumTrip[i] = np.sum(np.multiply(iso[i,:,:],T))

for i in range((nLink)):
    print'\nFor the segment defined by a starting node of "%s", and an end node of "%s"' % (gdf2['FNODE_'][i], gdf2['TNODE_'][i])
    print "Total number of cancelled trips:", isolate_sumTrip[i]
    print "Total number of disrupted journeys:", disrupt_sumTrip[i]
    print "Total social cost of disruption of link i:", disrupt_sumCost[i]

# criticality analysis with traffic for total user cost
col_1 = gdf2['id']
col_2 = disrupt_sumCost
col_3 = disrupt_sumTrip
col_4 = isolate_sumTrip
col_5=np.array(col_2)/np.array(col_3)  # this is the average disuption cost for a link being disrupted

criticality=np.column_stack((col_1,col_2,col_3, col_4,col_5))
criticality[:,0]=criticality[:,0].astype(int)
criticality_traffic = pd.DataFrame({'Road_id':criticality[:,0],'disrupt_sumCost':criticality[:,1],'disrupt_sumTrip':criticality[:,2],'isolated_sumTrip':criticality[:,3],'averageCost':criticality[:,4]})
criticality_traffic.to_csv(outpath+'\\criticality_output.csv', sep=',', index=False)
