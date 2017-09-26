# -*- coding: utf-8 -*-
"""
Created on Wed 20 Sep
@author: Charles Fox
"""
import pandas as pd
import numpy as np
import os

path = r'C:\Users\charl\Documents\Vietnam\Ronet'
f = 'RONET.xlsx'
tables = os.path.join(path, f)

# Construct DFs
print 'reading in'
Roads = pd.read_excel(tables, sheetname = 'Roads', header = 0, index_col = 0)
Control = pd.read_excel(tables, sheetname = 'Control', header = 0, index_col = 0)
VehicleFleet = pd.read_excel(tables, sheetname = 'VehicleFleet', header = 0, index_col = 0)
TrafficGrowth = pd.read_excel(tables, sheetname = 'TrafficGrowth', header = 0, index_col = 0)
RoadWorks = pd.read_excel(tables, sheetname = 'RoadWorks', header = [0,1], index_col = 0)
RecurrentMaintenance = pd.read_excel(tables, sheetname = 'RecurrentMaintenance', header = [0,1], index_col = 0)
Width = pd.read_excel(tables, sheetname = 'Width', header = 0, index_col = 0)
MCoeff = pd.read_excel(tables, sheetname = 'MCoeff', header = 0, index_col = 0)
SurfaceType = pd.read_excel(tables, sheetname = 'SurfaceType', header = 0, index_col = 0).squeeze()
ConditionData = pd.read_excel(tables, sheetname = 'ConditionData', header = 0, index_col = [0,1])
TrafficLevels = pd.read_excel(tables, sheetname = 'TrafficLevels', header = [0,1], index_col = 0)
VOC = pd.read_excel(tables, sheetname = 'VOC', header = 0, index_col = 0)
Speed = pd.read_excel(tables, sheetname = 'Speed', header = 0, index_col = [0,1,2])
RoadDeterioration = pd.read_excel(tables, sheetname = 'RoadDeterioration', header = 0, index_col = [0,1])
WorkEvaluated = pd.read_excel(tables, sheetname = 'WorkEvaluated', header = [0,1], index_col = [0,1,2])
WorkEvaluated.index.set_names(['SurfaceType','RoadClass','ConditionClass'],inplace = True)
VehicleTypes = ['Motor_cycle','Small_Car','Medium_Car','Delivery_Vehicle','4_Wheel_Drive','Light_Truck','Medium_Truck','Heavy_Truck','Articulated_Truck','Small_Bus','Medium_Bus','Large_Bus']

# Add defaults
Roads['SurfaceType'] = Roads['SurfaceType'].fillna(Roads.RoadType.replace(SurfaceType))
Roads['Width'] = Roads['Width'].fillna(Roads.Lanes.replace(Width['Default Carriageway Width (m)']))
Roads = Roads.set_index(['SurfaceType','ConditionClass']).fillna(ConditionData).reset_index()
foo = pd.DataFrame(TrafficLevels['Pavement Condition Class'].stack(), columns = ['StructuralNo'])
foo.index.set_names(['TrafficLevel','ConditionClass'],inplace = True)
Roads.loc[Roads.SurfaceType < 4] = Roads.set_index(['TrafficLevel','ConditionClass']).fillna(foo).reset_index()
bah = []
for vehicle in VehicleTypes:
    Roads['AADT_%s' % vehicle] = Roads['AADT_%s' % vehicle].fillna(Roads.TrafficLevel.replace((TrafficLevels['Default AADT (%)','%s' % vehicle]) * TrafficLevels['Default Total Traffic','(veh/day)']))
    bah.append(Roads['AADT_%s' % vehicle])
Roads['AADT_Total'] = sum(bah)
RoadsDict = Roads.to_dict(orient = 'index')

# Define Alternatives
def AlternativeFrame(road):
    a = WorkEvaluated['BaseCase','RoadWorkNumber'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    b = WorkEvaluated['AlternativesRoadWorks','Alternative1'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    c = WorkEvaluated['AlternativesRoadWorks','Alternative2'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    Work = pd.Series([a,b,b,b,b,b,b,c,c,c,c,c,c], index=[range(0,13)])
    Year = pd.Series(np.nan, index=[range(0,13)])
    Alternatives = pd.DataFrame.from_dict({'Work':Work,'Year':Year}, orient = 'columns')
    Alternatives.loc[0] = WorkEvaluated['BaseCase','RoadWorkYear'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    if pd.notnull(Alternatives['Work'].loc[1]):
        Alternatives['Year'].loc[1:6] = range(1,7)
    if pd.notnull(Alternatives['Work'].loc[7]):
        Alternatives['Year'].loc[7:13] = range(1,7)
    Alternatives.index.set_names('Alternative')
    if Alternatives['Work'].count() == 1:
        Alternatives.loc[1] = Alternatives.loc[0]
    iNoAlternatives = Alternatives['Work'].count()
    return Alternatives

# Main Call
for x in Roads.index:
    road = RoadsDict[x]
    Alternatives = AlternativeFrame(road)



Roads.to_csv(os.path.join(path,'Output_Road.csv'))
