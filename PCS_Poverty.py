# -*- coding: utf-8 -*-
import pandas as pd
import geopandas as gpd
import shapely.wkt
import os
#Global vars
district = str(raw_input('\nDistrict Code: (YD | TT) '))
path = r'C:\Users\charl\Documents\GitHub\RoadLabPro_Utils\\'
povpath = os.path.join(path,'PCS','Poverty','Poverty_Communes_2009.shp')
roadpath = os.path.join(path,'Runtime', '%s' % district, 'Network.csv')
dash = os.path.join(path, 'PCS',r'dashboard.xlsx')
crs = {'init': 'epsg:4326'}
admin = r'Poverty_Communes_2009.shp'
#read CMS output (roads) in as GDF
df = pd.read_csv(roadpath)
geometry = df['Line_Geometry'].map(shapely.wkt.loads)
gdf_road = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
#Read in Shapefile as GDF
gdf_adm = gpd.read_file(povpath)
gdf_adm = gdf_adm.to_crs(crs)
#read in weights as dictionary
weightsdf = pd.read_excel(dash, sheetname = "POV")
w8s = weightsdf.to_dict(orient='records')
for x in range(0,10):
    print "Doing calcs for: %s" % w8s[x]['B']
    plc = ((gdf_adm[w8s[x]['B']]-gdf_adm[w8s[x]['B']].min())/(gdf_adm[w8s[x]['B']].max()-gdf_adm[w8s[x]['B']].min()))
    if w8s[x]['D'] == False:
        gdf_adm['Povcomp_%d'% (x+1)] = (1 - plc)*float(w8s[x]['C'])
    else:
        gdf_adm['Povcomp_%d'% (x+1)] = plc*float(w8s[x]['C'])
gdf_adm['POV_SCORE'] = (
    gdf_adm['Povcomp_1'] +
    gdf_adm['Povcomp_2'] +
    gdf_adm['Povcomp_3'] +
    gdf_adm['Povcomp_4'] +
    gdf_adm['Povcomp_5'] +
    gdf_adm['Povcomp_6'] +
    gdf_adm['Povcomp_7'] +
    gdf_adm['Povcomp_8'] +
    gdf_adm['Povcomp_9'] +
    gdf_adm['Povcomp_10']
    )
gdf_adm = gdf_adm[['CCode04New','Ccode02_05','DISTCODE02','P_EName','D_EName','DISTCODE04','__Commun_1','COMNAME','COMPOPULA','__Commune','POV_SCORE','geometry']]
gdf_adm.to_file(os.path.join(path, 'PCS','Poverty','pov_layer.shp'), driver = 'ESRI Shapefile')
gdf_adm_join = gdf_adm[['geometry','POV_SCORE']]
gdf_adm = gdf_adm.drop('geometry', axis=1)
gdf_adm.to_excel(povpath+'Pov_all_communes.xlsx')

#Spatial Join
def LINEGROUPER(x2):    #EDIT
    y = pd.DataFrame()
    y['POV_SCORE'] = [x2.POV_SCORE.mean()]
    y['ID'] = [x2.ID.iloc[0]]
    return y

#create mini DF with just the average poverty score and VPROMMS_ID
gdf_road_join = gdf_road[['geometry','ID']]
join_pov = gpd.sjoin(gdf_road_join, gdf_adm_join, how = "inner", op = 'intersects')
join_pov = join_pov.groupby(['ID']).apply(lambda x: LINEGROUPER(x))
gdf_road = gdf_road.merge(join_pov, how = "inner", on = "ID")
gdf_road = gdf_road.drop('geometry',axis = 1)
gdf_road = pd.DataFrame(gdf_road)
gdf_road.to_csv(os.path.join(path,'Outputs','%s' % district,'poverty_output.csv'), index = False)
