# -*- coding: utf-8 -*-
###################################################################################################
# Calculate a poverty score for various linear features
# Charles Fox (and a little by Ben), September 2017
# Purpose: determine the poverty index for each linear feature in the network dataset
###################################################################################################
import os, sys, inspect, logging

import pandas as pd
import geopandas as gpd
import shapely.wkt

roadID = ''

def main(district="YD", admin="Poverty_Communes_2009.shp", curRoadID="ID"):
    #district = str(raw_input('\nDistrict Code: (YD | TT) '))
    #admin = r'Poverty_Communes_2009.shp'
    roadID = curRoadID
    path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    
    #Set up logging
    logging.basicConfig(filename = os.path.join(path, "PCS_Poverty.log"), level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")    
    logging.info("Starting PCS Poverty Process")
    
    povpath = os.path.join(path,'PCS','Poverty')
    povAdmin = os.path.join(povpath, admin)
    roadpath = os.path.join(path,'Runtime', '%s' % district, 'Network.csv')
    dash = os.path.join(path, 'PCS',r'dashboard.xlsm')
    crs = {'init': 'epsg:4326'}
    
    for curFile in [povpath, povAdmin, roadpath, dash]:
        if not os.path.exists(curFile):
            logging.error("No input found: %s" % curFile)
            raise ValueError("No input found: %s" % curFile)
            
    #read CMS output (roads) in as GDF
    df = pd.read_csv(roadpath)
    geometry = df['Line_Geometry'].map(shapely.wkt.loads)
    gdf_road = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    #Read in Shapefile as GDF
    gdf_adm = gpd.read_file(povAdmin)
    gdf_adm = gdf_adm.to_crs(crs)
    #read in weights as dictionary
    weightsdf = pd.read_excel(dash, sheetname = "POV")
    w8s = weightsdf.to_dict(orient='records')
    
    for x in range(0,len(w8s)):
        logging.debug("Doing calcs for: %s" % w8s[x]['B'])
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
    
    #gdf_adm = gdf_adm[['CCode04New','Ccode02_05','DISTCODE02','P_EName','D_EName','DISTCODE04','__Commun_1','COMNAME','COMPOPULA','__Commune','POV_SCORE','geometry']]
    gdf_adm.to_file(os.path.join(path, 'PCS','Poverty','pov_layer.shp'), driver = 'ESRI Shapefile')
    gdf_adm_join = gdf_adm[['geometry','POV_SCORE']]
    gdf_adm = gdf_adm.drop('geometry', axis=1)
    gdf_adm.to_excel(os.path.join(povpath, 'Pov_all_communes.xlsx'))
    
    #create mini DF with just the average poverty score and VPROMMS_ID        
    gdf_road_join = gdf_road[['geometry',roadID]]    
    join_pov = gpd.sjoin(gdf_road_join, gdf_adm_join, how = "inner", op = 'intersects')    
    join_pov = join_pov.groupby([roadID]).apply(lambda x: LINEGROUPER(x))    
    gdf_road = gdf_road.merge(join_pov, how = "inner", on = roadID)    
    gdf_road = gdf_road.drop('geometry',axis = 1)
    gdf_road = pd.DataFrame(gdf_road)
    gdf_road['POV_SCORE'] = ((gdf_road['POV_SCORE'] - gdf_road['POV_SCORE'].min()) / (gdf_road['POV_SCORE'].max() - gdf_road['POV_SCORE'].min()))
    gdf_road.to_csv(os.path.join(path,'Outputs','%s' % district,'poverty_output.csv'), index = False)
    logging.info("Finished PCS Poverty Calculations")
    
#Spatial Join
def LINEGROUPER(x2):    #EDIT
    y = pd.DataFrame()
    y['POV_SCORE'] = [x2.POV_SCORE.mean()]
    y['ID'] = [x2.ID.iloc[0]]
    return y

main()