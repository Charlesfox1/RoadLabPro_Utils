# coding: utf-8
import os, sys
import shapely
import rtree
import pandas as pd
import geopandas as gdp
district = str(raw_input('\nDistrict Code: (YD | TT) '))
path = r'C:\Users\charl\Documents\GitHub\RoadLabPro_Utils\\'
Outpath = os.path.join(path,'Outputs','%s' % district)
roadpath = os.path.join(path, 'runtime', '%s' % district,"Network.csv")
dash = os.path.join(path,'PCS','dashboard.xlsx')
weightsdf = pd.read_excel(dash, sheetname = "DASH_MIRROR", index_col = 0)
df = pd.read_csv(roadpath)

#Define component values
#Component1: Poverty
#Poverty - Vietnam poverty map. Tati has his name. African. (simple extract attributes)
poverty = pd.read_csv(os.path.join(Outpath,'Poverty_output.csv'))
poverty = poverty[['ID','POV_SCORE']]

#Component2: Disaster Risk
risk = pd.read_csv(os.path.join(Outpath,'risk_output.csv'))
risk = risk[['ID','RISK_SCORE']]

#Component3: Criticality
criticality = pd.read_csv(os.path.join(Outpath,'criticality_output.csv'))
criticality = criticality[['ID','CRIT_SCORE']]

#Component4: Access Index
#access index - 'keith megatool' KG
component4 = 1
df['component4'] = 1

#Component5: Roughness
iridf = pd.read_excel(dash, sheetname = "ROUGHNESS", index_col = 0)
roughness = pd.DataFrame({'ID':df['ID']})
roughness['a'] = (
        iridf['WEIGHT']['iri_med']*df['iri_med'] +
        iridf['WEIGHT']['iri_stdev']*df['iri_stdev'] +
        iridf['WEIGHT']['iri_min']*df['iri_min'] +
        iridf['WEIGHT']['iri_max']*df['iri_max'] +
        iridf['WEIGHT']['iri_mean']*df['iri_mean'])
roughness['ROUGHNESS_SCORE'] = ((roughness['a'] - roughness['a'].min()) / (roughness['a'].max() - roughness['a'].min()))
roughness = roughness.drop(['a'],axis = 1)

#Calculate PCS
for component in [poverty, risk, criticality, roughness]:
    df = df.merge(component, how = 'left', on = 'ID')

df["Pov_weight"] = weightsdf['Weight']['PCS_POV']
df["Risk_weight"] = weightsdf['Weight']['PCS_RISK']
df["Crit_weight"] = weightsdf['Weight']['PCS_CRIT']
df["Access_weight"] = weightsdf['Weight']['PCS_ACCESS']
df["Rough_weight"] = weightsdf['Weight']['PCS_ROUGH']
df["PCS"] = (
    df['POV_SCORE']*df["Pov_weight"] +
    df['RISK_SCORE']*df["Risk_weight"] +
    df['CRIT_SCORE']*df["Crit_weight"] +
    df['component4']*df["Access_weight"] +
    df['ROUGHNESS_SCORE']*df["Rough_weight"]
)*100

#df cleanup
df.to_csv(os.path.join(Outpath,"PCS.csv"))
