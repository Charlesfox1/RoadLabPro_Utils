import os, sys
import pandas as pd
path = r'C:\Users\charl\Documents\GitHub\RoadLabPro_Utils\\'
district = str(raw_input('\nDistrict Code: (YD | TT) '))
roadpath = os.path.join(path, 'runtime', '%s' % district,'Network.csv')
df = pd.read_csv(roadpath)
dash = os.path.join(path,'PCS','dashboard.xlsx')
iridf = pd.read_excel(dash, sheetname = 'ROUGHNESS', index_col = 0)

roughness = pd.read_csv(roadpath)
roughness['ROUGHNESS_SCORE'] = (
        iridf['WEIGHT']['iri_med']*roughness['iri_med'] +
        iridf['WEIGHT']['iri_min']*roughness['iri_min'] +
        iridf['WEIGHT']['iri_max']*roughness['iri_max'] +
        iridf['WEIGHT']['iri_mean']*roughness['iri_mean'])
roughness['ROUGHNESS_SCORE'] = ((roughness['ROUGHNESS_SCORE'] - roughness['ROUGHNESS_SCORE'].min()) / (roughness['ROUGHNESS_SCORE'].max() - roughness['ROUGHNESS_SCORE'].min()))
roughness.to_csv(os.path.join(path,'Outputs','%s' % district,'roughness_output.csv'), index = False)
