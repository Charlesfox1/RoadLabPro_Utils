# -*- coding: utf-8 -*-
###################################################################################################
# Calculate a roughness score for various linear features
# Charles Fox (and a little by Ben), September 2017
# Purpose: determine the roughness score for each linear feature in the network dataset
###################################################################################################

import os, sys, inspect, logging
import pandas as pd

def Main(district="test"):
    district = str(raw_input('\nDistrict Code: (YD | TT) '))
    path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))

    #Set up logging
    logging.basicConfig(filename = os.path.join(path, 'runtime', district, "PCS_Roughness_log.log"), level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
    logging.info("Starting PCS Roughness Process")

    #Check for existance of input files
    roadpath = os.path.join(path, 'runtime', '%s' % district,'Network.csv')
    dash = os.path.join(path,'PCS','dashboard.xlsm')
    for curFile in [dash, roadpath]:
        if not os.path.exists(curFile):
            logging.error("No input found: %s" % curFile)
            raise ValueError("No input found: %s" % curFile)

    #read in input files
    iridf = pd.read_excel(dash, sheetname = 'ROUGHNESS', index_col = 0)
    roughness = pd.read_csv(roadpath)

    #check to see if there are at least 2 roads in roads df
    if roughness.index.max() < 2:
        logging.error("roadfile contains fewer than 2 roads")
        raise ValueError("roadfile contains fewer than 2 roads")

    #Calculate Roughness Score
    logging.debug("Calculating Roughness Score")
    roughness['ROUGHNESS_SCORE'] = (
            iridf['WEIGHT']['iri_med']*roughness['iri_med'] +
            iridf['WEIGHT']['iri_min']*roughness['iri_min'] +
            iridf['WEIGHT']['iri_max']*roughness['iri_max'] +
            iridf['WEIGHT']['iri_mean']*roughness['iri_mean'])

    #Normalise Roughness Score
    logging.debug("Normalising Roughness Score")
    roughness['ROUGHNESS_SCORE'] = ((roughness['ROUGHNESS_SCORE'] - roughness['ROUGHNESS_SCORE'].min()) / (roughness['ROUGHNESS_SCORE'].max() - roughness['ROUGHNESS_SCORE'].min()))
    outpath = os.path.join(path, 'Outputs','%s' % district)
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    try:
        roughness.to_csv(os.path.join(outpath,'roughness_output.csv'), index = False)
    except:
        logging.error("error writing file to: %s" % outpath)
        raise ValueError("error writing file to: %s" % outpath)
    logging.info("Finished PCS Roughness Calculations, without issue.")
Main()
