import pandas as pd
import numpy as np
import glob2
import re
import os
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns
scrubber = raw_input('\nRemove all non-VPROMMS_ID roads? y = yes, n = no\n')
district = str(raw_input('\nDistrict Code: (YD | TT) '))
InPath = r'C:\Users\charl\Documents\Vietnam\Fieldwork\\roadlab_bin_%s' % district
OutPath = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_%s' % district
sdthresh = 0.5

# Line file
file_names_list = glob2.glob(InPath+'\\'+'**'+'\\*intervals*.csv')

for files in file_names_list:
    print files+'\n'
print len(file_names_list)

dataframes=[]
for files in file_names_list:
    df = pd.DataFrame(pd.read_csv(files))
    df['Input_file'] = str(files)
    dataframes.append(df)
X = pd.concat(dataframes, ignore_index=True)
X['Line_Geometry'] = 'LINESTRING ('+X['start_lon'].map(str)+' '+X['start_lat'].map(str)+', '+X['end_lon'].map(str)+' '+X['end_lat'].map(str)+')'
X['Part1'] = X['Input_file'].map(str).str.extract(('(bin.*\.csv)'), expand=False)
X['Part1'] = X['Part1'].map(str).str.replace('bin', '').str.replace('\.csv', '')
X['Part1'] = X['Part1'].str.split('\\').str.get(-3)
X['Part2'] = X['Input_file'].map(str).str.extract(('(bin.*\.csv)'),expand=False)
X['Part2'] = X['Part2'].map(str).str.replace('bin', '').str.replace('\.csv', '')
X['Part2'] = X['Part2'].str.split('\\').str.get(-2)
for y in X.Part1.unique():
    miniframe = X.loc[X['Part1'] == y]
    n = 1
    for y2 in miniframe.Part2.unique():
        X.loc[(X['Part2'] == y2)&(X['Part1']==y),'Part1'] = X.loc[(X['Part2'] == y2)&(X['Part1']==y),'Part1']+'_'+str(n)
        n +=1
X['VPROMMS_ID'] = X['Part1']
X = X.drop(['Part2','Part1','Input_file'], axis = 1)
if scrubber == 'y':
    X = X[X.VPROMMS_ID.str.contains('OTHER*') == False]
else:
   pass
X.to_csv(OutPath+'\\'+'FOR_JOIN_INTS.csv')

# Point file
file_names_list2 = glob2.glob(InPath+'\\'+'**'+'\\*RoadPath'+'*.csv')

for files in file_names_list2:
    print files+'\n'
print len(file_names_list2)

dataframes2=[]
for files in file_names_list2:
    df = pd.DataFrame(pd.read_csv(files))
    df['Input_file'] = str(files)
    dataframes2.append(df)
Y = pd.concat(dataframes2, ignore_index=True)
Y['Point_Geometry'] = 'POINT ('+Y['longitude'].map(str)+' '+Y['latitude'].map(str)+')'
Y['Part1'] = Y['Input_file'].map(str).str.extract(('(bin.*\.csv)'),expand=False)
Y['Part1'] = Y['Part1'].map(str).str.replace('bin', '').str.replace('\.csv', '')
Y['Part1'] = Y['Part1'].str.split('\\').str.get(-3)
Y['Part2'] = Y['Input_file'].map(str).str.extract(('(bin.*\.csv)'),expand=False)
Y['Part2'] = Y['Part2'].map(str).str.replace('bin', '').str.replace('\.csv', '')
Y['Part2'] = Y['Part2'].str.split('\\').str.get(-2)
for y in Y.Part1.unique():
    miniframe = Y.loc[Y['Part1'] == y]
    n = 1
    for y2 in miniframe.Part2.unique():
        Y.loc[(Y['Part2'] == y2)&(Y['Part1']==y),'Part1'] = Y.loc[(Y['Part2'] == y2)&(Y['Part1']==y),'Part1']+'_'+str(n)
        n +=1
Y['VPROMMS_ID'] = Y['Part1']
Y = Y.drop(['Part2','Part1','Input_file'], axis = 1)
if scrubber == 'y':
    Y = Y[Y.VPROMMS_ID.str.contains('OTHER*') == False]
else:
    pass

def main(X):
    timeconvert(X)
    deltas(X)
    X['disterror'] = outliers(X['distdiff'].as_matrix(), 'distance')
    X['timeerror'] = outliers(X['timediff'].as_matrix(), 'time')
    X['error'] = ((X['timeerror'] == True) | (X['disterror'] == True))
    Z = X.drop(['latdiff','longdiff'], axis = 1)
    Z.to_csv(os.path.join(OutPath,'FOR_JOIN_PATHS.csv'))
    #plt.figure()
    #plotter(Z['distdiff'].dropna(), name = 'dist_deltas')
    #plt.figure()
    #plotter(Z['timediff'].dropna().as_matrix(), name = 'time_deltas')

# Convert RLP timestamp to pandas time object
def timeconvert(X):
    X['time'] = pd.to_datetime(X['time'],infer_datetime_format = True)
    try:
        X = X.drop(['Unnamed: 0'], axis =1)
    except:
        pass
    return X

def deltas(X):
    for y in X.VPROMMS_ID.unique():
        MF = X.loc[X['VPROMMS_ID'] == y]
        X.loc[X['VPROMMS_ID'] == y, 'timediff'] = MF.time.diff() / np.timedelta64(1, 's')
        X.loc[X['VPROMMS_ID'] == y, 'latdiff'] = MF.latitude.diff()
        X.loc[X['VPROMMS_ID'] == y, 'longdiff'] = MF.longitude.diff()
    X['distdiff'] = np.sqrt(np.square(X['latdiff'])+np.square(X['longdiff']))

def outliers(points,name,thresh = sdthresh):
    print '\n---Summary of outlier detection: %s---' % name
    mean = np.nanmean(points, axis=0)
    print '\nmean =',mean
    squarediff = (points - mean)**2
    print '\nsquared differences = ',squarediff
    sumsquarediff = np.nansum(squarediff)
    print '\nsumsquarediff = ',sumsquarediff
    N = len(points)
    print '\nN = ',N
    sd_deviation = np.sqrt((sumsquarediff / N))
    print '\nsd_deviation = ',sd_deviation
    z_score = (points - mean) / sd_deviation
    print z_score
    return z_score > thresh

def plotter(Z, name):
    ax = sns.distplot(Z, bins = 20, hist = False, rug=True, label = name)
    fig1 = ax.get_figure()
    fig1.savefig("%s\\%s.png" % (OutPath, name))

main(Y)
