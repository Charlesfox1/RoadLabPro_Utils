import pandas as pd
import os
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns
import sys

points = 'IRIpoints_csv_format_n.csv'
lines = 'IRIlines_csv_format_n.csv'
district = str(raw_input('\nDistrict Code: (YD | TT)'))
Path = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_%s' % district

# File Management
try:
    os.mkdir(Path+'\\AnalysisOutputs',0755)
    os.mkdir(Path+'\\AnalysisOutputs'+'\\SpeedAnalysis',0755)
    os.mkdir(Path+'\\AnalysisOutputs'+'\\LongitudinalPlots',0755)
except:
    pass
OutPath = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_%s\\AnalysisOutputs' % district
OutPathSpeeds = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_%s\\AnalysisOutputs\\SpeedAnalysis' % district
OutPathLongs = r'C:\Users\charl\Documents\Vietnam\Analysis\Workflow_%s\\AnalysisOutputs\\LongitudinalPlots' % district
pointsdf = pd.read_csv(os.path.join(Path, points))
pnttotal = pointsdf['speed'].count()
linesdf = pd.read_csv(os.path.join(Path, lines))
lntotal = linesdf['speed_med'].count()
IDs = pointsdf['VPROMMS_ID'].unique()
IDLIST = IDs.tolist() #['212TH00008_1','212TH00029_1','212TH00047_1','212TH00047_2']
answer = int(raw_input('RelativeIRI or Absolute for Longitudinal plots? [1,2]: '))
if answer == 1:
    z = 'RelativeIRI'
elif answer == 2:
    z = 'iri'
else:
    sys.exit
def statisticalreportpoints(x,y):
    a = x[['iri','speed']].describe()
    b = y[['length','iri_med','speed_med']].describe()
    writer = pd.ExcelWriter(OutPath+'\\statoutput_%s.xlsx' % district, engine='xlsxwriter')
    a.to_excel(writer,'Points')
    b.to_excel(writer,'Lines')
    writer.save()

def plotterhist(Z,z,z1,R,r,x):
    sns.distplot(Z, label = z, bins = 10, color = z1)
    sns.distplot(R, label = r, bins = 10)
    plt.legend()
    plt.savefig("%s\\%s_%s_hist_%s.png" % (OutPathSpeeds, z, r,x))
    plt.clf()

def speedonecondition(a,b,d):
    if a == 'lessthan':
        c = (pointsdf['speed']<=b)
    elif a == 'morethan':
        c = (pointsdf['speed']>=b)
    else:
        print 'error: choose lessthan or morethan'
    df = pointsdf.loc[c]
    df.to_csv(os.path.join(OutPathSpeeds,'point_speed%s%s.csv' %(a,b)))
    count = df['speed'].count()
    fraction = (count*100)/float(pnttotal)
    print 'fraction of points (speed %s %d) = %f percent (%d points of %d)' % (a,b,fraction,count,pnttotal)
    plt.figure()
    plt.ylim(0,0.2), plt.xlim(0,30)
    plotterhist(df['iri'].dropna(),'speed_%s_%s' % (a,b),d, pointsdf['iri'].dropna(), 'all points','%s%s' % (a,b))
    plt.clf()

def twoconditions(a,b,c,d): #, pointsdf['iri']<=b]
    if a == 'lessthan':
        e = (pointsdf['speed']<=b)
    elif a == 'morethan':
        e = (pointsdf['speed']>=b)
    else:
        print 'error: choose lessthan or morethan'
    if c == 'lessthan':
        f = (pointsdf['iri']<=d)
    elif c == 'morethan':
        f = (pointsdf['iri']>=d)
    else:
        print 'error: choose lessthan or morethan'
    df = pointsdf.loc[e & f]
    df.to_csv(os.path.join(OutPathSpeeds,'point_speed%s%s_iri%s%s.csv' %(a,b,c,d)))
    count = df['speed'].count()
    fraction = (count*100)/float(pnttotal)
    print 'fraction of points (speed %s %d, iri %s %d) = %f percent (%d points of %d)' % (a,b,c,d,fraction,count,pnttotal)

def speedoneconditionline(a,b):
    if a == 'lessthan':
        c = (linesdf['speed_med']<=b)
    elif a == 'morethan':
        c = (linesdf['speed_med']>=b)
    else:
        print 'error: choose lessthan or morethan'
    df = linesdf.loc[c]
    df.to_csv(os.path.join(OutPathSpeeds,'lines_speed%s%s.csv' %(a,b)))
    count = df['speed_med'].count()
    fraction = (count*100)/float(lntotal)
    print 'fraction of lines (speed %s %d) = %f percent (%d points of %d)' % (a,b,fraction,count,lntotal)

def longitudinal(a,IDLIST):
    n = 1
    for ID in IDLIST:
        MF = a.loc[a['VPROMMS_ID'] == ID]
        x = MF['Cm_Distance_metres']
        y = MF['%s' % z]
        if len(y) > 5:
            print 'creating long graph number %d...' % (n)
            plot = sns.regplot(x,y,fit_reg = False, label = r'%s, %s' % (ID, district)) #xlim = (0,), ylim = (0,)
            if answer == 1:
                plot.set(ylim=(0,10), xlim=(0,None))
            else:
                plot.set(ylim=(0,40), xlim=(0,None))
            plt.legend()
            plt.savefig("%s\\%s_%s_%s.png" % (OutPathLongs, district, ID,   z))
            plt.clf()
            print
            n = n + 1
        else:
            pass

def plotterkdejoint(Z,R):
    sns.jointplot(Z, R, kind="kde", size=7, space=0, dropna = True, xlim=(0,60), ylim=(0,20))
    plt.legend()
    plt.savefig("%s\\%s_%s_jointplot_kde.png" % (OutPathSpeeds, Z.name, R.name))
    plt.clf()


def plotterscatter(Z,R):
    sns.jointplot(Z, R, kind="scatter", size=7, space=0, dropna = True, xlim=(0,100), ylim=(0,60))
    plt.legend()
    plt.savefig("%s\\%s_%s_jointplot_scatter.png" % (OutPathSpeeds, Z.name, R.name))
    plt.clf()

def plotterhex(Z,R):
    sns.jointplot(Z, R, kind="hex", size=7, space=0, dropna = True, ylim=(0,40))
    plt.legend()
    plt.savefig("%s\\%s_%s_jointplot_hex.png" % (OutPathSpeeds, Z.name, R.name))
    plt.clf()

speed = pointsdf['speed']
iri = pointsdf['iri']
plotterkdejoint(speed,iri)
plotterscatter(speed,iri)
plotterhex(speed, iri)
twoconditions('lessthan',20,'lessthan',7)
twoconditions('lessthan',20,'lessthan',5)
twoconditions('morethan',40,'morethan',15)
twoconditions('morethan',40,'morethan',12)
speedonecondition('lessthan',20,'r')
speedonecondition('morethan',40,'y')
speedoneconditionline('lessthan', 20)
speedoneconditionline('morethan', 40)
statisticalreportpoints(pointsdf,linesdf)
longitudinal(pointsdf, IDLIST)
