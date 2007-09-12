from rsfproj import *

# ------------------------------------------------------------
# lags to slant-stacks

# input:  z-h-x
# output: z-tan(a)-x
def xcig2ssk(np,op,dp):
    return '''
    slant adj=y np=%d p0=%g dp=%g verb=y
    ''' % (np,op,dp) 

# input: z-tan(a)-x
# output z-a-x
def tcig2ssk(np,op,dp):
    return '''
    slant adj=y np=%d p0=%g dp=%g verb=y
    ''' % (np,op,dp) 

# ------------------------------------------------------------
# slant-stacks to angle
def xssk2ang(na,oa,da):
    return '''
    tan2ang na=%d a0=%g da=%g
    ''' % (na,oa,da)

def tssk2ang(na,oa,da):
    return '''
    tshift cos=y na=%d a0=%g da=%g
    velocity=${SOURCES[1]} dip=${SOURCES[2]}
    ''' % (na,oa,da)

# ------------------------------------------------------------
def agrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100 gainpanel=a
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=2
    label2="\F10 q\F3 " unit2="\^o\_"
    %s
    ''' % (
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['labelattr']+' '+custom)

def xgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100 gainpanel=a
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=2
    label2="\F10 l\F3 " unit2=%s
    %s
    ''' % (
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['ux'],
           par['labelattr']+' '+custom)

def tgrey(custom,par):
    return '''
    grey labelrot=n wantaxis=y title=""
    pclip=100 gainpanel=a
    min1=%g max1=%g label1=%s unit1=%s
    screenratio=2
    label2="\F10 t\F3 " unit2=%s
    %s
    ''' % (
           par['zmin'],par['zmax'],par['lz'],par['uz'],
           par['ut'],
           par['labelattr']+' '+custom)
