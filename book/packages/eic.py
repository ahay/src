from rsfproj import *


# ------------------------------------------------------------
def cic(imag,swfl,rwfl,custom,par):
    par['ciccustom'] = custom

    Flow(imag,[swfl,rwfl],
         'xcor2d uu=${SOURCES[1]} axis=3 verb=y nbuf=10 ompnth=%(ompnth)d' % par)


# ------------------------------------------------------------
def eic(cip,swfl,rwfl,cc,custom,par):
    par['eiccustom'] = custom
    
    Flow(cip,[swfl,rwfl,cc],
         '''
         laps2d verb=y
         nhx=%(nhx)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
         ur=${SOURCES[1]}
         cc=${SOURCES[2]}
         %(eiccustom)s
         ''' %par)

