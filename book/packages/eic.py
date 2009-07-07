from rsfproj import *

# ------------------------------------------------------------
def eic(cip,swfl,rwfl,cc,custom,par):
    par['iccustom'] = custom
    
    Flow(cip,[swfl,rwfl,cc],
         '''
         laps2d verb=y
         nhx=%(nhx)d nhz=%(nhz)d nht=%(nht)d dht=%(dht)g
         ur=${SOURCES[1]}
         cc=${SOURCES[2]}
         %(iccustom)s
         ''' %par)

