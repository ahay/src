from rsfproj import *
import spmig, sgmig

# ------------------------------------------------------------
# SURVEY-SINKING migration
# ------------------------------------------------------------
def sinking(par):
    # surface data
    sgmig.wflds('d0','cmps',par)

    # datuned data
    sgmig.datum('d1','sd','d0',par)
    
    for k in ('0','1'):
        s = 's' + k # slowness
        d = 'd' + k # prestack data
        i = 'i' + k # prestack image
        e = 'e' + k # inner offset data
        z = 'z' + k # inner offset image

        # prestack migration
        sgmig.image(i,s,d,par)
        
        # near offset migration
        Flow(e,d,'window squeeze=n n3=8')
        sgmig.image(z,s,e,par)


# ------------------------------------------------------------
# SHOT-PROFILE MIGRATION
# ------------------------------------------------------------
def profile(par):

    # surface wavefields
    spmig.wflds('d0s','d0r','wave','shot',par)
    
    # datumed wavefields
    spmig.datum('d1s','d1r','sd','d0s','d0r',par)
    
    for k in ('0','1'):
        s = 's' + k       # slowness
        j = 'j' + k       # image
        ds= 'd' + k + 's' # source   wavefield
        dr= 'd' + k + 'r' # receiver wavefield
        
        spmig.image(j,s,ds,dr,par)
        
# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
def igrey(custom,par):
    return '''
    grey labelrot=n title="" %s min2=%g max2=%g min1=%g max1=%g
    ''' % (custom,par['xmin'],par['xmax'],par['zmin'],par['zmax'])

def result(par):

    for k in ('0','1','2','d','o'):
        s = 's' + k
        Result(s,s,'window      | transp |'
               +igrey('title=s pclip=100 color=j allpos=y',par))

        for k in ('0','1'):
            z = 'z' + k
            i = 'i' + k
            j = 'j' + k
            
            Result(z,z,'window n3=1 | transp |'
                   +igrey('title=z pclip=99',par))
            Result(i,i,'window n3=1 | transp |'
                   +igrey('title=i pclip=99',par))
            Result(j,j,'window      | transp |'
                   +igrey('title=j pclip=99',par))

# ------------------------------------------------------------
