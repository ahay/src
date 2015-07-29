from rsf.proj import *
import spmig, sgmig, zomig,fdmod


# ------------------------------------------------------------
# RTM
# ------------------------------------------------------------

# WR: forward in time
def fWRrtm(data,wfld,cccc,dens,coor,custom,par):
    iwindow = ' ' + \
        '''
        nqz=%(nqz)d oqz=%(oqz)g
        nqx=%(nqx)d oqx=%(oqx)g
        jsnap=%(jdata)d jdata=%(jdata)d
        ''' % par + ' '

    fdmod.ewefd(wfld+'_out',wfld,
                 data,cccc,dens,coor,coor,iwindow+custom,par)
   



    # WR: backward in time
def bWRrtm(data,wfld,velo,dens,coor,custom,par):
    iwindow = ' ' + \
        '''
        nqz=%(nqz)d oqz=%(oqz)g
        nqx=%(nqx)d oqx=%(oqx)g
        jsnap=%(jdata)d jdata=%(jdata)d
        ''' % par + ' '

    Flow(data+'_rev',data,'reverse which=2 opt=i verb=y')
    fdmod.anifd2d(wfld+'_out',wfld+'_rev',
                  data+'_rev',velo,dens,coor,coor,iwindow+custom,par)
    Flow(wfld,wfld+'_rev','reverse which=4 opt=i verb=y')
