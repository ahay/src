from rsf.proj import *
import rsf.gallery 

method = rsf.gallery.method()

Fetch('french.asc','french')
Flow('french','french.asc','dd form=native | transp | scale dscale=2')

def get_refl(refl):
    Flow(refl,'french',
         '''
         remap1 n1=161 o1=0 d1=51.325 | transp |
         remap1 n1=161 o1=0 d1=51.325 | transp
         ''')
    
def get_zodata(data):
    refl = data+'-refl'
    dipx = data+'-dipx'
    dipy = data+'-dipy'

    get_refl(refl)
    Flow(dipx,refl,'smooth rect1=5 | deriv scale=y')
    Flow(dipy,refl,'transp | smooth rect1=5 | deriv scale=y | transp')


    Flow(data,[refl,dipx,dipy],
         '''
         kirmod3 nt=601 dt=0.010265 vel=2000
         dipx=${SOURCES[1]} dipy=${SOURCES[2]}
         s0x=0 dsx=51.325 nsx=161 verb=n
         s0y=0 dsy=51.325 nsy=161
         nhx=1 h0x=0 dhx=51.325
         nhy=1 h0y=0 dhy=51.325 freq=5 |
         window |
         put label1=Time unit1=s
         label2=North-South unit2=m
         label3=West-East   unit3=m
         ''',split=[2,161],reduce='add')

