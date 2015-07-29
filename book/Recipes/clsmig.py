import os

# prepare for shot-profile migration
def sppre(ngroup):

    allsou = map(lambda x: '_s%03d.rsf' % x,range(ngroup))
    allrec = map(lambda x: '_r%03d.rsf' % x,range(ngroup))
    alljmg = map(lambda x: '_j%03d.rsf' % x,range(ngroup))
    
    silent = ' -s '
    #silent = ' -Q '
    
    for sou in allsou:
        mycom = 'scons' + silent + sou
        os.system(mycom)
        
    for rec in allrec:
        mycom = 'scons' + silent + rec
        os.system(mycom)
            
    for jmg in alljmg:
        mycom = 'scons -n -Q ' + jmg
        os.system(mycom)

# prepare for survey sinking migration
def sgpre(ngroup):

    alldat = map(lambda x: '_e%03d.rsf' % x,range(ngroup))
    allimg = map(lambda x: '_i%03d.rsf' % x,range(ngroup))
    
    silent = ' -s '
    #silent = ' -Q '

    for dat in alldat:
        mycom = 'scons' + silent + dat
        os.system(mycom)

    for img in allimg:
        mycom = 'scons -n -Q ' + img
        os.system(mycom)
