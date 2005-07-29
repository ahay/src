from rsfproj import *

# generate launch file
def launch(edir,job,nodes,time,queue,project):
    launcher = job + '.launch'
    launcher = edir + '/' + launcher
    mycom    = 'rm '      + launcher
    os.system(mycom)
    
    echoline(launcher,'#!/bin/csh')
    echoline(launcher,'#BSUB -J ' + job)
    echoline(launcher,'#BSUB -o ' + job + '.%J.out')
    echoline(launcher,'#BSUB -e ' + job + '.%J.err')
    echoline(launcher,'#BSUB -n ' + str(nodes))
    echoline(launcher,'#BSUB -W ' + str(time))
    echoline(launcher,'#BSUB -q ' + queue)
    echoline(launcher,'#BSUB -P ' + project)
    echoline(launcher,'#BSUB -u paul.sava@gmail.com')
    echoline(launcher,'#BSUB -B')
    echoline(launcher,'#BSUB -N')
#    echoline(launcher,'#BSUB -R "span[ptile=1]"')
#    echoline(launcher,'#BSUB -x')	
    echoline(launcher,'setenv EXECUTABLE     launcher')
    echoline(launcher,'setenv WORKDIR        .')
    echoline(launcher,'setenv CONTROL_FILE   ' + job)

    mycom = 'cat ../../Scons/LAUNCH >> ' + launcher
    os.system(mycom)

def echoline(launcher,line):
    mycom ='echo \'' + line + '\' >>' + launcher
    os.system(mycom)

# ------------------------------------------------------------

def submit(EDIR,JOB,nodes,time,queue,project):
    launcher = EDIR + '/' + JOB + '.launch'

    f = open(launcher,'w')

    f.write('#!/bin/csh' + '\n')
    f.write('#BSUB -J '  + JOB + '\n')
    f.write('#BSUB -o '  + JOB + '.%J.out' + '\n')
    f.write('#BSUB -e '  + JOB + '.%J.err' + '\n')
    f.write('#BSUB -n '  + str(nodes) + '\n')
    f.write('#BSUB -W '  + str(time)  + '\n')
    f.write('#BSUB -q '  + queue      + '\n')
    f.write('#BSUB -P '  + project    + '\n')
    f.write('setenv  EXECUTABLE     launcher' + '\n')
    f.write('setenv  WORKDIR               .' + '\n')
    f.write('setenv  CONTROL_FILE '     + JOB + '\n')

    f.close()

    mycom = 'cat ../../Scons/LAUNCH >> ' + launcher
    os.system(mycom)

# ------------------------------------------------------------

def sprd(file,EDIR,axis,ngroup,nelems,jelems):
    allk = map(lambda x: '%03d' % x,range(ngroup))
    for k in allk:
        _f =  '_' + file + '.' + k + '.rsf'

        oelems = nelems * int(k)*jelems

        axiswin = ' n'+str(axis)+'='+str(nelems)
        axiswin+= ' f'+str(axis)+'='+str(oelems)
        axiswin+= ' j'+str(axis)+'='+str(jelems)

        Flow(EDIR+_f,file,'window squeeze=n %s out=stdout' % axiswin )

def copy(file,EDIR):
    _f = '_' + file + '.rsf'
    Flow(EDIR + _f,file,'window squeeze=n out=stdout')

# add image files
def summ(file,EDIR,ngroup):
    one = EDIR + '_' + file + '.%03d.rsf' 
    all = map(lambda x: one % x,range(ngroup))
    Flow(file,all,'add ${SOURCES[:%d]}' % len(all), stdin=0)

# ------------------------------------------------------------

# split input for S-G migration
def SGsplit(edir,wfld,ngroup,nfreqs,jfreqs):

    allk = map(lambda x: '%03d' % x,range(ngroup))
    for k in allk:
        _w = '_' + wfld + '.' + k + '.rsf'

        ofreq = nfreqs * int(k)*jfreqs

        Flow(edir + _w,
             wfld+'.rsf',
             'window squeeze=n n4=%d f4=%d j4=%d out=stdout' % (nfreqs,ofreq,jfreqs) )

# split input for S-R migration
def SRsplit(edir,swfl,rwfl,ngroup,nshots,jshots):

    allk = map(lambda x: '%03d' % x,range(ngroup))
    for k in allk:
        _s = '_' + swfl + '.' + k + '.rsf'
        _r = '_' + rwfl + '.' + k + '.rsf'

        oshot = nshots * int(k)*jshots

        Flow(edir + _s,swfl+'.rsf',
             'window squeeze=n n4=%d f4=%d j4=%d out=stdout' % (nshots,oshot,jshots) )
        Flow(edir + _r,rwfl+'.rsf',
             'window squeeze=n n4=%d f4=%d j4=%d out=stdout' % (nshots,oshot,jshots) )

# add image files
def outadd(edir,out,ngroup):
    oneout = edir + '_' + out + '.%03d.rsf' 
    allout = map(lambda x: oneout % x,range(ngroup))
    Flow(out,allout,'add ${SOURCES[:%d]}' % len(allout), stdin=0)

