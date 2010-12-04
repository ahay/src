from rsf.proj import *
import os,pbs


def processes(nodes=None,single=None):
    numThreads = int(os.environ.get('OMP_NUM_THREADS','8'))
    if nodes:
        return nodes*numThreads
    try:
        nodefile = os.environ['PBS_NODEFILE']
        ff = open(nodefile,'r')
        nodes = ff.readlines()
        ff.close()
        unique = []
        for node in nodes:
            if node not in unique:
                unique.append(node)
        return len(unique)*numThreads
    except:
        if single:
            return 1
        else:
            return numThreads

def _find(np,command,custom=''):
    ''' Find the mpiexec command, and the command to execute '''
    return '%s -np %d %s %s' % (WhereIs('mpiexec'),np,custom,WhereIs(command))

def encode(encodings, shotGathers, encoding, 
           np, 
           eprefix, dprefix,
           nx,ox,dx,ny,oy,dy,custom):
    ''' encode using sfbigmpiencode
    encodings - list of produced encoding files
    shotGathers - lsit of shotgathers to encode
    encoding - encoding produced by sfencodemaker
    np - number of processes to use
    eprefix - encoding prefix
    dprefix - data prefix
    nx,ox,dx,ny,oy,dy - output coordinates for encodings
    '''
    
    if not '.rsf' in eprefix:
        eprefix +='.rsf'
    if not '.rsf' in dprefix:
        dprefix +='.rsf'
    shotGathers.insert(0,encoding) 
    Flow(encodings, shotGathers,
        '''
        %s
        ''' % (_find(np,'sfbigmpiencode',custom)) + 
        ''' eprefix=''' + eprefix + 
        ''' dprefix=''' + dprefix + 
        '''
        encode=${SOURCES[0]}
        nx=%d ox=%f dx=%f
        ny=%d oy=%f dy=%f
        ''' % (nx,ox,dx,ny,oy,dy) ,stdin=0, stdout=-1)
        
def gridandstack(stack,files,np,
                fprefix,
                nx,ox,dx,
                ny,oy,dy,
                nz,oz,dz,
                nf=None,of=None,jf=None,
                shots=None):
    ''' stack files using sfbigencode, does not require files
    to be on the same cube, will relocate them in the cube
    
    stack - output file
    files - input files
    fprefix - input file prefix (dprefix)
    oprefix - output file prefix (eprefix)
    nx,ox,dx,ny,oy,dy - output stacked file dimensions
    either specify shots or nf,of,jf:
        shots - list of shot indices to be put into a file
        nf,of,jf - number of files in sequential order
    '''

    nfiles = len(files)

    oprefix=str(stack)

    if not '.rsf' in fprefix:
        fprefix += '.rsf'
    if not '.rsf' in oprefix:
        oprefix += '.rsf'

    if shots:

        def buildShotlist(target,source,env):
            shots = env['SHOTLIST']
            Fshotlist = open(str(target[0]),'w')
            for shot in shots:
                Fshotlist.write(str(shot)+'\n')
            Fshotlist.write('data_format=ascii_int\n')
            Fshotlist.write('esize=4\n')
            Fshotlist.write('in=./%s\n' % str(target[0]))
            Fshotlist.write('n1=%d\n' % len(shots))
            Fshotlist.write('o1=0\n')
            Fshotlist.write('d1=1\n')
            Fshotlist.close()

        bld = Builder(action= buildShotlist,
                      suffix='.asc')

        proj = Project()
        proj.Append(BUILDERS = {'ShotList' : bld})

        shotfile = stack+'-shots'
        proj.ShotList(shotfile,None,SHOTLIST=shots)

        if not '.asc' in shotfile:
            shotfile += '.asc'
        files.append(shotfile)
        Flow(stack,files,
            '''
            %s 
            ''' % (_find(np,'sfbigmpistack')) +
            '''
            nx=%d ny=%d nz=%d
            ox=%f oy=%f oz=%d
            dx=%f dy=%f dz=%f
            shots="%s"
            ''' % (nx,ny,nz,ox,oy,oz,dx,dy,dz,shotfile) + 
            '''
            prefix="'''+fprefix+'''" oname="'''+oprefix+'''"''',stdin=0,stdout=-1)

    else:
        if (not of) or (not jf) or (not nf):
            raise Exception('must specify either shots or nf,of,jf')
        Flow(stack,files,
            '''
            %s 
            ''' % (_find(np,'sfbigmpistack')) +
            '''
            nx=%d ny=%d nz=%d
            ox=%f oy=%f oz=%d
            dx=%f dy=%f dz=%f
            nf=%d of=%d jf=%d
            ''' % (nx,ny,nz,ox,oy,oz,dx,dy,dz,nf,of,jf) + 
            '''
            prefix="'''+fprefix+'''" oname="'''+oprefix+'''"''',stdin=0,stdout=-1)
     
        
def stack(stack,np,fprefix,nf,of,jf):
    ''' stack files using sfmpistack
    
    stack - output file
    np - number of processes
    fprefix - input file name prefix
    nf - number of files
    of - origin of files
    jf - delta of files
    see self-doc for more info
    '''
    

    filerange = range(of,of+nf*jf,jf)
    
    files = [ fprefix % f for f in filerange]
    
    
    if not '.rsf' in fprefix:
        fprefix +='.rsf'
    oname = stack
    if not '.rsf' in oname:
        oname += '.rsf'
        
    Flow(stack,files,
        '''
        %s
        nf=%d
        of=%d
        jf=%d
        ''' % (_find(np,'sfmpistack'),nf,of,jf) + 
        ''' fprefix="'''+fprefix + '''" oname="'''+oname+'''"''',stdin=0, stdout=-1)

