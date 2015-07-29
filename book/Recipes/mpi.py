cluster = False
try:
    from rsf.cluster import *
    if CSCONS:
        cluster = True
except:
    from rsf.proj import *
import subprocess


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

def WhereIs(command):
    p = subprocess.Popen('which %s' % command, shell=True, stdout=subprocess.PIPE)
    stdout,stderr = p.communicate()
    if len(stdout) == 0: raise Exception('could not find %s' % command)
    return stdout.strip('\n')

def _find(np,command,custom=''):
    ''' Find the mpiexec command, and the command to execute '''
    return '%s -np %d %s %s' % (WhereIs('mpiexec'),np,custom,WhereIs(command))


def encode(encodings,shotGathers,encoding,
           np, eprefix,dprefix,mpiopts='--bynode',**kw):
    ''' 
    encode using sfmpiencode.

    '''

    if not '.rsf' in eprefix: eprefix += '.rsf'

    if not '.rsf' in dprefix: dprefix += '.rsf'

    shotGathers.insert(0,encoding)

    command = _find(np,'sfmpiencode',mpiopts)

    global cluster
    if not cluster:
        kw = {}

    Flow(encodings,shotGathers,
        '''
        %s
        ''' % (command) + 
        '''eprefix='''+eprefix+''' dprefix=''' + dprefix + 
        '''
        encode=${SOURCES[0]}
        verb=y
        ''',stdin=0, stdout=-1,**kw)

def gridandencode(encodings, shotGathers, encoding, 
           np, 
           eprefix, dprefix,
           nx,ox,dx,ny,oy,dy,custom=None,
           mpiopts='--bynode',**kw):
    ''' 
    encode using sfbigmpiencode
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

    command = _find(np,'sfbigmpiencode',mpiopts)

    global cluster
    if not cluster:
        kw = {}

    Flow(encodings, shotGathers,
        '''
        %s
        ''' % (command) + 
        ''' eprefix=''' + eprefix + 
        ''' dprefix=''' + dprefix + 
        '''
        encode=${SOURCES[0]}
        nx=%d ox=%f dx=%f
        ny=%d oy=%f dy=%f
        ''' % (nx,ox,dx,ny,oy,dy), stdin=0,stdout=-1,**kw)
            
def gridandstack(stack,files,np,
                fprefix,
                nx,ox,dx,
                ny,oy,dy,
                nz,oz,dz,
                nf=None,of=None,jf=None,
                shots=None,
                mpiopts='--bynode', **kw):
    ''' stack files using sfbigencode, does not require files
    to be on the same cube, will relocate them in the cube
    
    stack - output file
    files - input files
    fprefix - input file prefix (dprefix)
    oprefix - output file prefix (eprefix)
    nx,ox,dx,ny,oy,dy - output stacked file dimensions
    either specify shots or nf,of,jf:
    -> shots - list of shot indices to be put into a file
    -> nf,of,jf - number of files in sequential order
    '''

    nfiles = len(files)

    oprefix=str(stack)

    if not '.rsf' in fprefix:
        fprefix += '.rsf'
    if not '.rsf' in oprefix:
        oprefix += '.rsf'

    global cluster
    if not cluster:
        kw = {}
    if shots:
        shotfile = stack+'-shots'
        Flow(shotfile,None,'points type=i out=${TARGETS[0]} x=%s' % reduce(lambda x,y: str(x)+','+str(y),shots))
       
        shotfile += '.rsf'
        files.append(shotfile)
        Flow(stack,files,
                '''
                %s 
                ''' % (_find(np,'sfbigmpistack',mpiopts)) + 
                '''
                nx=%d ny=%d nz=%d
                ox=%f oy=%f oz=%d
                dx=%f dy=%f dz=%f
                shots="%s"
                verb=y
                ''' % (nx,ny,nz,ox,oy,oz,dx,dy,dz,shotfile) + 
                '''
                prefix="'''+fprefix+'''" oname="'''+oprefix+'''"''',
                stdin=0, stdout=-1, **kw)
                    
    else:
        if of == None or jf == None or nf == None:
            raise Exception('must specify either shots or nf,of,jf')

        Flow(stack,files,
            '''
            %s 
            ''' % (_find(np,'sfbigmpistack',mpiopts)) +
            '''
            nx=%d ny=%d nz=%d
            ox=%f oy=%f oz=%d
            dx=%f dy=%f dz=%f
            nf=%d of=%d jf=%d
            ''' % (nx,ny,nz,ox,oy,oz,dx,dy,dz,nf,of,jf) + 
            '''
            prefix="'''+fprefix+'''" oname="'''+oprefix+'''"''',stdin=0,stdout=-1,**kw)
    
def stack(stack,np,fprefix=None,files=None,nf=None,of=None,jf=None,shots=None,mpiopts='--bynode',**kw):
    ''' stack files using sfmpistack
    
    stack - output file
    np - number of processes
    fprefix - input file name prefix
    nf - number of files
    of - origin of files
    jf - delta of files
    see self-doc for more info
    '''
   
    if not files:
        assert fprefix != None

        if not '.rsf' in fprefix:
            fprefix +='.rsf'

        if shots:
            files = [fprefix % x for x in shots]
        else:
            files = [fprefix % x for x in range(of,of+nf*jf,jf)]

    oname = stack
    if not '.rsf' in oname:
        oname += '.rsf'

    if shots:
        Flow(stack+'-shots',None,
            '''
            spike n1=%d nsp=%d mag=%s k1=%s | transp | sfdd type=int
            ''' % (len(shots),len(shots),','.join([str(shot) for shot in shots]),
                ','.join([str(i) for i in range(1,len(shots)+1)])
                    ))

    global cluster
    if not cluster:
        kw = {}

    if shots:
        files.insert(0,stack+'-shots')
        Flow(stack,files,
            '''
            %s shots=${SOURCES[0]}
            ''' % (_find(np,'sfmpistack',mpiopts)) + 
            ''' prefix="'''+fprefix + 
            '''" oname="'''+oname+'''"''',stdin=0,stdout=-1,**kw)

    elif nf != None and of != None and jf != None:
        Flow(stack,files,
            '''
            %s
            nf=%d
            of=%d
            jf=%d
            seq=y
            ''' % (_find(np,'sfmpistack',mpiopts),nf,of,jf) + 
            ''' prefix="'''+fprefix + 
            '''" oname="'''+oname+'''"''',stdin=0,stdout=-1,**kw)
    elif files:
        Flow(stack,files,
            '''
            %s ${SOURCES[:]}
            ''' % (_find(np,'sfmpistack',mpiopts)) + 
            ''' oname="'''+oname+'''"''',stdin=0,stdout=-1,**kw)
    else:
        raise Exception('must specify either files, nf of and jf or shots')
