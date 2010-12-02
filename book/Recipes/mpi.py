from rsf.proj import *
import os

def processes(nodes):
    return nodes*int(os.environ.get('OMP_NUM_THREADS','8'))

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
        
def gridandstack(stack,files,np,fprefix,oprefix,nx,ox,dx,ny,oy,dy):
    ''' stack files using sfbigencode, does not require files
    to be on the same cube, will relocate them in the cube
    
    stack - output file
    files - input files
    fprefix - input file prefix (dprefix)
    oprefix - output file prefix (eprefix)
    nx,ox,dx,ny,oy,dy - output stacked file dimensions
    '''
    
    nfiles = len(files)
    Flow(stack+'-encode-amp',None,
        '''
        spike n1=%d o1=0 d1=1 n2=1 o2=0 d2=1 |
        math output="1"
        ''' % nfiles)
    Flow(stack+'-encode-pha',None,
        '''
        spike n1=%d o1=0 d1=1 n2=1 o2=0 d2=1 |
        math output="0"
        ''' % nfiles)
    Flow(stack+'-encode',[stack+'-encode-amp',stack+'-encode-pha'],
        '''
        cat axis=3 ${SOURCES[1]}
        ''')
    encode(stack,files,stack+'-encode',
        np,oprefix,fprefix,nx,ox,dx,ny,oy,dy)
        
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
        
    Flow(stack,files,
        '''
        %s
        nf=%d
        of=%d
        jf=%d
        ''' % (_find(np,'sfmpistack'),nf,of,jf) + 
        ''' fprefix='''+fprefix + ''' oname='''+stack,stdin=0, stdout=-1)
        
