from rsf.proj import *

def _find(np,command):
    ''' Find the mpiexec command, and the command to execute '''
    return '%s -np %d %s' % (WhereIs('mpiexec'),np,WhereIs(command))

def encode(encodings, shotGathers, encoding, 
           np, 
           eprefix, dprefix,
           nx,ox,dx,ny,oy,dy):
    ''' encode using sfbigmpiencode
    encodings - list of produced encoding files
    shotGathers - lsit of shotgathers to encode
    encoding - encoding produced by sfencodemaker
    np - number of processes to use
    eprefix - encoding prefix
    dprefix - data prefix
    nx,ox,dx,ny,oy,dy - output coordinates for encodings
    '''
    shotGathers.insert(0,encoding) 
    Flow(encodings, shotGathers,
        '''
        %s
        eprefix=%s
        dprefix=%s
        ''' % (_find(np,'sfbigmpiencode'),eprefix,dprefix) + 
        '''
        encode=${SOURCES[0]}
        nx=%d ox=%f dx=%f
        ny=%d oy=%f dy=%f
        ''' % (nx,ox,dx,ny,oy,dy) ,stdin=0, stdout=-1)
        
def stack(stack,files,np,fprefix,oprefix,nx,ox,dx,ny,oy,dy):
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
        
