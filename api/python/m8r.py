##   Copyright (C) 2010 University of Texas at Austin
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals
'''
Python interface to Madagascar. Serves dual purpose:
- an ability to manipulate RSF file objects (File class), figures (Vplot class) and Madagascar programs (Filter class) in a script environment.
- an ability to write Madagascar programs in Python (using Par, Input, and Output classes).
'''

import os, sys, tempfile, re, subprocess, urllib, datetime
import numpy as np

import rsf.doc
import rsf.prog
import rsf.path

# There are two parallel implementations:
# - one is using calls to C functions from api/c with SWIG interface
# - the other is pure Python
# Find out which one we have
try:
    import c_m8r as c_rsf
    _swig_ = True
    # to switch off SWIG, use no_swig() function defined below
except:
    _swig_ = False

python2 = sys.version_info[0] < 3
    
def view(name):
    'for use in Jupyter notebooks'
    try:
        from IPython.display import Image
    except:
        print ('No IPython Image support.')
        return None
    error_message = 'Failed to generate image.'
    try:
        png = name+'.png'
        makefile = os.path.join(rsf.prog.RSFROOT,'include','Makefile')
        proc = subprocess.Popen('make -f %s %s' % (makefile,png),
                                shell=True,
                                stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE)
        error_message = proc.stderr.read().decode()
        return Image(filename=png)
    except:
        print(error_message,file=sys.stderr)
        return None

class _Simtab(dict):
    'symbol table (emulates api/c/simtab.c)'
    # inherited from dict, because it is basically a dictionary
    def enter(self,key,val):
        'add key=val to the table'
        self[key] = val
    def getint(self,key):
        val = self.get(key)
        if val:
            return True,int(val)
        return False,None
    def getints(self,key,n):
        val = self.get(key)
        if val:
            vals = val.split(',')
            nval = len(vals)
            # set array to length n
            if n < nval:
                vals = vals[:n]
            elif n > nval:
                vals.extend((n-nval)*[vals[nval-1]])
            return True,[int(v) for v in vals]
        return False,None
    def getfloat(self,key):
        val = self.get(key)
        if val:
            return True,float(val)
        return False,None
    def getfloats(self,key,n):
        val = self.get(key)
        if val:
            vals = val.split(',')
            nval = len(vals)
            # set array to length n
            if n < nval:
                vals = vals[:n]
            elif n > nval:
                vals.extend((n-nval)*[vals[nval-1]])
            return True,[float(v) for v in vals]
        return False,None
    def getstring(self,key):
        val = self.get(key)
        if None == val:
            return None
        # strip quotes
        if val[0] == '"' and val[-1] == '"':
            val = val[1:-1]
        return val
    def getbool(self,key):
        val = self.get(key)
        if val:
            if val[0] == 'y' or val[0] == 'Y' or val[0] == '1':
                return True,True
            else:
                return True,False
        return False,None
    def getbools(self,key,n):
        val = self.get(key)
        if val:
            vals = val.split(',')
            nval = len(vals)
            # set array to length n
            if n < nval:
                vals = vals[:n]
            elif n > nval:
                vals.extend((n-nval)*[vals[nval-1]])
            bools = []
            for val in vals:
                bools.append(val[0] == 'y' or val[0] == 'Y' or val[0] == '1')
            return True,bools
        return False,None
    def put(self,keyval):
        if '=' in keyval:
            key,val = keyval.split('=')
            self.enter(key,val)
    def string(self,string):
        'extract parameters from a string'
        for word in string.split():
            self.put(word)
    def input(self,filep,out=None):
        'extract parameters from header file'
        # Special code b'\x0c\x0c\x04', if encountered, signifies
        # the end of the header and the start of the data.
        # With each new line, we will try to read the first three
        # bytes and compare them to the code before reading the
        # rest of the line.
        if python2:
            stream = filep
        else:
            stream = filep.buffer
        while True:
            try:
                line3 = stream.read(3)
                # skip new lines
                while line3[:1] == b'\n':
                    line3 = line3[1:] + stream.read(1)
                # check code for the header end
                if line3 == b'\x0c\x0c\x04':
                    break
                line = line3+stream.readline()
                if not python2:
                    # bytes to string
                    line = str(line,'utf-8')
                if len(line) < 1:
                    break
                if out:
                    out.write(line)
                # extract parameters
                self.string(line)
            except:
                break
        if out:
            out.flush()
    def output(self,filep):
        'output parameters to a file'
        for key in self.keys():
            filep.write('\t%s=%s\n' % (key,self[key]))
 
class Par(object):
    '''command-line parameter table'''
    def __init__(self,argv=sys.argv):
        global par
        if _swig_:
            # c_rsf.sf_init needs 'utf-8' (bytes) not unicode
            c_rsf.sf_init(len(argv),
                            [arg.encode('utf-8') for arg in argv])
            self.prog = c_rsf.sf_getprog()
        else:
            self.pars = _Simtab()
            self.prog = argv[0]
            
            for arg in argv[1:]:
                if arg[:4] == 'par=':
                    # extract parameters from parameter file
                    parfile = open(arg[4:],'r')
                    self.pars.input(parfile)
                    parfile.close()
                else:
                    self.pars.put(arg)
            
        for type in ('int','float','bool'):
            # create methods int, float, bool
            setattr(self,type,self.__get(type))
            # create methods ints, floats, bools
            setattr(self,type+'s',self.__gets(type))
        par = self
    def getprog(self):
        return self.prog
    def __get(self,type):
        # factory for parameter extraction
        if _swig_:
            func = getattr(c_rsf,'sf_get'+type)
        else:
            func = getattr(self.pars,'get'+type)
        def _get(key,default=None):
            if python2 and _swig_:
                # c function only knows utf-8 (ascii).  translate the unicode
                key = key.encode('utf-8')
            get,par = func(key)
            if get:
                return par
            elif default != None:
                return default
            else:
                return None
        return _get
    def __gets(self,type):
        # factory for parameter extraction for lists
        if _swig_:
            func = getattr(c_rsf,'sf_get'+type+'s')
        else:
            func = getattr(self.pars,'get'+type+'s')
        def _gets(key,num,default=None):
            if python2 and _swig_:
                # c function only knows utf-8 (ascii).  translate the unicode
                key = key.encode('utf-8')
            pars = func(key,num)
            if pars:
                return pars
            elif default:
                return default
            else:
                return None
        return _gets
    def string(self,key,default=None):        
        if _swig_:
            if python2:
                # c function only knows utf-8 (ascii).  translate the unicode
                key = key.encode('utf-8')
            val = c_rsf.sf_getstring(key)
        else:
            val = self.pars.getstring(key)
    
        if val:
            return val
        elif default:
            return default
        else:
            return None

# set default parameters for interactive scripts
par = Par(['python','-'])

def no_swig():
    'disable swig for testing'
    global _swig_, par
    _swig_ = False
    par = Par(['python','-'])
        
class Temp(str):
    'Temporaty file name'
    # inherits from string
    datapath = rsf.path.datapath()
    tmpdatapath = os.environ.get('TMPDATAPATH',datapath)
    def __new__(cls):
        return str.__new__(cls,tempfile.mktemp(dir=Temp.tmpdatapath))

class File(object):
    'generic RSF file object'
    attrs = ['rms','mean','norm','var','std','max','min','nonzero','samples']
    def __init__(self,tag,temp=False,name=''):
        'Constructor'
        self.temp = temp
        if isinstance(tag,File):
            # copy file (name is ignored)
            self.__init__(tag.tag)
        elif isinstance(tag,np.ndarray):
            # numpy array
            dtype = tag.dtype.kind
            if dtype=='f':
                dformat='native_float'
            elif dtype=='i':
                dformat='native_int'
            elif dtype=='u':
                dformat='native_uchar'
            else:
                raise TypeError('Unsupported type: %s' % dtype)
            if not name:
                name = Temp()
            out = Output(name,data_format=dformat)
            shape = tag.shape
            dims = len(shape)
            for axis in range(1,dims+1):
                # in numpy, axis order is reversed
                out.put('n%d' % axis,shape[dims-axis])
            out.write(tag)
            out.close()
            self.__init__(out,temp=True)
        elif isinstance(tag,list):
            # convert list to numpy array before importing
            if type(tag[0]) == int:
                dtype='int32'
            elif type(tag[0]) == float:
                dtype='float32'
            else:
                raise TypeError('Unsupported format %s' % type(tag[0]))
            array = np.array(tag,dtype)
            self.__init__(array,True,name)
        else:
            self.tag = tag
        self.narray = []
        for filt in Filter.plots + Filter.diagnostic:
            # run things like file.grey() or file.attr()
            setattr(self,filt,self.__filter(filt))
        for attr in File.attrs:
            # extract atrributes 
            setattr(self,attr,self.want(attr))
    def __filter(self,filt):
        'apply self-filter'
        def _filter(**kw):
            command = Filter(filt)(**kw)
            return command.apply(self)
        return _filter
    def __str__(self):
        'String representation'
        if self.tag:
            tag = str(self.tag)
            if os.path.isfile(tag):
                return tag
            else:
                raise TypeError('Cannot find "%s" ' % tag)
        else:
            raise TypeError('Cannot find tag')
    def sfin(self):
        'Output of sfin'
        return Filter('in',stdout=False).apply(self)
    def want(self,attr):
        'Attributes from sfattr'
        def wantattr():
            try:
                val = os.popen('%s want=%s < %s' %
                               (Filter('attr'),attr,self)).read()
            except:
                raise RuntimeError('trouble running sfattr')
            m = re.search('=\s*(\S+)',val)
            if m:
                val = float(m.group(1))
            else:
                raise RuntimeError('no match')
            return val
        return wantattr
    def real(self):
        'Take real part'
        return Filter('real').apply(self)
    def cmplx(self,im):
        return Filter('cmplx').apply(self,im)
    def imag(self):
        'Take imaginary part'
        return Filter('imag').apply(self)
    def __add__(self,other):
        'Overload addition'
        return Filter('add').apply(self,other)
    def __sub__(self,other):
        'Overload subtraction'
        sub = Filter('add')(scale=[1,-1])
        return sub.apply(self,other)
    def __mul__(self,other):
        'Overload multiplication'
        try:
            mul = Filter('scale')(dscale=float(other))
            return mul.apply(self)
        except:
            mul = Filter('mul')
            return mul.apply(self,other)
    def __div__(self,other):
        'Overload division'
        try:
            div = Filter('scale')(dscale=1.0/float(other))
            return div.apply(self)
        except:
            div = Filter('div')
            return div.apply(self,other)
    def __neg__(self):
        neg = Filter('scale')(dscale=-1.0)
        return neg.apply(self)
    def dot(self,other):
        'Dot product'
        # incorrect for complex numbers
        prod = self * other
        stack = Filter('stack')(norm=False,axis=0).apply(prod)
        return stack[0]
    def cdot2(self):
        'Dot product with itself'
        filt = Filter('math')(output="\"input*conj(input)\"").real.stack(norm=False,axis=0)
        stack = filt.apply(self)
        return stack[0]
    def dot2(self):
        'Dot product with itself'
        return self.dot(self)
    def __array__(self,context=None):
        'create narray'
        if [] == self.narray:
            if _swig_:
                if not hasattr(self,'file'):
                    f = c_rsf.sf_input(self.tag)
                else:
                    f = self.file
                self.narray = c_rsf.rsf_array(f)
                if not hasattr(self,'file'):
                    c_rsf.sf_fileclose(f)
            else:
                # gets only the real part of complex arrays ##kls
                if not hasattr(self,'file'):
                    f=Input(self.tag)
                else:
                    f=self.file
                self.narray=np.memmap(f.string('in'),dtype=f.datatype,
                                      mode='r+',
                                      shape=f.shape())
        return self.narray
    def __array_wrap__(self,array,context=None):
        inp = Input(self)
        inp.read(array)
        return inp
    def __getitem__(self,i):
        array = self.__array__()
        return array[i]
    def __setitem__(self,index,value):
        array = self.__array__()
        array.__setitem__(index,value)
    def size(self,dim=0):
        return File.leftsize(self,dim)

    def leftsize(self,dim=0):
        if _swig_:
            if hasattr(self,'file'):
                f = self.file
            else:
                f = c_rsf.sf_input(self.tag)
            s = c_rsf.sf_leftsize(f,dim)
            if not hasattr(self,'file'):
                c_rsf.sf_fileclose(f)
            return s
        else:
            s = 1
            for axis in range(dim+1,10):
                n = self.int("n%d" % axis)
                if n:
                    s *= n
                else:
                    break
            return s
    def int(self,key,default=None):
        get = self.get(key)
        if get:
            val = int(get)
        elif default:
            val = default
        else:
            val = None
        return val
    def float(self,key,default=None):
        get = self.get(key)
        if get:
            val = float(get)
        elif default:
            val = default
        else:
            val = None
        return val
    def get(self,key):
        'returns a string'
        try:
            p = subprocess.Popen('%s %s parform=n < %s' %
                                 (Filter('get'),key,self),
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 close_fds=True)
            result = p.stdout.read()
        except:
            raise RuntimeError('trouble running sfget')
        return result
    def shape(self):
        s = []
        dim = 1
        for i in range(1,10):
            ni = self.int('n%d' % i)
            if ni:
                dim = i
            s.append(ni)
        s = s[:dim]
        # the trailing members of s that are 1 ie fix situations like
        # s=(1500,240,1,1)
        while s[-1]==1 and len(s)>1:
            s=s[:-1]
        # axes are reversed for consistency with numpy
        s.reverse()
        return tuple(s)
    def reshape(self,shape=None):
        if not shape:
            shape = self.size()
        try:
            shape = list(shape)
        except:
            shape = [shape]
        old = list(self.shape())
        old.reverse()
        shape.reverse()
        lold = len(old)
        lshape = len(shape)
        puts = {}
        for i in range(max(lold,lshape)):
            ni = 'n%d' % (i+1)
            if i < lold:
                if i < lshape:
                    if old[i] != shape[i]:
                        puts[ni] = shape[i]
                else:
                    puts[ni] = 1
            else:
                puts[ni] = shape[i]
        put = Filter('put')
        put.setcommand(puts)
        return put.apply(self)
    def __len__(self):
        return self.size()
    def __del__(self):
        # remove temporary files
        if hasattr(self,'temp') and self.temp:
            Filter('rm',stdin=0,stdout=0).apply(self)

class _File(File):
    'umbrella class for RSF files that can read and write' 
    types = ['uchar','char','int','float','complex']
    forms = ['ascii','xdr','native']
    def __init__(self,tag):
        if not self.file:
            # can only be used by Input and Output
            # that define self.file
            raise TypeError('Use Input or Output instead of File')
        self.type_form()
        if self.type=='float':
            self.datatype=np.float32
        elif self.type=='complex':
            self.datatype=np.complex64
        elif self.type=='int':
            self.datatype=np.int32
        elif self.type=='uchar':
            self.datatype=np.uint8
        else:
            raise TypeError('Unsupported type %s' % self.type)
        File.__init__(self,tag)

        # create methods line int, float, and ints
        # for extracting parameters
        for type in ('int','float'):
            setattr(self,type,self.__get(type))
        for type in ('int',):
            setattr(self,type+'s',self.__gets(type))
    def type_form(self):
        if _swig_:
            self.type = _File.types[c_rsf.sf_gettype(self.file)]
            self.form = _File.forms[c_rsf.sf_getform(self.file)]
        else:
            self.type = self.file.gettype()
            self.form = self.file.getform()
    def tell(self):
        if _swig_:
            return c_rsf.sf_tell(self.file)
        else:
            return self.file.tell()
    def close(self):
        if _swig_:
            c_rsf.sf_fileclose(self.file)
        else:
            self.file.fileclose()
    def __del__(self):
        if self.file:
            if _swig_:
                c_rsf.sf_fileclose(self.file)
            else:
                self.file.fileclose()
        File.__del__(self) # this removes file if it is temporary
    def settype(self,type):
        if _swig_:
            for i,filetype in enumerate(_File.types):
                if type == filetype:
                    self.type = type
                    c_rsf.sf_settype(self.file,i)
                    break
        else:
            self.file.settype(type)
            self.type = type
    def setformat(self,format):
        if _swig_:
            c_rsf.sf_setformat(self.file,format)
        else:
            self.file.setformat(format)
        self.type_form()
    def __get(self,type):
        if _swig_:
            func = getattr(c_rsf,'sf_hist'+type)
        else:
            func = getattr(self.file,type)
        def _get(key,default=None):
            if _swig_:
                if python2:
                    # c function only knows utf-8 (ascii).  translate the unicode
                    key = key.encode('utf-8')
                get,par = func(self.file,key)
            else:
                get,par = func(key)   
            if get:
                return par
            elif default:
                return default
            else:
                return None
        return _get
    def __gets(self,type):
        if _swig_:
            func = getattr(c_rsf,'sf_hist'+type+'s')
        else:
            func = getattr(self.file,type)
        def _gets(key,num,default=None):
            if _swig_:
                if python2:
                    # c function only knows utf-8 (ascii).  translate the unicode
                    key = key.encode('utf-8')
                pars = func(self.file,key,num)
            else:
                pars = func(key,num)
 
            if pars:
                return pars
            elif default:
                return default
            else:
                return None
        return _gets
    def string(self,key):
        if _swig_:
            if python2:
                # c function only knows utf-8 (ascii).  translate the unicode
                key = key.encode('utf-8')
            return c_rsf.sf_histstring(self.file,key)
        else:
            return self.file.string(key)
    def bytes(self):
        if _swig_:
            return c_rsf.sf_bytes(self.file)
        else:
            return self.file.bytes()
    def put(self,key,val):
        if _swig_:
            if python2:
                # c function only knows utf-8 (ascii).  translate the unicode
                key = key.encode('utf-8')
            if isinstance(val,int):
                c_rsf.sf_putint(self.file,key,val)
            elif isinstance(val,np.int32) or isinstance(val,np.int64):
                c_rsf.sf_putint(self.file,key,int(val))
            elif isinstance(val,float):
                c_rsf.sf_putfloat(self.file,key,val)
            elif isinstance(val,str):
                c_rsf.sf_putstring(self.file,key,val)
            elif isinstance(val,unicode):
                c_rsf.sf_putstring(self.file,key,str(val))
            elif isinstance(val,list):
                if isinstance(val[0],int):
                    c_rsf.sf_putints(self.file,key,val)
            else:
                raise TypeError('Unsupported type in put: %s' % type(val))
        else:
            if isinstance(val,int):
                self.file.putint(key,val)
            elif isinstance(val,np.int32) or isinstance(val,np.int64):
                self.file.putint(key,int(val))
            elif isinstance(val,float):
                self.file.putfloat(key,val)
            elif isinstance(val,str):
                self.file.putstring(key,val)
            elif isinstance(val,unicode):
                self.file.putstring(key,str(val))
            elif isinstance(val,list):
                if isinstance(val[0],int):
                    self.file.putints(key,val)
            else:
                raise TypeError('Unsupported type in put: %s' % type(val))
    def axis(self,i):
        ax = {}
        ax['n'] = self.int("n%d"   % i)
        ax['d'] = self.float("d%d" % i)
        ax['o'] = self.float("o%d" % i)
        ax['l'] = self.string("label%d" % i)
        ax['u'] = self.string("unit%d" % i)
        return ax
    def putaxis(self,ax,i):
        self.put("n%d" % i,ax['n'])
        self.put("d%d" % i,ax['d'])
        self.put("o%d" % i,ax['o'])
        if ax['l']:
            self.put("label%d" % i,ax['l'])
        if ax['u']:
            self.put("unit%d" % i,ax['u'])

class Input(_File):
    '''
    RSF file for reading.  
    '''
    def __init__(self,tag='in'):
        self.file = None
        if isinstance(tag,File):
            # copy file
            self.__init__(tag.tag)
        else:
            if _swig_:
                if python2:
                    # c function only knows utf-8 (ascii).  translate the unicode
                    tag = tag.encode('utf-8')
                self.file = c_rsf.sf_input(tag)
            else:
                self.file = _RSF(True,tag)
        _File.__init__(self,tag)
    def read(self,data=[],shape=None,datatype=None):
        if len(data) == 0:
            allocate = True
            if shape==None:
                shape=self.shape()
            if datatype==None:
                datatype=self.datatype
            data=np.zeros(shape,dtype=datatype)
        else:
            allocate = False
        shape=data.shape
        datacount=data.size
        if _swig_:
            if self.type == 'float':
                c_rsf.sf_floatread(np.reshape(data,(data.size,)),self.file)
            elif self.type == 'complex':
                c_rsf.sf_complexread(np.reshape(data,(data.size)),self.file)
            elif self.type == 'int':
                c_rsf.sf_intread(np.reshape(data,(data.size,)),self.file)
            else:
                raise TypeError('Unsupported file type %s' % self.type)
        else:
            if self.type == 'float':
                self.file.floatread(np.reshape(data,(data.size,)))
            elif self.type == 'int':
                self.file.intread(np.reshape(data,(data.size,)))
            else:
                raise TypeError('Unsupported file type %s' % self.type)
        if allocate:
            return data
            
class Output(_File):
    def __init__(self,tag='out',data_format=None):
        self.file = None
        if _swig_:
            if python2:
                # c function only knows utf-8 (ascii).  translate the unicode
                tag = tag.encode('utf-8')
            self.file = c_rsf.sf_output(tag)
        else:
            self.file = _RSF(False,tag)

        _File.__init__(self,tag)
        if data_format:
                self.setformat(data_format)
    def write(self,data):
        if _swig_:
            if self.type == 'float':
                c_rsf.sf_floatwrite(np.reshape(data.astype(np.float32),(data.size,)),self.file)
            elif self.type == 'complex':
                c_rsf.sf_complexwrite(np.reshape(data,(data.size,)),
                                      self.file)
            elif self.type == 'int':
                c_rsf.sf_intwrite(np.reshape(data.astype(np.int32),(data.size,)),self.file)
            elif self.type == 'uchar':
                c_rsf.sf_ucharwrite(np.reshape(data.astype(np.uint8),(data.size,)),self.file)
            else:
                raise TypeError('Unsupported file type %s' % self.type)
        else:
            if self.type == 'float':
                self.file.floatwrite(np.reshape(data.astype(np.float32),(data.size,)))
            elif self.type == 'int':
                self.file.intwrite(np.reshape(data.astype(np.int32),(data.size,)))
            elif self.type == 'uchar':
                self.file.ucharwrite(np.reshape(data.astype(np.uint8),(data.size,)))
            else:
                raise TypeError('Unsupported file type %s' % self.type)
                
dataserver = os.environ.get('RSF_DATASERVER',
                            'http://www.reproducibility.org')

def Fetch(directory,filename,server=dataserver,top='data'):
    'retrieve a file from remote server'
    if server == 'local':
        remote = os.path.join(top,
                            directory,os.path.basename(filename))
        try:
            os.symlink(remote,filename)
        except:
            print ('Could not link file "%s" ' % remote)
            os.unlink(filename)
    else:
        rdir =  os.path.join(server,top,
                             directory,os.path.basename(filename))
        try:
            urllib.urlretrieve(rdir,filename)
        except:
            try:
                urllib.request.urlretrieve(rdir,filename)
            except:
                print ('Could not retrieve file "%s" from "%s"' % (filename,rdir))

class Filter(object):
    'Madagascar filter'
    plots = ('grey','contour','graph','contour3',
             'dots','graph3','thplot','wiggle','grey3')
    diagnostic = ('attr','disfil','headerattr')
    def __init__(self,name,prefix='sf',stdin=True,stdout=True):
        self.plot = False
        self.stdin = stdin
        self.stdout = stdout
        self.checkpar = False
        self.prog = None
        rsfroot = rsf.prog.RSFROOT
        if rsfroot:
            lp = len(prefix)
            if name[:lp] != prefix:
                name = prefix+name
            self.prog = rsf.doc.progs.get(name)
            prog = os.path.join(rsfroot,'bin',name)
            if os.path.isfile(prog):
                self.plot   = name[lp:] in Filter.plots
                if self.stdout and name[lp:] in Filter.diagnostic:
                    self.stdout = False
                name = prog
        self.command = name
        if self.prog:
            # self documentation
            self.__doc__ =  self.prog.text(None)
    def getdoc():
        '''for IPython'''
        return self.__doc__
    def _sage_argspec_():
        '''for Sage'''
        return None
    def __wrapped__():
        '''for IPython'''
        return None
    def __str__(self):
        'convert to string'
        return self.command
    def __or__(self,other):
        'pipe overload'
        new = self.copy()
        new.command = '%s | %s' % (self,other)
        new.stdout = other.stdout
        return new
    def setcommand(self,kw,args=[]):
        'set parameters'
        parstr = []
        for (key,val) in kw.items():
            if key[:2] == '__': # convention to handle -- parameters
                key = '--'+key[2:]
            if self.checkpar and self.prog and not self.prog.pars.get(key):
                sys.stderr.write('checkpar: No %s= parameter in %s\n' %
                                 (key,self.prog.name))
            if isinstance(val,str):
                val = '\''+val+'\''
            elif isinstance(val,File):
                val = '\'%s\'' % val
            elif isinstance(val,bool):
                if val:
                    val = 'y'
                else:
                    val = 'n'
            elif isinstance(val,list):
                val = ','.join(map(str,val))
            else:
                val = str(val)
            parstr.append('='.join([key,val]))
        self.command = ' '.join([self.command,
                                 ' '.join(map(str,args)),
                                 ' '.join(parstr)])
    def apply(self,*srcs):
        'Apply to data'

        # crude handle of stdin=0, fix later
        first = srcs[0]
        if isinstance(first,np.ndarray) or isinstance(first,File):
            mysrcs = list(srcs)
        else:
            mysrcs = []

        # handle numpy input
        numpy = False
        for n, src in enumerate(mysrcs):
            if isinstance(src,np.ndarray):
                mysrcs[n] = File(src)
                numpy = True

        if self.stdout:
            if isinstance(self.stdout,str):
                out = self.stdout
            else:
                out = Temp()
            command = '%s > %s' % (self.command,out)
        else:
            command = self.command

        (first,pipe,second) = command.partition('|')
        if mysrcs:
            if self.stdin:
                command = ' '.join(['< ',str(mysrcs[0]),first] +
                                   [str(x) for x in mysrcs[1:]] +
                                   [pipe,second])
            else:
                command = ' '.join([first]+
                                   [str(x) for x in mysrcs] +
                                   [pipe,second])

        fail = os.system(command)
        if fail:
            raise RuntimeError('Could not run "%s" ' % command)

        if self.stdout:
            if self.plot:
                return Vplot(out,temp=True)
            else:
                outfile = File(out,temp=True)
                if numpy:
                    return outfile[:]
                else:
                    return outfile
    def __getitem__(self,srcs):
        'overload square brackets'
        return self.apply(srcs)
    def __call__(self,*args,**kw):
        'brackets set parameters'
        if args:
            self.stdout = args[0]
            if len(args) > 1:
                self.stdin = args[1]
        self.setcommand(kw)
        return self
    def pipe(self,other):
        new = self.copy()
    def __getattr__(self,attr):
        'overload dot: piping'
        new = Filter(attr)
        new.command = '%s | %s' % (self,new)
        new.stdin = self.stdin
        return new

def Vppen(plots,args):
    'combining plots'
    name = Temp()
    os.system('vppen %s %s > %s' % (args,
                                    ' '.join([str(plot) for plot in plots]),
                                    name))
    return Vplot(name,temp=True)

def Overlay(*plots):
    return Vppen(plots,'erase=o vpstyle=n')

def Movie(*plots):
    return Vppen(plots,'vpstyle=n')

def SideBySide(*plots,**kw):
    n = len(plots)
    iso = kw.get('iso')
    if iso:
        return Vppen(plots,'size=r vpstyle=n gridnum=%d,1' % n)
    else:
        return Vppen(plots,'yscale=%d vpstyle=n gridnum=%d,1' % (n,n))

def OverUnder(*plots,**kw):
    n = len(plots)
    iso = kw.get('iso')
    if iso:
        return Vppen(plots,'size=r vpstyle=n gridnum=1,%d' % n)
    else:
        return Vppen(plots,'xscale=%d vpstyle=n gridnum=1,%d' % (n,n))

class Vplot(object):
    def __init__(self,name,temp=False,penopts=''):
        'Constructor'
        self.name = name
        self.temp = temp
        self.img = None
        self.penopts = penopts+' '
    def __del__(self):
        'Destructor'
        if self.temp:
            try:
                os.unlink(self.name)
            except:
                raise RuntimeError('Could not remove "%s" ' % self)
    def __str__(self):
        return self.name
    def __mul__(self,other):
        'Overload multiplicaton'
        return Overlay(self,other)
    def __add__(self,other):
        'Overload addition'
        return Movie(self,other)
    def show(self):
        'Show on screen'
        os.system('sfpen %s' % self.name)
    def hard(self,printer='printer'):
        'Send to printer'
        os.system('PRINTER=%s pspen %s' % (printer,self.name))
    def image(self):
        'Convert to PNG in the current directory (for use with IPython and SAGE)'
        self.img = os.path.basename(self.name)+'.png'
        self.export(self.img,'png',args='bgcolor=w')
    def _repr_png_(self):
        'return PNG representation'
        if not self.img:
            self.image()
        img = open(self.img,'rb')
        guts = img.read()
        img.close()
        return guts

    try:
        from IPython.display import Image

        @property
        def png(self):
            return Image(self._repr_png_(), embed=True)
    except:
        pass

    def movie(self):
        'Convert to animated GIF in the current directory (for use with SAGE)'
        self.gif = os.path.basename(self.name)+'.gif'
        self.export(self.gif,'gif',args='bgcolor=w')
    def export(self,name,format=None,pen=None,args=''):
        'Export to different formats'
        from rsf.vpconvert import convert
        if not format:
            if len(name) > 3:
                format = name[-3:].lower()
            else:
                format = 'vpl'
        convert(self.name,name,format,pen,self.penopts+args,verb=False)

class _Wrap(object):
    'helper class to wrap all Madagascar programs as Filter objects'
    def __init__(self, wrapped):
        self.wrapped = wrapped
    def __getattr__(self, name):
        try:
            return getattr(self.wrapped, name)
        except AttributeError:
            if name in rsf.doc.progs.keys() or \
              'sf'+name in rsf.doc.progs.keys():
                return Filter(name)
            else:
                raise

sys.modules[__name__] = _Wrap(sys.modules[__name__])

little_endian = (sys.byteorder == 'little')
    
class _RSF(object):
    'rsf file (emulates api/c/file.c)'
    infiles = [None]
    def __init__(self,input,tag=None):
        global par
        self.pars = _Simtab()
        for type in ('int','float','bool'):
            setattr(self,type, getattr(self.pars,'get'+type))
            setattr(self,type+'s',getattr(self.pars,'get'+type+'s'))
        if input: # input file
            # set data stream
            if tag==None or tag=='in':
                self.stream = sys.stdin
                filename = None
            else:
                filename = par.string(tag)
                if filename==None:
                    filename = tag
                self.stream = open(filename,'r+')
            # temporary file for header
            self.headname = Temp() 
            self.head = open(self.headname,'w+')
            # read parameters
            self.pars.input(self.stream,self.head)
            # keep track of input files
            global infiles
            if filename==None:
                _RSF.infiles[0] = self
            else:
                _RSF.infiles.append(self)
            # get dataname
            filename = self.string('in')
            if filename==None:
                raise TypeError('No in= in file "%s" ' % tag)
            self.dataname = filename
            # keep stream in the special case of in=stdin
            if filename != 'stdin':
                self.stream = open(filename,'r+b')
                
            # set format
            data_format = self.string('data_format')
            if not data_format:
                data_format = 'ascii_float'
            self.setformat(data_format)
        else: # output file
            if tag==None or tag=='out':
                self.stream = sys.stdout
                headname = None
            else:
                headname = par.string(tag)
                if headname==None:
                    headname = tag
                self.stream = open(headname,'w+')
            self.headname = None
            # check if piping
            try:
                t = self.stream.tell()
                self.pipe = False
            except:
                self.pipe = True
            if self.stream == sys.stdout:
                dataname = par.string('--out') or par.string('out')
            else:
                dataname = None
            if self.pipe:
                self.dataname = 'stdout'
            elif dataname==None:
                path = rsf.path.datapath()
                name = self.getfilename()
                if name != None:
                    if name=='/dev/null':
                        self.dataname = 'stdout'
                    else:
                        self.dataname = os.path.join(path,name+'@')
                else:
                    # invent a name
                    self.dataname = Temp()
            else:
                self.dataname = dataname

            self.putstring('in',self.dataname)

            if None != _RSF.infiles[0]:
                data_format = _RSF.infiles[0].string('data_format','native_float')
            else:
                data_format = 'native_float' 
                
            self.setformat(data_format)
            
            self.rw = par.bool('--readwrite',False)
            self.dryrun = par.bool('--dryrun',False)                    
    def getfilename(self):
        'find the name of the file to which we are writing'
        found_stdout = False

        f = '/dev/null'
        if os.fstat(1)[1] == os.stat(f)[1]:
            found_stdout = True
        else:        
            for f in os.listdir('.'):
                # Comparing the unique file ID stored by the OS for the file stream
                # stdout with the known entries in the file table:
                if os.path.isfile(f) and os.fstat(1)[1] == os.stat(f)[1]:
                    found_stdout = True
                    break

        if found_stdout:
            return f
        else:
            return None
    def gettype(self):
        return self.type
    def getform(self):
        return self.form
    def settype(self,type):
        self.type = type
    def setform(self,form):
        self.form = form
        if form == 'ascii':
            if None != self.dataname:
                self.put('esize',0) # for compatibility with SEPlib
            self.aformat = None
            self.eformat = None
            self.aline = 8
    def setformat(self,dataformat):
        done = False
        for type in ('float','int','complex','uchar','short','long','double'):
            if type in dataformat:
                self.settype(type)
                done = True
                break
        if not done:
            if 'byte' in dataformat:
                self.settype('uchar')
            else:
                self.settype('char')
        if dataformat[:6]=='ascii_':
            self.setform('ascii')
        elif dataformat[:4]=='xdr_':
            self.setform('xdr')
        else:
            self.setform('native')
    def string(self,key,default=None):
        get = self.pars.getstring(key)
        if get:
            return get
        else:
            return default
    def fileflush(self,src):
        import pwd, socket, time
        if None==self.dataname:
            return
        if None != src and None != src.head:
            src.head.seek(0)
            for line in src.head:
                self.stream.write(line)

        user = os.getuid()
        username = pwd.getpwuid(user)[0]
        self.stream.write("%s\t%s:\t%s@%s\t%s\n" % \
                       (par.getprog(),
                         os.getcwd(),
                        username,
                        socket.gethostname(),
                        time.ctime()))
        self.putstring('data_format','_'.join([self.form,self.type]))
        self.pars.output(self.stream)
        self.stream.flush()

        if self.dataname == 'stdout':
            # keep stream, write the header end code
            self.stream.write('\tin="stdin"\n\n\x0c\x0c\x04')
            self.stream.flush()
        else:                   
            self.stream = open(self.dataname,'w+b')

        self.dataname = None
        if self.dryrun:
            sys.exit(0)
    def fflush(self):
        self.stream.flush()
    def putint(self,key,par):
        if None==self.dataname:
            raise TypeError('putint to a closed file')
        val = '%d' % par
        self.pars.enter(key,val)
    def putints(self,key,par,n):
        if None==self.dataname:
            raise TypeError('putints to a closed file')
        val = ''
        for i in range(n-1):
            val += '%d,' % par[i]
        val += '%d' % par[n-1]
        self.pars.enter(key,val)
    def putfloat(self,key,par):
        if None==self.dataname:
            raise TypeError('putfloat to a closed file')
        val = '%g' % par
        self.pars.enter(key,val)
    def putfloats(self,key,par,n):
        if None==self.dataname:
            raise TypeError('putfloats to a closed file')
        val = ''
        for i in range(n-1):
            val += '%g,' % par[i]
        val += '%g' % par[n-1]
        self.pars.enter(key,val)
    def putstring(self,key,par):
        if None==self.dataname:
            raise TypeError('putstring to a closed file')
        val = '\"%s\"' % par
        self.pars.enter(key,val)
    def ucharwrite(self,arr):
        if None != self.dataname:
            self.fileflush(_RSF.infiles[0])
                
        try:
            self.stream.buffer.write(arr.tobytes())
        except:
            if python2:
                self.stream.write(arr.tostring())
            else:
                self.stream.write(arr.tobytes())
    def intwrite(self,arr):
        if None != self.dataname:
            self.fileflush(_RSF.infiles[0])
                
        if self.form=='ascii':
            if self.aformat == None:
                aformat = '%d '
            else:
                aformat = self.aformat
            if self.eformat == None:
                eformat = '%d '
            else:
                eformat = self.aformat
            size = arr.size
            farr = arr.flatten()
            left = size    
            while left > 0:
                if self.aline < left:
                    nbuf = self.aline
                else:
                    nbuf = left
                last = size-left+nbuf-1
                for i in range(size-left,last):
                    self.stream.write(aformat % farr[i])
                self.stream.write(eformat % farr[last])
                self.stream.write("\n")
                left -= nbuf
        else:
            try:
                self.stream.buffer.write(arr.tobytes())
            except:
                if python2:
                    self.stream.write(arr.tostring())
                else:
                    self.stream.write(arr.tobytes())
    def intread(self,arr):
        if self.form=='ascii':
            arr[:] = np.loadtxt(self.stream,dtype='int32',count=arr.size)
        else:
            try:
                data = self.stream.buffer.read(arr.size*4)
            except:
                data = self.stream.read(arr.size*4)
                if not python2 and type(data) == str:
                    data = data.encode()
            arr[:] = np.frombuffer(data,dtype='int32')
    def floatwrite(self,arr):
        if None != self.dataname:
            self.fileflush(_RSF.infiles[0])
                
        if self.form=='ascii':
            if self.aformat == None:
                aformat = '%g '
            else:
                aformat = self.aformat
            if self.eformat == None:
                eformat = '%g '
            else:
                eformat = self.aformat
            size = arr.size
            farr = arr.flatten()
            left = size    
            while left > 0:
                if self.aline < left:
                    nbuf = self.aline
                else:
                    nbuf = left
                last = size-left+nbuf-1
                for i in range(size-left,last):
                    self.stream.write(aformat % farr[i])
                self.stream.write(eformat % farr[last])
                self.stream.write("\n")
                left -= nbuf
        else:
            try:
                self.stream.buffer.write(arr.tobytes())
            except:
                if python2:
                    self.stream.write(arr.tostring())
                else:
                    self.stream.write(arr.tobytes())
    def floatread(self,arr):
        if self.form=='ascii':
            arr[:] = np.loadtxt(self.stream,dtype='float32',count=arr.size)
        else:
            try:
                data = self.stream.buffer.read(arr.size*4)
            except:
                data = self.stream.read(arr.size*4)
                if not python2 and type(data) == str:
                    data = data.encode()
            arr[:] = np.frombuffer(data,dtype='float32')
    def tell(self):
        return self.stream.tell()
    def bytes(self):
        if self.dataname=='stdin':
            return -1
        if self.dataname==None:
            st = os.fstat(self.stream.fileno())
        else:
            st = os.stat(self.dataname)
        return st.st_size
    def fileclose(self):
        if self.stream != sys.stdin and \
          self.stream != sys.stdout and \
          self.stream != None:
            self.stream.close()
            self.stream = None
        if self.headname != None:
            os.unlink(self.headname)
            self.headname = None


if __name__ == "__main__":

#      a=100 Xa=5
#      float=5.625 cc=fgsg
#      dd=1,2x4.0,2.25 true=yes false=2*no label="Time (sec)"
   
#    no_swig()
    # Testing getpar
    par = Par(["prog","a=5","b=as","a=100","float=5.625",
               "true=y","false=n"]) #,"par=%s" % sys.argv[0]])
    assert 100 == par.int("a")
    assert not par.int("c")
    assert 10 == par.int("c",10)
    assert 5.625 == par.float("float")
    assert par.bool("true")
    assert not par.bool("false")
    #assert "Time (sec)" == par.string("label")
    #assert "Time (sec)" == par.string("label","Depth")
    assert not par.string("nolabel")
    assert "Depth" == par.string("nolabel","Depth")
    # no function for this   par.close()
    # Testing file
    # Redirect input and output
    inp = os.popen("sfspike n1=100 d1=0.25 nsp=2 k1=1,10 label1='Time'")
    out = open("junk.rsf","w")
    os.dup2(inp.fileno(),sys.stdin.fileno())
    os.dup2(out.fileno(),sys.stdout.fileno())
    # Initialize
    par = Par()
    input = Input()
    output = Output()
    # Test
    assert 'float' == input.type
    assert 'native' == input.form
    n1 = input.int("n1")
    assert 100 == n1
    assert 0.25 == input.float("d1")
    assert 'Time' == input.string("label1")
    n2 = 10
    output.put('n2',n2)
    output.put('label2','Distance (kft)')
    trace = np.zeros(n1,'f')
    input.read(trace)
    for i in range(n2):
        output.write(trace)
    os.system("sfrm junk.rsf")
