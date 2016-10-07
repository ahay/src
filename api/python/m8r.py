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
import os, sys, tempfile, re, subprocess, urllib

try:
    import c_m8r as c_rsf
    import numpy
    _swig_ = True
except:
    _swig_ = False

import rsf.doc
import rsf.prog
import rsf.path

###

def view(name):
    try:
        from IPython.display import Image
        png = name+'.png'
        makefile = os.path.join(rsf.prog.RSFROOT,'include','Makefile')
        os.system('make -f %s %s' % (makefile,png))
        return Image(filename=png)
    except:
        print 'No IPython Image support'
        return None

if _swig_:
    class Par(object):
        '''parameter table'''
        def __init__(self,argv=sys.argv):
            c_rsf.sf_init(len(argv),argv)
            self.prog = c_rsf.sf_getprog()
            for type in ('int','float','bool'):
                setattr(self,type,self.__get(type))
                setattr(self,type+'s',self.__gets(type))
        def __get(self,type):
            func = getattr(c_rsf,'sf_get'+type)
            def _get(key,default=None):
                get,par = func(key)
                if get:
                    return par
                elif default != None:
                    return default
                else:
                    return None
            return _get
        def __gets(self,type):
            func = getattr(c_rsf,'get'+type+'s')
            def _gets(key,num,default=None):
                pars = func(key,num)
                if pars:
                    return pars
                elif default:
                    return default
                else:
                    return None
            return _gets
        def string(self,key,default=None):
            par = c_rsf.sf_getstring(key)
            if par:
                return par
            elif default:
                return default
            else:
                return None
else:
    # from apibak
    class Par(object):
        '''parameter table'''
        def __init__(self,argv=sys.argv):
            self.noArrays = True
            self.prog = argv[0]
            self.__args = self.__argvlist2dict(argv[1:])
        def __argvlist2dict(self,argv):
            """Eliminates duplicates in argv and makes it a dictionary"""
            argv = self.__filter_equal_sign(argv)
            args = {}
            for a in argv:
                key = a.split('=')[0]
                args[key] = a.replace(key+'=','')
            return args

        def __filter_equal_sign(self,argv):
            """Eliminates "par = val", "par= val" and "par =val" mistakes."""
            argv2 = []
            # Could not use the simpler 'for elem in argv'...argv.remove because
            # lonely '=' signs are treated weirdly. Many things did not work as
            # expected -- hence long and ugly code. Test everything.
            for i in range( len(argv) ):
                if argv[i] != '=':
                    if argv[i].find('=') != 0:
                        if argv[i].find('=') != -1:
                            if argv[i].find('=') != len(argv[i])-1:
                                argv2.append(argv[i])
            return argv2

        def __get(self, key, default):
            """Obtains value of argument from dictionary"""
            if self.__args.has_key(key):
                return self.__args[key]
            elif str(default):
                return default
            else:
                return None

        def string(self, key, default=None):
            """Returns string argument given to program"""
            return self.__get(key, default)

        def int(self,key,default=None):
            """Returns integer argument given to program"""
            try:
                return int( self.__get(key, default) )
            except:
                return None

        def float(self,key,default=None):
            """Returns float argument given to program"""
            try:
                return float( self.__get(key, default) )
            except:
                return None

        def bool(self,key,default=None):
            """Returns bool argument given to program"""
            val = self.__get(key, default)
            val = str(val).lower()
            if val[0] == 'y' or val == 'true':
                return True
            elif val[0] =='n' or val == 'false':
                return False
            else:
                return None

# default parameters for interactive runs
par = Par(['python','-'])

class Temp(str):
    'Temporaty file name'
    datapath = rsf.path.datapath()
    tmpdatapath = os.environ.get('TMPDATAPATH',datapath)
    def __new__(cls):
        return str.__new__(cls,tempfile.mktemp(dir=Temp.tmpdatapath))

class File(object):
    attrs = ['rms','mean','norm','var','std','max','min','nonzero','samples']
    def __init__(self,tag,temp=False,name=''):
        'Constructor'
        if isinstance(tag,File):
            # copy file (name is ignored)
            self.__init__(tag.tag)
            tag.close()
        elif _swig_ and isinstance(tag,numpy.ndarray):
            # numpy array
            if not name:
                name = Temp()
            out = Output(name)
            shape = tag.shape
            dims = len(shape)
            for axis in range(1,dims+1):
                out.put('n%d' % axis,shape[dims-axis])
            out.write(tag)
            out.close()
            self.__init__(out,temp=True)
        elif _swig_ and isinstance(tag,list):
            self.__init__(numpy.array(tag,'f'),name)
        else:
            self.tag = tag
        self.temp = temp
        self.narray = None
        for filt in Filter.plots + Filter.diagnostic:
            # run things like file.grey() or file.sfin()
            setattr(self,filt,Filter(filt,srcs=[self],run=True))
        for attr in File.attrs:
            setattr(self,attr,self.want(attr))
    def __str__(self):
        'String representation'
        if self.tag:
            tag = str(self.tag)
            if os.path.isfile(tag):
                return tag
            else:
                raise TypeError, 'Cannot find "%s" ' % tag
        else:
            raise TypeError, 'Cannot find tag'
    def sfin(self):
        'Output of sfin'
        return Filter('in',run=True)(0,self)
    def want(self,attr):
        'Attributes from sfattr'
        def wantattr():
            try:
                val = os.popen('%s want=%s < %s' % 
                               (Filter('attr'),attr,self)).read()
            except:
                raise RuntimeError, 'trouble running sfattr'
            m = re.search('=\s*(\S+)',val)
            if m:
                val = float(m.group(1))
            else:
                raise RuntimeError, 'no match'
            return val
        return wantattr
    def real(self):
        'Take real part'
        re = Filter('real')
        return re[self]
    def cmplx(self,im):
        c = Filter('cmplx')
        return c[self,im]
    def imag(self):
        'Take imaginary part'
        im = Filter('imag')
        return im[self]
    def __add__(self,other):
        'Overload addition'
        add = Filter('add')
        return add[self,other]
    def __sub__(self,other):
        'Overload subtraction'
        sub = Filter('add')(scale=[1,-1])
        return sub[self,other]
    def __mul__(self,other):
        'Overload multiplication'
        try:
            mul = Filter('scale')(dscale=float(other))
            return mul[self]
        except:
            mul = Filter('mul')(mode='product')
            return mul[self,other]
    def __div__(self,other):
        'Overload division'
        try:
            div = Filter('scale')(dscale=1.0/float(other))
            return div[self]
        except:
            div = Filter('add')(mode='divide')
            return div[self,other]
    def __neg__(self):
        neg = Filter('scale')(dscale=-1.0) 
        return neg[self]
    def dot(self,other):
        'Dot product'
        # incorrect for complex numbers
        prod = self.__mul__(other)
        stack = Filter('stack')(norm=False,axis=0)[prod]
        return stack[0]
    def cdot2(self):
        'Dot product with itself'
        abs2 = Filter('math')(output="abs(input)").real[self]
        return abs2.dot(abs2)
    def dot2(self):
        'Dot product with itself'
        return self.dot(self)
    def __array__(self,context=None):
        'numpy array'
        if _swig_:
            if None == self.narray:
                if not hasattr(self,'file'):
                    f = c_rsf.sf_input(self.tag)
                else:
                    f = self.file
                self.narray = c_rsf.rsf_array(f)
                if not hasattr(self,'file'):
                    c_rsf.sf_fileclose(f)
            return self.narray
        else:
            # gets only the real part of complex arrays
            val = os.popen('%s < %s' % 
                           (Filter('disfil')(number=False),self)).read()

            return map(lambda x: float(x.rstrip(',')),val.split())
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
        if _swig_:
            if hasattr(self,'file'):
                f = self.file
            else:
                f = c_rsf.sf_input(self.tag)
            s = c_rsf.sf_leftsize(f,dim)
            if not hasattr(self,'file'):
                c_rsf.sf_fileclose(f)
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
        try:
            p = subprocess.Popen('%s %s parform=n < %s' % 
                                 (Filter('get'),key,self),
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 close_fds=True)
            get = p.stdout.read()
        except:
            raise RuntimeError, 'trouble running sfget'
        if get:
            val = int(get)
        elif default:
            val = default
        else:
            val = None
        return val
    def shape(self):
        # axes are reversed for consistency with numpy
        s = []
        dim = 1
        for i in range(1,10):
            ni = self.int('n%d' % i,1)
            if ni > 1:
                dim = i
            s.append(ni)
        s = s[:dim]
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
        return put[self]
    def __len__(self):
        return self.size()
    def close(self):
        if self.temp:
            Filter('rm',run=True)(0,self)
    def __del__(self):
        print 'Closing File'
        self.close()

if _swig_:
    class _File(File):
        type = ['uchar','char','int','float','complex']
        form = ['ascii','xdr','native']
        def __init__(self,tag):
            if not self.file:
                raise TypeError, 'Use Input or Output instead of File'
            File.__init__(self,tag)
            self.type = _File.type[c_rsf.sf_gettype(self.file)]
            self.form = _File.form[c_rsf.sf_getform(self.file)]
        def tell(self):
            return c_rsf.sf_tell(self.file)
        def close(self):
            c_rsf.sf_fileclose(self.file)
        def __del__(self):
            print 'Closing ', self.file
            self.close()
            File.close(self)
        def settype(self,type):
            for i in xrange(len(_File.type)):
                if type == _File.type[i]:
                    self.type = type
                    c_rsf.sf_settype (self.file,i)
        def setformat(self,format):
            c_rsf.sf_setformat(self.file,format)
        def __get(self,func,key,default):
            get,par = func(self.file,key)
            if get:
                return par
            elif default:
                return default
            else:
                return None
        def __gets(self,func,key,num,default):
            pars = func(self.file,key,num)
            if pars:
                return pars
            elif default:
                return default
            else:
                return None
        def string(self,key):
            return c_rsf.sf_histstring(self.file,key)
        def int(self,key,default=None):
            return self.__get(c_rsf.sf_histint,key,default)
        def float(self,key,default=None):
            return self.__get(c_rsf.sf_histfloat,key,default)
        def ints(self,key,num,default=None):
            return self.__gets(c_rsf.histints,key,num,default)    
        def bytes(self):
            return c_rsf.sf_bytes(self.file)
        def put(self,key,val):
            if isinstance(val,int):
                c_rsf.sf_putint(self.file,key,val)
            elif isinstance(val,float):
                c_rsf.sf_putfloat(self.file,key,val)
            elif isinstance(val,str):
                c_rsf.sf_putstring(self.file,key,val)
            elif isinstance(val,list):
                if isinstance(val[0],int):
                    c_rsf.sf_putints(self.file,key,val)
        
    class Input(_File):
        def __init__(self,tag='in'):
            if isinstance(tag,File):
                # copy file
                self.__init__(tag.tag)
                self.copy = True
            else:
                self.file = c_rsf.sf_input(tag)
                _File.__init__(self,tag)
                self.copy = False
        def read(self,data):
            if self.type == 'float':
                c_rsf.sf_floatread(numpy.reshape(data,(data.size,)),self.file)
            elif self.type == 'complex':
                c_rsf.sf_complexread(numpy.reshape(data,(data.size)),self.file)
            else:
                raise TypeError, 'Unsupported file type %s' % self.type
        def close(self):
            if not self.copy:
                c_rsf.sf_fileclose(self.file)
            _File.close(self)

    class Output(_File):
        def __init__(self,tag='out',src=None):
            if not tag:
                self.tag = Temp()
                self.temp = True
            else:
                self.tag = tag
                self.temp = False
            self.file = c_rsf.sf_output(self.tag)
            if src: # clone source file
                if hasattr(src,'file'):
                   srcfile = src.file
                   srctype = src.type
                else:
                   srcfile = c_rsf.sf_input(self.tag)
                   srctype = c_rsf.sf_gettype(srcfile)
                c_rsf.sf_settype(self.file,_File.type.index(srctype))
                c_rsf.sf_fileflush(self.file,srcfile)
                if not hasattr(src,'file'):
                    c_rsf.sf_fileclose(srcfile)
            _File.__init__(self,self.tag)
        def write(self,data):
            if self.type == 'float':
                c_rsf.sf_floatwrite(numpy.reshape(data.astype(numpy.float32),(data.size,)),self.file)
            elif self.type == 'complex':
                c_rsf.sf_complexwrite(numpy.reshape(data,(data.size,)),
                                      self.file)
            elif self.type == 'int':
                c_rsf.sf_intwrite(numpy.reshape(data.astype(numpy.int32),(data.size,)),self.file)
            else:
                raise TypeError, 'Unsupported file type %s' % self.type
        def close(self):
            c_rsf.sf_fileclose(self.file)
            _File.close(self)
            
else:

    class Input(object):
        def __create_variable_dictionary(self, header):
            'Parse RSF header into a dictionary of variables'
            self.vd={} # variable dictionary
            ilist = header.split()
            pos = 0
            squot = "'"
            dquot = '"'
            while pos < len(ilist):
                if '=' in ilist[pos]:
                    tokenlist = ilist[pos].split('=')
                    lhs = tokenlist[0]
                    rhs = tokenlist[1]
                    quotmark = None
                    if rhs[0] in (squot, dquot):
                        if rhs[0] == squot:
                            quotmark = squot
                        else:
                            quotmark = dquot
                        if rhs[-1] == quotmark:
                            rhs_out = rhs.strip(quotmark)
                            pos += 1
                        else:
                            rhs_out = rhs.lstrip(quotmark)
                            while pos < len(ilist):
                                pos += 1
                                rhs_out += ' '
                                if ilist[pos][-1] == quotmark:
                                    rhs_out += ilist[pos][:-1]
                                    break
                                else:
                                    rhs_out += ilist[pos]
                    else:
                        rhs_out = rhs
                        pos += 1
                    self.vd[lhs] = rhs_out
                else:
                    pos += 1

        def __init__(self,tag='in'):
            # Temporary solution. Need to scan for \EOL\EOL\EOT, else this will
            # choke on a .HH file!
            if tag == 'in':
                self.__create_variable_dictionary(sys.stdin.read())
            else:
                try:
                    f = open(str(tag),'r')
                    self.__create_variable_dictionary(f.read())
                    f.close()
                except:
                    sys.stderr.write("Cannot read from \"%s\"\n" % tag)
                    sys.exit(1)

        def int(self, nm):
            return int(self.vd[nm])

        def float(self, nm):
            return float(self.vd[nm])

    class Output(object):
        def __init__(self,tag='out',src=None):
            print 'To use Output, you need to install SWIG http://www.swig.org/'
            sys.exit(1)

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
            print 'Could not link file "%s" ' % remote
            os.unlink(filename)
    else:
        rdir =  os.path.join(server,top,
                             directory,os.path.basename(filename))
        try:
            urllib.urlretrieve(rdir,filename)
        except:
            print 'Could not retrieve file "%s" from "%s"' % (filename,rdir)
        
class Filter(object):
    'Madagascar filter'
    plots = ('grey','contour','graph','contour3',
             'dots','graph3','thplot','wiggle','grey3')
    diagnostic = ('attr','disfil')
    def __init__(self,name,prefix='sf',srcs=[],
                 run=False,checkpar=False,pipe=False):
        rsfroot = rsf.prog.RSFROOT
        self.plot = False
        self.stdout = True
        self.prog = None
        if rsfroot:
            lp = len(prefix)
            if name[:lp] != prefix:
                name = prefix+name
            self.prog = rsf.doc.progs.get(name)   
            prog = os.path.join(rsfroot,'bin',name)
            if os.path.isfile(prog):
                self.plot   = name[lp:] in Filter.plots
                self.stdout = name[lp:] not in Filter.diagnostic
                name = prog
        self.srcs = srcs
        self.run=run
        self.command = name
        self.checkpar = checkpar
        self.pipe = pipe
        if self.prog:
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
        return self.command
    def __or__(self,other):
        'pipe overload'
        self.command = '%s | %s' % (self,other) 
        return self
    def setcommand(self,kw,args=[]):
        parstr = []
        for (key,val) in kw.items():
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
    def __getitem__(self,srcs):
        'Apply to data'
        mysrcs = self.srcs[:]
        if isinstance(srcs,tuple):
            mysrcs.extend(srcs)
        elif srcs:
            mysrcs.append(srcs)

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
            command = ' '.join(['< ',str(mysrcs[0]),first]+
                               map(str,mysrcs[1:])+[pipe,second])  
                
        fail = os.system(command)
        if fail:
            raise RuntimeError, 'Could not run "%s" ' % command

        if self.stdout:
            if self.plot:
                return Vplot(out,temp=True)
            else:
                return File(out,temp=True)
    def __call__(self,*args,**kw):
        if args:
            self.stdout = args[0]
            self.run = True
        elif not kw and not self.pipe:
            self.run = True
        self.setcommand(kw,args[1:])
        if self.run:
            return self[0]
        else:
            return self
    def __getattr__(self,attr):
        'Making pipes'
        other = Filter(attr)
        self.pipe = True
        self.command = '%s | %s' % (self,other)
        return self

def Vppen(plots,args):
    name = Temp()
    os.system('vppen %s %s > %s' % (args,' '.join(map(str,plots)),name))
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
                raise RuntimeError, 'Could not remove "%s" ' % self
    def __str__(self):
        return self.name
    def __mul__(self,other):
        return Overlay(self,other)
    def __add__(self,other):
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
     def __init__(self, wrapped):
         self.wrapped = wrapped
     def __getattr__(self, name):
         try:
             return getattr(self.wrapped, name)
         except AttributeError:
             if name in rsf.doc.progs.keys() or 'sf'+name in rsf.doc.progs.keys():
                 return Filter(name)
             else:
                 raise

sys.modules[__name__] = _Wrap(sys.modules[__name__])


if __name__ == "__main__":
    import numpy

#      a=100 Xa=5
#      float=5.625 cc=fgsg
#      dd=1,2x4.0,2.25 true=yes false=2*no label="Time (sec)"
    
    # Testing getpar
    par = Par(["prog","a=5","b=as","a=100","par=%s" % sys.argv[0]])
    assert 100 == par.int("a")
    assert not par.int("c")
    assert 10 == par.int("c",10)
    assert 5.625 == par.float("float")
    assert [1.0, 4.0, 4.0, 2.25] == par.floats("dd",4)
    assert par.bool("true")
    no = par.bools("false",2)
    assert no and not no[0] and not no[1]
    assert "Time (sec)" == par.string("label")
    assert "Time (sec)" == par.string("label","Depth")
    assert not par.string("nolabel")
    assert "Depth" == par.string("nolabel","Depth")
    par.close()
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
    assert 10 == output.int('n2')
    output.put('label2','Distance (kft)')
    input.put("n",[100,100])
    assert [100,100] == input.ints("n",2)
    trace = numpy.zeros(n1,'f')
    input.read(trace)
    for i in xrange(n2):
        output.write(trace)
    os.system("sfrm junk.rsf")
