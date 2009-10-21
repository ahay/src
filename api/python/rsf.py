import os, sys, tempfile, re
import c_rsf
import numpy

import rsfdoc
import rsfprog

###

class Par(object):
    def __init__(self,argv=sys.argv):
        c_rsf.sf_init(len(argv),argv)
        self.prog = c_rsf.sf_getprog()
        for type in ('int','float','bool'):
            setattr(self,type,self.__get(type))
            setattr(self,type+'s',self.__gets(type))
    def close(self):
        c_rsf.sf_parclose()
    def __get(self,type):
        func = getattr(c_rsf,'sf_get'+type)
        def _get(key,default=None):
            get,par = func(key)
            if get:
                return par
            elif default:
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

# default parameters for interactive runs
par = Par(['python','-'])

class Temp(str):
    'Temporaty file name'
    datapath = os.environ.get('DATAPATH')
    if not datapath:
        try:
            pathfile = open('.datapath','r')
        except:
            try:
                pathfile = open(os.path.join(os.environ.get('HOME'),
                                             '.datapath'),'r')
            except:
                pathfile = None
        if pathfile:
            for line in pathfile.readlines():
                check = re.match("(?:%s\s+)?datapath=(\S+)" % os.uname()[1],
                                 line)
                if check:
                    datapath = check.group(1)
            pathfile.close()
        if not datapath:
            datapath = './' # the ultimate fallback
    tmpdatapath = os.environ.get('TMPDATAPATH',datapath)
    def __new__(cls):
        return str.__new__(cls,tempfile.mktemp(dir=Temp.tmpdatapath))

class File(object):
    attrs = ['rms','mean','norm','var','std','max','min','nonzero','samples']
    def __init__(self,tag,temp=False):
        'Constructor'
        if isinstance(tag,File):
            # copy file
            self.__init__(tag.tag)
        elif isinstance(tag,numpy.ndarray):
            # numpy array
            out = Output(Temp())
            shape = tag.shape
            dims = len(shape)
            for axis in range(1,dims+1):
                out.put('n%d' % axis,shape[dims-axis])
            out.write(tag)
            out.close()
            self.__init__(out,temp=True)
        elif isinstance(tag,list):
            self.__init__(numpy.array(tag,'f'))
        else:
            self.tag = tag
        self.temp = temp
        self.narray = None
        for filt in Filter.plots + Filter.diagnostic:
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
            mul = Filter('add')(mode='product')
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
        prod = self.__mul__(other).reshape()
        stack = Input(Filter('stack')(norm=False,axis=1)[prod])
        dp = numpy.zeros(1,'f')
        stack.read(dp)
        stack.close()
        return dp[0]
    def __array__(self,context=None):
        'numpy array'
        # danger: dangling open file descriptor
        if None == self.narray:
            if not hasattr(self,'file'):
                self.file = c_rsf.sf_input(self.tag)
            self.narray = c_rsf.rsf_array(self.file)
        return self.narray
    def __array_wrap__(self,array,context=None):
        return Input(array)
    def __getitem__(self,i):
        array = self.__array__()
        return array[i]
    def __setitem__(self,index,value):
        array = self.__array__()
        array.__setitem__(index,value)
    def size(self,dim=0):
        if hasattr(self,'file'):
            f = self.file
        else:
            f = c_rsf.sf_input(self.tag)
        s = c_rsf.sf_leftsize(f,dim)
        if not hasattr(self,'file'):
            c_rsf.sf_fileclose(f)
        return s
    def int(self,key,default=None):
        try:
            inp, out, err = os.popen3('%s %s parform=n < %s' % 
                                      (Filter('get'),key,self))
            get = out.read()
        except:
            raise RuntimeError, 'trouble running sfget'
        inp.close()
        out.close()
        err.close()
        if get:
            val = int(get)
        elif default:
            val = default
        else:
            val = None
        return val
    def shape(self):
        s = []
        dim = 1
        for i in range(1,10):
            ni = self.int('n%d' % i,1)
            if ni > 1:
                dim = i
            s.append(ni)
        return tuple(s[:dim])
    def reshape(self,shape=None):
        if not shape:
            shape = self.size()
        try:
            shape = tuple(shape)
        except:
            shape = (shape,)
        old = self.shape()
        lold = len(old)
        lshape = len(shape)
        puts = {}
        for i in range(max(lold,lshape)):
            ni = 'n%d' % (i+1)
            if i < lold and i < lshape:
                if old[i] != shape[i]:
                    puts[ni] = shape[i]
            elif i < lold:
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
        self.close()

class _File(File):
    type = ['uchar','char','int','float','complex']
    form = ['ascii','xdr','native']
    def __init__(self,tag):
        if not self.file:
            raise TypeError, 'Use Input or Output instead of File'
        File.__init__(self,tag)
        self.type = _File.type[c_rsf.sf_gettype(self.file)]
        self.form = _File.form[c_rsf.sf_getform(self.file)]
    def close(self):
        c_rsf.sf_fileclose(self.file)
    def __del__(self):
        self.close()
        File.close(self)
    def settype(self,type):
        for i in xrange(len(File.type)):
            if type == File.type[i]:
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
        else:
            self.file = c_rsf.sf_input(tag)
            _File.__init__(self,tag)
    def read(self,data):
        if self.type == 'float':
            c_rsf.sf_floatread(numpy.reshape(data,(data.size,)),self.file)
        elif self.type == 'complex':
            c_rsf.sf_complexread(numpy.reshape(data,(data.size)),self.file)
        else:
            raise TypeError, 'Unsupported file type %s' % self.type

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
            c_rsf.sf_settype(self.file,_File.type.index(src.type))
            c_rsf.sf_fileflush(self.file,src.file)
        _File.__init__(self,self.tag)
    def write(self,data):
        if self.type == 'float':
            c_rsf.sf_floatwrite(numpy.reshape(data,(data.size,)),self.file)
        elif self.type == 'complex':
            c_rsf.sf_complexwrite(numpy.reshape(data,(data.size,)),self.file)
        else:
            raise TypeError, 'Unsupported file type %s' % self.type

class Filter(object):
    'Madgagascar filter'
    plots = ('grey','contour','graph','contour3',
             'dots','graph3','thplot','wiggle')
    diagnostic = ('attr','disfil')
    def __init__(self,name,prefix='sf',srcs=[],run=False,checkpar=False):
        rsfroot = os.environ.get('RSFROOT')
        self.plot = False
        self.stdout = True
        self.prog = None
        if rsfroot:
            lp = len(prefix)
            if name[:lp] != prefix:
                name = prefix+name
            self.prog = rsfdoc.progs.get(name)   
            prog = os.path.join(rsfroot,'bin',name)
            if os.path.isfile(prog):
                self.plot   = name[lp:] in Filter.plots
                self.stdout = name[lp:] not in Filter.diagnostic
                name = prog
        self.srcs = srcs
        self.run=run
        self.command = name
        self.checkpar = checkpar
        if self.prog:
            self.__doc__ =  self.prog.docstring()
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
            
        if mysrcs:    
            command = '< %s %s %s' % \
                (mysrcs[0],command,' '.join(map(str,mysrcs[1:])))  
                
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
        elif not kw:
            self.run = True
        self.setcommand(kw,args[1:])
        if self.run:
            return self[0]
        else:
            return self
    def __getattr__(self,attr):
        'Making pipes'
        other = Filter(attr)
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
    def __init__(self,name,temp=False):
        'Constructor'
        self.name = name
        self.temp = temp
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
        os.system('xtpen %s' % self.name)
    def hard(self,printer='printer'):
        'Send to printer'
        os.system('PRINTER=%s pspen %s' % (printer,self.name))
    def image(self):
        'Convert to PNG in the current directory (for use with SAGE)'
        self.png = os.path.basename(self.name)+'.png'
        self.export(self.png,'png')
    def export(self,name,format=None):
        'Export to different formats'
        if not format:
            if len(name) > 3:
                format = name[-3:].lower()
            else:
                format = 'vpl'
        if format in ('eps','gif','avi','png'):
            os.system('vplot2%s %s %s' % (format,self.name,name))
        else:
            os.system('cp %s %s' % (self.name,name))

class _Wrap(object):
     def __init__(self, wrapped):
         self.wrapped = wrapped
     def __getattr__(self, name):
         try:
             return getattr(self.wrapped, name)
         except AttributeError:
             if name in rsfdoc.progs.keys() or 'sf'+name in rsfdoc.progs.keys():
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
#    output.settype('int')
#    assert 'int' == output.type
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
    
# 	$Id: rsf.py 3148 2007-11-13 00:13:20Z sfomel $	
