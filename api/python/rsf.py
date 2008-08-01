import os, sys, types, tempfile, re
import c_rsf
import numpy

class Par(object):
    def __init__(self,argv=sys.argv):
        c_rsf.sf_init(len(argv),argv)
        self.prog = c_rsf.sf_getprog()
    def close(self):
        c_rsf.sf_parclose()
    def __get(self,func,key,default):
        get,par = func(key)
        if get:
            return par
        elif default:
            return default
        else:
            return None
    def __gets(self,func,key,num,default):
        pars = func(key,num)
        if pars:
            return pars
        elif default:
            return default
        else:
            return None
    def string(self,key):
        return c_rsf.sf_getstring(key)
    def int(self,key,default=None):
        return self.__get(c_rsf.sf_getint,key,default)
    def float(self,key,default=None):
        return self.__get(c_rsf.sf_getfloat,key,default)
    def bool(self,key,default=None):
        return self.__get(c_rsf.sf_getbool,key,default)
    def ints(self,key,num,default=None):
        return self.__gets(c_rsf.getints,key,num,default)
    def floats(self,key,num,default=None):
        return self.__gets(c_rsf.getfloats,key,num,default)
    def bools(self,key,num,default=None):
        return self.__gets(c_rsf.getbools,key,num,default)

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
                pathfile = open(os.path.join(os.environ.get('HOME'),'.datapath'),'r')
            except:
                pathfile = None
        if pathfile:
            for line in pathfile.readlines():
                check = re.match("(?:%s\s+)?datapath=(\S+)" % os.uname()[1],line)
                if check:
                    datapath = check.group(1)
            pathfile.close()
        if not datapath:
            datapath = './' # the ultimate fallback
    tmpdatapath = os.environ.get('TMPDATAPATH',datapath)
    def __new__(cls):
        return str.__new__(cls,tempfile.mktemp(dir=Temp.tmpdatapath))

class File(object):
    type = ['uchar','char','int','float','complex']
    form = ['ascii','xdr','native']
    def __init__(self):
        'Constructor'
        if not self.file:
            raise TypeError, 'Use Input or Output instead of File'
        self.type = File.type[c_rsf.sf_gettype(self.file)]
        self.form = File.form[c_rsf.sf_getform(self.file)]
        self.narray = None
        for plot in Filter.plots:
            setattr(self,plot,Filter(plot,inp=self)) 
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
    def __del__(self):
        'Destructor'
        if self.temp:
            os.system('sfrm %s' % self)
    def sfin(self):
        'Output of sfin'
        p = os.popen('sfin %s' % self)
        p_str = p.read()
        p.close()
        print p_str
    def attr(self):
        'Output of sfattr'
        p = os.popen('sfattr < %s' % self)
        p_str = p.read()
        p.close()
        print p_str
    def __add__(self,other):
        'Overload addition'
        add = Temp()
        if os.system('sfadd %s %s > %s' % (self,other,add)):
            raise TypeError, 'Could not run sfadd'
        return Input(add,temp=True)
    def __mul__(self,other):
        'Overload multiplication'
        mul = Temp()
        if isinstance(other,File):
            if os.system('sfadd mode=p %s %s > %s' % 
                         (self,other,mul)):
                raise TypeError, 'Could not run sfadd'
        else:
            try:
                fmul = float(other)
            except:
                raise TypeError, 'Cannot cast to float'

            if os.system('sfscale < %s dscale=%g > %s' % 
                         (self,fmul,mul)):
                raise TypeError, 'Could not run sfscale'
        return Input(mul,temp=True)
    def __div__(self,other):
        'Overload division'
        mul = Temp()
        if isinstance(other,File):
            if os.system('sfadd mode=d %s %s > %s' % 
                         (self,other,mul)):
                raise TypeError, 'Could not run sfadd'
        else:
            try:
                fmul = float(other)
            except:
                raise TypeError, 'Cannot cast to float'
            if os.system('sfscale < %s dscale=%g > %s' % 
                         (self,1.0/fmul,mul)):
                raise TypeError, 'Could not run sfscale'
        return Input(mul,temp=True)
    def __array__(self,context=None):
        if None == self.narray:
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
        return c_rsf.sf_leftsize(self.file,dim)
    def __len__(self):
        return self.size()
    def settype(self,type):
        for i in xrange(len(File.type)):
            if type == File.type[i]:
                self.type = type
                c_rsf.sf_settype (self.file,i)
    def setformat(self,format):
        c_rsf.sf_setformat(self.file,format)
        File.__init__(self)
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
    def close(self):
        return c_rsf.sf_fileclose(self.file)
    def put(self,key,val):
        what = type(val)
        if what is types.IntType:
            c_rsf.sf_putint(self.file,key,val)
        elif what is types.FloatType:
            c_rsf.sf_putfloat(self.file,key,val)
        elif what is types.StringType:
            c_rsf.sf_putstring(self.file,key,val)
        elif what is types.ListType:
            what = type(val[0])
            if what is types.IntType:
                c_rsf.sf_putints(self.file,key,val)

class Input(File):
    def __init__(self,tag='in',temp=False):
        if isinstance(tag,File):
            # copy file
            self.__init__(tag.tag)
            self.old = tag # to increase refcount
        elif isinstance(tag,numpy.ndarray):
            # numpy array
            out = Output(None)
            shape = tag.shape
            for axis in range(len(shape)):
                out.put('n%d' % (axis+1),shape[axis])
            out.write(tag)
            out.close()
            self.__init__(out,temp=True)
        else:
            self.tag = tag
            self.file = c_rsf.sf_input(self.tag)
            File.__init__(self)
        self.temp=temp
    def read(self,data):
        if self.type == 'float':
            c_rsf.sf_floatread(numpy.reshape(data,(data.size,)),self.file)
        elif self.type == 'complex':
            c_rsf.sf_complexread(numpy.reshape(data,(data.size)),self.file)
        else:
            raise TypeError, 'Unsupported file type %s' % self.type

class Output(File):
    def __init__(self,tag='out',src=None):
        if not tag:
            self.tag = Temp()
            self.temp = True
        else:
            self.tag = tag
            self.temp = False
        self.file = c_rsf.sf_output(self.tag)
        if src: # clone source file
            c_rsf.sf_settype(self.file,File.type.index(src.type))
            c_rsf.sf_fileflush(self.file,src.file)
        File.__init__(self)
    def write(self,data):
        if self.type == 'float':
            c_rsf.sf_floatwrite(numpy.reshape(data,(data.size,)),self.file)
        elif self.type == 'complex':
            c_rsf.sf_complexwrite(numpy.reshape(data,(data.size,)),self.file)
        else:
            raise TypeError, 'Unsupported file type %s' % self.type

class Filter(object):
    plots = ('grey','contour','graph','contour3','dots','graph3','thplot','wiggle')
    def __init__(self,name,prefix='sf',inp=None):
        self.prog = name
        rsfroot = os.environ.get('RSFROOT')
        self.plot = False
        if rsfroot:
            lp = len(prefix)
            if name[:lp] != prefix:
                name = prefix+name
            prog = os.path.join(rsfroot,'bin',name)
            if os.path.isfile(prog):
                self.prog = prog
                self.plot = name[lp:] in Filter.plots
        self.inp = inp
    def __call__(self,*inp,**kw):
        out = Temp()
        pars = []
        for (key,val) in kw.items():
            if isinstance(val,str):
                val = '\''+val+'\''
            else:
                val = str(val)
            pars.append('='.join([key,val]))
        params = ' '.join(pars)
        if inp:
            command = '< %s %s %s > %s' % (inp[0],self.prog,params,out)
        elif self.inp:
            command = '< %s %s %s > %s' % (self.inp,self.prog,params,out)
        else:
            command = '%s %s > %s' % (self.prog,params,out)
        if os.system(command):
            raise TypeError, 'Could not run %s' % self.prog
        if self.plot:
            return Vplot(out,temp=True)
        else:
            return Input(out,temp=True)


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
    def export(self,name):
        'Export to different formats'
        if len(name) > 3:
            suffix = name[-3:].lower()
        else:
            suffix = 'vpl'
        if suffix in ('eps','gif','avi'):
            os.system('vplot2%s %s %s' % (suffix,self.name,name))
        else:
            os.system('cp %s %s' % (self.name,name))

class _Wrap(object):
     def __init__(self, wrapped):
         self.wrapped = wrapped
     def __getattr__(self, name):
         try:
             return getattr(self.wrapped, name)
         except AttributeError:
             return Filter(name)

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
