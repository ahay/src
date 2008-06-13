import os, sys, types, tempfile, re
import c_rsf
import numpy

#      a=100 Xa=5
#      float=5.625 cc=fgsg
#      dd=1,2x4.0,2.25 true=yes false=2*no label="Time (sec)"

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

class Par:
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
par = Par(['self','-'])

datatype = ['uchar','char','int','float','complex']
dataform = ['ascii','xdr','native']

class File(object):
    def __init__(self):
        self.type = datatype[c_rsf.sf_gettype(self.file)]
        self.form = dataform[c_rsf.sf_getform(self.file)]
    def __add__(self,other):
        'Overload addition'
        add = Output(None)
        if os.path.isfile(self.tag) and os.path.isfile(other.tag):
            if os.system('sfadd %s %s > %s' % (self.tag,other.tag,add.tag)):
                raise TypeError, 'Could not run sfadd'
        else:
            raise TypeError, 'Cannot add %s and %s' % (self.tag,other.tag)
        return add
    def __mul__(self,other):
        'Overload multiplication'
        mul = Output(None)
        if isinstance(other,File):
            if os.path.isfile(self.tag) and os.path.isfile(other.tag):
                if os.system('sfadd mode=p %s %s > %s' % (self.tag,other.tag,mul.tag)):
                    raise TypeError, 'Could not run sfadd'
            else:
                raise TypeError, 'Cannot multiply %s and %s' % (self.tag,other.tag)
        else:
            try:
                fmul = float(other)
            except:
                raise TypeError, 'Cannot cast to float'
            if os.path.isfile(self.tag):
                if os.system('sfscale < %s dscale=%g > %s' % (self.tag,fmul,mul.tag)):
                    raise TypeError, 'Could not run sfscale'
            else:
                raise TypeError, 'No file "%s" ' % self.tag
        return mul
    def size(self,dim=0):
        return c_rsf.sf_leftsize(self.file,dim)
    def settype(self,type):
        for i in xrange(len(datatype)):
            if type == datatype[i]:
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
    def fileclose(self):
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
    def __init__(self,tag='in'):
        if isinstance(tag,File):
            self.__init__(tag.tag)
        else:
            self.tag = tag
            self.file = c_rsf.sf_input(self.tag)
            File.__init__(self)
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
            self.tag = tempfile.mktemp(dir=tmpdatapath)
        else:
            self.tag = tag
        self.file = c_rsf.sf_output(self.tag)
        if src: # clone source file
            c_rsf.sf_settype(self.file,datatype.index(src.type))
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
    def __init__(self,name,prefix='sf'):
        self.prog = name
        rsfroot = os.environ.get('RSFROOT')
        if rsfroot:
            lp = len(prefix)
            if name[:lp] != prefix:
                name = prefix+name
            prog = os.path.join(rsfroot,'bin',name)
            if os.path.isfile(prog):
                self.prog = prog
    def __call__(self,inp=None,pars={}):
        out = Output(None)
        params = ' '.join([key+'='+str(val) for (key,val) in pars.items()])
        if inp:
            if os.isfile(inp.tag):
                command = '< %s %s %s > %s' % (inp.tag,self.prog,params,out.tag)
            else:
                raise TypeError, 'No file "%s" ' % inp.tag
        else:
            command = '%s %s > %s' % (self.prog,params,out.tag)
        if os.system(command):
            raise TypeError, 'Could not run %s' % self.prog
        return out                                   

if __name__ == "__main__":
    import numpy
    
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
