import os, sys, types
import numarray
import c_rsf

#      a=100 Xa=5
#      float=5.625 cc=fgsg
#      dd=1,2x4.0,2.25 true=yes false=2*no label="Time (sec)"

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

datatype = ('char','int','float','complex')
dataform = ('ascii','xdr','native')

class File:
    def __init__(self):
        self.type = datatype[c_rsf.sf_gettype(self.file)]
        self.form = dataform[c_rsf.sf_getform(self.file)]
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
    def __init__(self,tag="in"):
        self.file = c_rsf.sf_input(tag)
        File.__init__(self)
    def read(self,data,n=-1):
        if n<0:
            n = numarray.size(data)
        c_rsf.sf_read(data,n,self.file)

class Output(File):
    def __init__(self,tag="out"):
        self.file = c_rsf.sf_output(tag)
        File.__init__(self)
    def write(self,data,n=-1):
        if n<0:
            n = numarray.size(data)       
        c_rsf.sf_write(data,n,self.file)

if __name__ == "__main__":
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
    trace = numarray.zeros(n1,'f')
    input.read(trace)
    for i in xrange(n2):
        output.write(trace)
    os.system("sfrm junk.rsf")
    
# 	$Id: rsf.py,v 1.4 2004/03/22 05:43:24 fomels Exp $	
