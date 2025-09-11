import matplotlib.ticker as ticker
import matplotlib.colorbar as colorbar
from datetime import datetime
import os, re, sys, io
from subprocess import Popen,PIPE, SubprocessError
from subprocess import run as Run
import numpy as np
import matplotlib.pyplot as plt
from SCons.Util import WhereIs
import copy
from matplotlib.font_manager import fontManager
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.ticker import MaxNLocator as MaxNLocator1
from mpl_toolkits.axisartist import Axes
from mpl_toolkits.axisartist.floating_axes import \
    GridHelperCurveLinear
from matplotlib import transforms
from matplotlib import colors as mcolors
import matplotlib as mpl
from matplotlib import colorbar as mbar

# fontManager.addfont("./times.ttf")
# fontManager.addfont("./timesbd.ttf")
# fontManager.addfont("./timesbi.ttf")
# fontManager.addfont("./timesi.ttf")
# fontManager.addfont("./arial.ttf")
# FONTDIR=os.path.dirname(__file__)
# fontManager.addfont(FONTDIR+"/NimbusRomNo9L-Reg.otf")
# fontManager.addfont(FONTDIR+"/NimbusRomNo9L-RegIta.otf")
# fontManager.addfont(FONTDIR+"/NimbusRomNo9L-Med.otf")
# fontManager.addfont(FONTDIR+"/NimbusRomNo9L-MedIta.otf")

# plt.rcParams['font.family'] = 'arial'
# plt.rcParams['font.family'] = 'Nimbus Roman'
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 15
from numpy import ndarray

# Check RSFROOT
RSFROOT = os.getenv("RSFROOT")
if RSFROOT is None:
    raise RuntimeWarning("RSFROOT environment variable not set.")
# if WhereIs("sfin") is None or WhereIs("sfpwd") is None:
#     raise RuntimeError("Madagascar not installed correctly. Check https://www.reproducibility.org/wiki/Installation for more information.")

# Check userid
try :user_sign = os.getlogin()
except: user_sign = os.getenv("USER")
try:
    import socket
    user_sign += "@"+socket.gethostname()
except:user_sign += "@localhost"

# RSF formats
RSFHSPLITER = b"\x0c\x0c\x04"

DATATOLSIZE=1000*1000


rsfforms = ["native", "ascii", "xdr"]
rsftypes = ["short", "int", "long", "float", "complex"]
forms = {"short": np.int16,
         "int": np.int32,
         "long": np.int64,
         "float": np.float32,
         "complex": np.complex64
             }
# alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890'

HANDLED_FUNCTIONS = {}

def implements(numpy_function):
    """Register an __array_function__ implementation for MyArray objects."""
    def decorator(func):
        HANDLED_FUNCTIONS[numpy_function] = func
        return func
    return decorator

class Rsfdata(np.ndarray):
    def __new__(cls,
                input_data=None,
                header=None, history="",
                dtype=np.float64,
                *args, **kwargs):
        '''
            Create the ndarray instance of rsf type, given the usual
            ndarray input arguments.  This will call the standard
            ndarray constructor, but return an object of rsf type.

            :param args: same as np.ndarray
            :param kwargs: same as np.ndarray
            '''
        if input_data is None:
            # obj = super().__new__(cls, *args, **kwargs)
            obj = np.zeros(*args, **kwargs).view(cls)
            obj.header, obj.history = {}, ""

        elif isinstance(input_data, io.IOBase):
            input_data.seek(0)
            obj = readrsf(input_data)

        elif isinstance(input_data, cls):
            obj = input_data.copy()

        elif type(input_data) == str:
            
            try:
                obj = readrsf(open(input_data, 'rb'))
            except:
                try:
                    suffix = kwargs.get("suffix", ".rsf")
                    obj = readrsf(open(input_data + suffix, 'rb'))
                except:
                    print("Error in reading file %s" % input_data, file=sys.stderr)
                    obj = np.zeros(*args, **kwargs).view(cls)
                    obj.header, obj.history = {}, ""
        else:
            obj = np.asarray(input_data, dtype=dtype,
                             *args, **kwargs).view(cls)

        if header is not None:
            obj.header.update(header)

        if history is not None:
            obj.history += history

        # print("__new__")
        return obj


    def __array_finalize__(self, obj):
        if obj is None: return
        if isinstance(obj, np.ndarray):
            try:
                self.header = copy.deepcopy(obj.header)
                self.history = obj.history
            except:
                self.header = {}
                self.history = ""
            self.update()

    def __array_function__(self, func, types, *args, **kwargs):
        '''
         Do not want to lose header and history after a few methods (like
         transpose). Almost the same as ndarray excepet a few attributes,
         so no need for NOT_IMPLEMENTS.
        '''
        if func not in HANDLED_FUNCTIONS:
            out = super(Rsfdata,self).__array_function__(func, types, *args, **kwargs)
            out = Rsfdata(out,self.header,self.history)
            out.history +="Python-%s Numpy-%s %s:\t%s\t%s \n\n" \
                     % (sys.version.split()[0],np.__version__,
                        func,
                        user_sign,
                        datetime.now().strftime('%a %b %d  %H:%M:%S %Y'))
            out.update()
            return out
        else: return HANDLED_FUNCTIONS[func](self, *args, **kwargs)




    def __array_wrap__(self, out_arr, context=None, return_scalar=False):

        # then just call the parent
        super(Rsfdata,self).__array_wrap__(out_arr, context, return_scalar)
        # print('In __array_wrap__:')
        # print('   self is %s' % repr(self))
        # print('   arr is %s' % repr(out_arr))
        if isinstance(out_arr, np.ndarray):
            out_arr = out_arr.view(Rsfdata)
            out_arr.header = self.header
            out_arr.history = self.history
            if context is not None:
                out_arr.history += "Python-%s Numpy-%s %s:\t%s\t%s \n\n" \
                               % (sys.version.split()[0], np.__version__,
                                  context[0],
                                  user_sign,
                                  datetime.now().strftime('%a %b %d  %H:%M:%S %Y'))
            out_arr.update()
        return out_arr



    def update(self):
        dim = len(self.data.shape)
        self.header.update({"dim":dim})
        for idim in range(dim):
            self.header.update({
                "n%d" % (idim + 1): self.data.shape[idim]
            })
            if "o%d" % (idim + 1) not in self.header.keys():
                self.header.update({"o%d" % (idim + 1): 0.})
            if "d%d" % (idim + 1) not in self.header.keys():
                self.header.update({"d%d" % (idim + 1): 1.})


    # def transpose(self, axes=None):
    #     out = super().transpose(axes)
    #     out = out.view(Rsfdata)
    #     out.header = self.header
    #     out.history = self.history
    #     out.update()
    #     return out
    # def fromfile(self,infile):
    #     '''
    #     read data from rsf file
    #     :param infile: TextIO | str
    #     :return: RsfData
    #     '''
    #     if type(infile) == str:
    #         infile = open(infile,'rb')
    #         return readrsf(infile)
    #     elif type(infile) == io.IOBase:
    #         return readrsf(infile)
    #     else:
    #         raise TypeError("Wrong type for infile %s" % type(infile))
    def transaxis(self,axis=None):
        return transaxis(self,axis)

    def write(self,outfile=None):
        '''
        write data to rsf file or BytesIO
        :param outfile: None | str
        :return: BytesIO | None
        '''
        if type(outfile) == str:
            writersf(self,outfile)
            return None
        return writersf(self)

    def axis(self, dict=None):
        '''
        Get axis parameters from headdict
        :param headdict: dict | None
        :return: List [Ndarray]
        '''
        axes = []
        dict1 = {}
        dict1.update(self.header)
        if dict is not None: dict1.update(dict)
        for idim in range(1, 10):
            if "n%d" % idim in dict1.keys():
                n = int(dict1["n%d" % idim])
                if "d%d" % idim in dict1.keys():
                    d = float(dict1["d%d" % idim])
                else:
                    d = 1
                if "o%d" % idim in dict1.keys():
                    o = float(dict1["o%d" % idim])
                else:
                    o = 0
                axes.append(np.linspace(o, o + n * d, n, False))
        return axes

    def put(self, headstr):
        self.header.update(sfput(headstr,self.header))

    def flow(self, cmd, **kargs):
        return flow(cmd,self.write(), **kargs)

    def copy(self):
        return copy.deepcopy(self)


@implements(np.transpose)
def transpose(*args, **kwargs):
    ...  # implementation of transpose for Rsfdata objects
    rdata = args[1][0]
    # print(type(rdata))
    axis = args[1][1]
    # print(axis)
    out = np.transpose(np.asarray(rdata),axis)
    out = Rsfdata(out,rdata.header,rdata.history)
    ndim = len(out.shape)
    # print(ndim)
    if axis is None:
        axis = list(range(ndim))
        axis = axis[::-1]
    for i in range(ndim):
        out.header["d%d"%(i+1)]=rdata.header["d%d"%(axis[i]+1)]
        out.header["o%d"%(i+1)]=rdata.header["o%d"%(axis[i]+1)]
        try: out.header["label%d"%(i+1)]=rdata.header["label%d"%(axis[i]+1)]
        except: pass
        try: out.header["unit%d"%(i+1)]=rdata.header["unit%d"%(axis[i]+1)]
        except: pass
    out.update()
    return out

def transaxis(inarr,axis=None):
    out = inarr
    ndim = len(out.shape)
    if axis is None:
        axis = list(range(ndim))
        axis = axis[::-1]
    for i in range(ndim):
        out.header["d%d"%(i+1)]=out.header["d%d"%(axis[i]+1)]
        out.header["o%d"%(i+1)]=out.header["o%d"%(axis[i]+1)]
        try: out.header["label%d"%(i+1)]=out.header["label%d"%(axis[i]+1)]
        except: pass
        try: out.header["unit%d"%(i+1)]=out.header["unit%d"%(axis[i]+1)]
        except: pass
    out.update()
    return out

def flow(cmd, source=None, dataout=True, verb=False):
    '''
    RSF Flow
    Usage:
      > sfin = flow("sfmath n1=1 output='1' | sfin")\n
      > print(sfin.read().decode())\n
        in:
            in="stdin"\n
            esize=4 type=float form=native\n
            n1=1           d1=1           o1=0\n
            1 elements 4 bytes\n
    :param source: str | BaseIO | None
        Input data file
    :param cmd: str
    :return: BytesIO
    '''
    out = io.BytesIO()

    if source is None:
        fsrc = None
    elif type(source) == str:
        fsrc = open(source,'rb')
    elif isinstance(source,io.IOBase):
        fsrc = source
    else: raise TypeError("Wrong type of source: %s"%type(source))

    if fsrc is not None:
        fsrc.seek(0)
        subprc = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=(None if verb else PIPE), shell=True)
        [sout, serr] = subprc.communicate(fsrc.read())
        if subprc.returncode != 0:
            raise SubprocessError("In Command: '%s':\n%s"%(cmd,serr.decode()))
    else :
        subprc = Run(cmd, stdin=None, stdout=PIPE, stderr=(None if verb else PIPE), shell=True, check=True)
        sout = subprc.stdout
        serr = subprc.stderr
        if subprc.returncode != 0:
            raise SubprocessError("In Command: '%s':\n%s"%(cmd,serr.decode()))

    if dataout:
        out.write(sout)
        out.seek(0)
        return Rsfdata(out)
    else: return sout

def sfput(fhead, headdict={}, wantdim=False):
    if type(fhead) is not str:
        raise TypeError("fhead must be a string")
    hdict = {}
    if headdict is not None: hdict.update(headdict)
    pattern = re.compile(r"\S+=\"{1}[^\"]+\"{1}|\S+=\'[.+]\'{1}|\S+=[\S]+|\S+=")
    fhead1 = pattern.findall(fhead)
    for record in fhead1:
        key0, val0 = record.split("=")
        if len(val0) > 1 and val0[0] in ["\"","'"] and val0[-1] in ["\"","'"] and val0[0]==val0[-1]:
            val0 = val0[1:-1]
        if val0 is None: val0 = ""
        hdict.update({key0:val0})
    dim = 1
    dims = [int(hdict["n1"])]
    for idim in range(2, 10):
        if "n%d" % idim in hdict.keys():
            dim += 1
            dims.append(int(hdict["n%d" % idim]))
    hdict["dim"] = dim
    if wantdim: return hdict, dims
    else:return hdict

def readrsf(fp):
    '''
    Read data and header from RSF data file

    :param fp:  TextIOWrapper | BytesIO
        file pointer to rsf
    :param rsfdata: RsfData | None
    :return: Tuple of (Ndarray, Dict)
        data and header
    '''
    fp.seek(0)

    fstream = fp.read()
    fheadend = fstream.find(RSFHSPLITER)

    fhead = str(fstream[:fheadend].decode())
    headdict,dims = sfput(fhead,wantdim=True)
    if headdict["in"] not in ["stdin","stdout"] :
        fpd = open(headdict["in"].replace("\'",""),'rb')
        datastream = fpd.read()
        fpd.close()
    else:
        datastream = fstream[fheadend + 3:]
    dim = headdict["dim"]
    while dim > 1 and dims[-1] == 1:
        dim -= 1
        dims.pop()
    if "data_format" in headdict.keys():
        format0 = headdict["data_format"]
    else:
        format0 = "native_float"
    fform, ftype = format0.split("_")
    if fform not in rsfforms:
        raise ValueError("Invalid format %s. Supported format: %s" % (fform, rsfforms))
    if ftype not in rsftypes:
        raise ValueError("Invalid type %s. Supported type: %s" % (ftype, rsftypes))
    if "xdr" == fform:
        try:
            fp1 = flow("sfdd form=native", fp)
        except Exception as e:
            raise ChildProcessError(str(e)+" Please manually use sfdd first.")
        fstream = fp1.read()
        fheadend = fstream.find(RSFHSPLITER)
    if "ascii"==fform:
        fdata = np.fromstring(datastream, dtype=forms[ftype]).copy()
    else:fdata = np.frombuffer(datastream, dtype=forms[ftype]).copy()
    dims.reverse()
    fdata = fdata.reshape(dims)
    fdata = fdata.transpose()
    return Rsfdata(input_data=fdata,
                   header=headdict,
                   history=str(fstream[:fheadend], encoding='utf-8').replace("\'","\""))


def strfun(str1):
    '''
    return str1.encode("utf-8")
    :param str1: str, to encode into utf-8
    :return:
        return str1.encode("utf-8")
    '''
    return str1.encode("utf-8")



def writersf(data, filename=None, dict=None, stdout=True):
    '''

    :param data: Rsfdata | np.ndarray[Any]
    :param filename: str
    :param dict: dict
        Overwrite data.header if type(data) == Rsfdata
    :return: out: io.BytesIO | None
        return None if filename is not None
    '''
    dims = []
    dict1 = {
        "o1":"0",
        "o2":"0",
        "d1":"1",
        "d2":"1",
        "data_format":"native_float"
    }
    outheader=""
    if type(data) == Rsfdata:
        if data.header: dict1.update(data.header)
    darray = data
    if dict: dict1.update(dict)
    if not stdout: dict1.update({"in":get_datapath() + filename + "@"})
    else:dict1.update({"in":"stdin"})
    # dtype check and transform
    try:
        fform,ftype = dict1["data_format"].split("_")
    except: fform,ftype = "native","float"
    if fform not in rsfforms:
        raise ValueError("Invalid format %s. Supported format: %s" % (fform, rsfforms))
    if ftype not in rsftypes:
        raise ValueError("Invalid type %s. Supported type: %s" % (ftype, rsftypes))
    darray = darray.astype(forms[ftype])
    for idim in range(len(darray.shape)):
        dims.append(darray.shape[idim])

    if filename == None:
        fout = io.BytesIO()
    else:
        fout = open(filename, "wb")
        if not stdout: fdata = open(get_datapath() + filename + "@", "wb")
    outheader += data.history
    try:
        oldheader = sfput(outheader)
    except:
        oldheader = {}
    if len(outheader)==0 or (not outheader[-1]=="\n"): outheader+="\n"
    outheader += "Python %s:\t%s\t%s \n\n\t" \
                 % (sys.version.split()[0],
                    user_sign,
                    datetime.now().strftime('%a %b %d  %H:%M:%S %Y'))

    for key,value in dict1.items():
        try: float(value)
        except: value = "\"%s\""%value
        outstr = "{key}={value}".format(key = key, value = value)
        if key in oldheader.keys():
            if dict1[key] != oldheader[key]: outheader += outstr + "\n\t"
            elif key in ["in","data_format","esize"]: outheader += outstr + "\n\t"
        else:
            outheader += outstr + "\n\t"
    if fform == "ascii":
        fmt = {"short": "%d",
               "int": "%d",
               "long": "%d",
               "float": "%.5f",
               "complex": "%.5f %.5fi",
               }[ftype]
        if stdout:
            np.savetxt(fout, darray.ravel(), fmt=fmt, delimiter=" ",
                   newline=" ", header=outheader + "\n" + RSFHSPLITER.decode(),
                   comments="")
        else:
            fout.write(strfun(outheader + "\n"))
            np.savetxt(fdata, darray.ravel(), fmt=fmt, delimiter=" ",
                       newline=" ", header="",
                       comments="")
    else:
        fout.write(strfun(outheader + "\n"))
        if fform == "xdr":
            print("Warning: Form \'xdr_%s\' not supported. Use \'native_%s\' instead. Please try sfdd form=xdr manually."%(ftype,ftype),file=sys.stderr)
            fout.write(strfun("\tdata_format=\"native_%s\"\n\n"%ftype))
        if stdout:
            fout.write(RSFHSPLITER)
            fout.write(darray.transpose().tobytes())
        else:
            fdata.write(darray.transpose().tobytes())

    if filename == None:
        fout.seek(0)
        return fout
    else:
        fout.close()
        if not stdout: fdata.close()
        return None


def get_datapath(cwd=os.getcwd()):

    top = datapath()
    path = os.path.dirname(top)
    if top[:2] != './':
        # create a hierarchical structure
        tree = dirtree(cwd)
        for level in tree:
            if level:
                path = os.path.join(path, level)
        mkdir(path)
    path = os.path.join(path, os.path.basename(top))
    return path

def datapath():
    '''
         Path for binary and temporary files.

           Datapath rules:
           1. check datapath= on the command line (Not applicable)
           2. check .datapath file in the current directory
           3. check DATAPATH environmental variable
           4. check .datapath in the home directory
           5. use '.' (not a SEPlib behavior)

       Note that option 5 results in:
       if you move the header to another directory then you must also move
       the data to the same directory. Confusing consequence.

        :return: str
           datapath
    '''
    path = os.environ.get('DATAPATH')
    if not path:
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
                check = re.match(r"(?:%s\s+)?datapath=(\S+)" % os.uname()[1],
                                 line)
                if check:
                    path = check.group(1)
            pathfile.close()
    if not path:
        path = './' # the ultimate fallback
    return path

def dirtree(cwd=None):
    'hierarcical directory structure'
    if not cwd:
        cwd = os.getcwd()
    tree = (os.path.basename(os.path.dirname(os.path.dirname(cwd))),
            os.path.basename(os.path.dirname(cwd)),
            os.path.basename(cwd))
    return tree

################################################################################

def mkdir(dir):
    'Recursive directory making'
    while os.path.basename(dir) == '.':
        dir = os.path.dirname(dir)
    if dir and not os.path.isdir(dir):
        mkdir(os.path.dirname(dir))
        os.mkdir(dir)
    return dir
 

def grey(rsfdata,
         vmin=None, vmax=None, pclip=100.,clip=None,allpos=False,
         color="seismic",
         figure=None,
         **kargs):

    if vmin is not None: vmin = float(vmin)
    if vmax is not None: vmax = float(vmax)
    if clip is not None: clip = float(clip)
    pclip = float(pclip)
    if allpos in ["True", "true", "T", "t", "1", "yes", "Yes", "YES", True]: allpos = True
    else: allpos = False

    n3 = 0 
    if "frame3" in kargs.keys(): n3 = int(kargs["frame3"])
    inputshape = rsfdata.shape
    if len(inputshape) > 2: 
        nn3 = rsfdata.size//inputshape[0]//inputshape[1]
        rsfdata = rsfdata.reshape(inputshape[0],inputshape[1],nn3)
        inputshape = inputshape[:2]
        try:rsfdata = rsfdata[:,:,n3-1]
        except: 
            rsfdata = rsfdata[:,:,np.newaxis]
            rsfdata = rsfdata[:,:,n3-1]
    if len(inputshape) == 1: raise IndexError("Need at least 2 dimensions. Input data size: %s"%str(inputshape))

    
    maxelem = np.max(np.abs(rsfdata.ravel()))
    if clip is not None:
        vmin = -clip
        vmax = clip
    if vmin == None and vmax==None:
        vmin = -maxelem
        vmax = maxelem
    elif vmin==None: vmin = -vmax
    elif vmax==None: vmax = -vmin
    if pclip>100 or pclip<=0: pclip=100.
    if pclip < 100.:
        vmax = np.percentile(np.abs(rsfdata),pclip)
        # vmax = np.percentile(rsfdata,pclip)
        # vmin = np.percentile(rsfdata,100.-pclip)
        vmin = -vmax
    if type(color) == str:
        try: usecmap = mpl.colormaps[color]
        except: usecmap = mpl.colormaps["gray"]
    elif isinstance(color,mcolors.Colormap):
        usecmap = color.copy()
    else: usecmap = mpl.colormaps["gray"]
    
    if allpos:
        usecmap = mcolors.LinearSegmentedColormap.from_list(
                 'allpos', usecmap(np.linspace(0.5, 1, 256))
                  )
        vmin=0.

    fhead0 = {
            "label1":"Time","label2":"Distance",
            "unit1":"s","unit2":"km",
            "title" : "",
            "n1":inputshape[0],
            "n2":inputshape[1],
            "d1":0.002,"o1":0.,
            "d2":1,"o2":0.,
            "fontsize":15,
            "fullsample": False,
            "colorbar": False,
            "bartitle": ""
        }
    fhead0.update(kargs)

    if not fhead0["fullsample"]:
        while rsfdata.size > DATATOLSIZE:
            if rsfdata.shape[0] > np.sqrt(DATATOLSIZE):
                rsfdata = rsfdata[::2,:]
                rsfdata.header["d1"] = float(rsfdata.header["d1"])*2
            if rsfdata.shape[1] > np.sqrt(DATATOLSIZE):
                rsfdata = rsfdata[:,::2]
                rsfdata.header["d2"] = float(rsfdata.header["d2"])*2

    if rsfdata.header is not None:
        fhead0.update(rsfdata.header)
    
    fhead0["label1"] = fhead0["label1"] + " (%s)" % fhead0["unit1"]
    fhead0["label2"] = fhead0["label2"] + " (%s)" % fhead0["unit2"]
    fhead0.update(kargs)
    plt.rcParams["font.size"] = fhead0["fontsize"]
    axes1 = rsfdata.axis()
    axist = axes1[0]
    axisx = axes1[1]

    if figure is not None:
        fig = figure
        axbase = fig.gca()
    else:
        fig = plt.gcf()
        axbase = plt.gca()
    axbase.axis("off")
    right_margin = 0.
    if fhead0["colorbar"]:
        right_margin = 0.15
        ax_bar = axbase.inset_axes([1.0-right_margin+0.05, 0, 0.05, 1])
    ax = axbase.inset_axes([0, 0, 1.0-right_margin, 1])
    
    im = ax.imshow(rsfdata.data, aspect="auto", cmap=usecmap, interpolation=kargs.get("interpolation", "nearest"),
              vmin=vmin, vmax=vmax, extent=[axisx[0], axisx[-1], axist[-1], axist[0]])

    ax.set_xlabel(fhead0["label2"],fontsize=fhead0["fontsize"])
    ax.set_ylabel(fhead0["label1"],fontsize=fhead0["fontsize"])
    ax.set_title(fhead0["title"],fontsize=fhead0["fontsize"]+5)
    lim1, lim2 = None, None
    
    if  "lim1" in fhead0.keys(): 
        if type(fhead0["lim1"]) is str: 
            lim1 = fhead0["lim1"][1:-1]
            lim1 = [float(lim1.split(",")[0]), float(lim1.split(",")[1])]
        else: lim1 = fhead0["lim1"]
        lim1 = lim1[::-1]
        ax.set_ylim(lim1)
    if  "lim2" in fhead0.keys(): 
        if type(fhead0["lim2"]) is str: 
            lim2 = fhead0["lim2"][1:-1]
            lim2 = [float(lim2.split(",")[0]), float(lim2.split(",")[1])]
        else: lim2 = fhead0["lim2"]
        ax.set_xlim(lim2)
 
    

    if fhead0["colorbar"]:
        if fhead0["colorbar"]:
            cbar = mbar.Colorbar(ax_bar, im, cmap=usecmap, label=fhead0["bartitle"])
            cbar.ax.get_yaxis().set_label_position('right')
    fig.add_axes(ax)
    fig.sca(ax)
    ax.set_label("Grey plot")
    return im,ax

def grey3flat(rsfdata, 
              frame1=0, frame2=0, frame3=0, point1=0.75, point2=0.5, 
              vmin=None, vmax=None, clip=None, pclip=100., allpos=False,
              color="seismic", 
              figure=None,
              **kargs):
    frame1 = int(frame1)
    frame2 = int(frame2)
    frame3 = int(frame3)
    point1 = float(point1)
    point2 = float(point2)
    if vmin is not None: vmin = float(vmin)
    if vmax is not None: vmax = float(vmax)
    if clip is not None: clip = float(clip)
    pclip = float(pclip)
    # if colorbar in ["True", "true", "T", "t", "1", "yes", "Yes", "YES"]: colorbar = True
    # else: colorbar = False
    # if wanttitle in ["True", "true", "T", "t", "1", "yes", "Yes", "YES"]: wanttitle = True
    # else: wanttitle = False


    # if rsfdata.header["dim"] < 3:
    #     raise IndexError("Need at least 3 dimensions. Input: %d"%rsfdata.header["dim"])

    inputshape = rsfdata.shape
    if len(inputshape) != 3: raise IndexError("Need 3 dimensions. Input data size: %s"%str(inputshape))
    
    if not isinstance(rsfdata, Rsfdata):
        rsfdata = Rsfdata(rsfdata, header=kargs)
    sfaxes = rsfdata.axis()
    frame1 = min(int(rsfdata.header["n1"])-1, frame1)
    frame2 = min(int(rsfdata.header["n2"])-1, frame2)
    frame3 = min(int(rsfdata.header["n3"])-1, frame3)
    slice1 = rsfdata[:,:,frame3]
    slice2 = rsfdata[:,frame2,:]
    slice3 = rsfdata[frame1,:,:]

    allslice = np.concatenate([slice1.ravel(), slice2.ravel(), slice3.ravel()],axis=0)
    # maxelems = [np.max(np.abs(slice1)),
    #             np.max(np.abs(slice2)),
    #             np.max(np.abs(slice3))]
    maxelem = np.max(np.abs(allslice))

    if clip is not None:
        vmin = -clip
        vmax = clip
    if vmin == None and vmax==None:
        vmin = -maxelem
        vmax = maxelem
    elif vmin==None: vmin = -vmax
    elif vmax==None: vmax = -vmin
    if pclip>100 or pclip<=0: pclip=100.
    if pclip < 100.:
        vmax = np.percentile(np.abs(allslice),pclip)
        # vmax = np.percentile(allslice,pclip)
        # vmin = np.percentile(allslice,100.-pclip)
        vmin = -vmax
    if type(color) == str:
        try: usecmap = mpl.colormaps[color]
        except: usecmap = mpl.colormaps["gray"]
    elif isinstance(color,mcolors.Colormap):
        usecmap = color.copy()
    else: usecmap = mpl.colormaps["gray"]
    
    # if allpos:
    #     usecmap = mcolors.LinearSegmentedColormap.from_list(
    #              'allpos', usecmap(np.linspace(0.5, 1, 256))
    #               )
    #     vmin=0.
    if allpos:
        if vmin < 0: vmin = 0.

    fhead0 = {
            "label1":"Time","label2":"Distance X","label3":"Distance X",
            "unit1":"s","unit2":"km","unit3":"km",
            "title" : "",
            "n1":inputshape[0],
            "n2":inputshape[1],
            "n3":inputshape[2],
            "d1":0.002,"o1":0.,
            "d2":1,"o2":0.,
            "d3":1,"o3":0.,
            "fontsize":15,
            "colorbar": False,
            "bartitle": "",
            "fullsample": False
        }
    
    if rsfdata.header is not None:
        fhead0.update(rsfdata.header)
    fhead0.update(kargs)

    nn1, nn2, nn3 = inputshape[0], inputshape[1], inputshape[2]
    j1, j2, j3 = 1, 1, 1
    if not fhead0["fullsample"]:
        while nn1*nn2*nn3 > DATATOLSIZE*np.sqrt(DATATOLSIZE):
            if nn1 > np.sqrt(DATATOLSIZE):
                nn1 = nn1//2
                j1 *= 2
            if nn2 > np.sqrt(DATATOLSIZE):
                nn2 = nn2//2
                j2 *= 2
            if nn3 > np.sqrt(DATATOLSIZE):
                nn3 = nn3//2
                j3 *= 2
    slice1 = slice1[::j1,::j2]
    slice2 = slice2[::j1,::j3]
    slice3 = slice3[::j2,::j3]
    sfaxes[0] = sfaxes[0][::j1]
    sfaxes[1] = sfaxes[1][::j2]
    sfaxes[2] = sfaxes[2][::j3]
    fhead0["d1"] = float(fhead0["d1"])*j1
    fhead0["d2"] = float(fhead0["d2"])*j2
    fhead0["d3"] = float(fhead0["d3"])*j3
    frame1 //= j1
    frame2 //= j2
    frame3 //= j3

    
    fhead0["label1"] = fhead0["label1"] + " (%s)" % fhead0["unit1"]
    fhead0["label2"] = fhead0["label2"] + " (%s)" % fhead0["unit2"]
    fhead0["label3"] = fhead0["label3"] + " (%s)" % fhead0["unit3"]
    fhead0.update(kargs)
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams["font.size"] = fhead0["fontsize"]

    # if label1.endswith("()"): label1 = label1[:-2]
    # if label2.endswith("()"): label2 = label2[:-2]
    # if label3.endswith("()"): label3 = label3[:-2]

    if color == 'gray':
        linecol = 'white'
        textcol = 'black'
    elif color == 'jet':
        linecol = 'black'
        textcol = 'black'
    else:
        linecol = kargs.get("linecol", "blue")
        textcol = kargs.get("textcol", "blue")

    fmt = ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-1, 3))
    current_map = usecmap

    if figure is not None:
        fig = figure
        ax0 = fig.gca()
    else:
        fig = plt.gcf()
        ax0 = plt.gca()
    ax0.set_label("Base axes")

    # width = fig.get_figwidth()
    # height = fig.get_figheight()

    btm_margin = 0.  # for labels and axis
    top_margin = 0.05  # for titles
    left_margin = 0  # for labels and axis
    right_margin = 0.2  # for labels and axis
    # ax0 = fig.add_axes([0, 0, 1, 1], facecolor='w')  # The base axes
    ax0.axis('off')
    ax0.set_frame_on(False)
    ax0.get_xaxis().set_visible(False)
    ax0.get_yaxis().set_visible(False)

    if fhead0["title"]!="":
        title_pos = [0, 1 - top_margin, 1, top_margin]
        ax_title = ax0.inset_axes(title_pos, facecolor='w')
        ax_title.set_frame_on(False)
        ax_title.get_xaxis().set_visible(False)
        ax_title.get_yaxis().set_visible(False)
        ax_title.set_xlim([0, 1])
        ax_title.set_ylim([0, 1])
        ax_title.text(0.5, 0.5, fhead0["title"], va='center', ha='center', fontsize=fhead0["fontsize"] + 5)
        # ax0.add_child_axes(ax_title)
    else:
        top_margin = 0.

    if fhead0["colorbar"]:
        # bar_pos = [1-right_margin,btm_margin,right_margin,1]
        # if newfig:
        #     bar_pos1 = [1 - right_margin + 0.1, btm_margin, right_margin - 0.2, 0.9 * (1 - top_margin - btm_margin)]
        #     ax_bar = fig.add_axes(bar_pos1, facecolor='none')
        # else:
        # if barpos == None: bar_pos1 = [0.9,btm_margin,right_margin - 0.2, 0.9 * (1 - top_margin - btm_margin)]
        # else :bar_pos1 = barpos
        bar_pos1 = [0.9,btm_margin,right_margin - 0.15, (1 - top_margin - btm_margin)]
        ax_bar = ax0.inset_axes(bar_pos1, facecolor='none')
        # ax_bar.set_title(fhead0["bartitle"])
    else: right_margin = 0.

    h_ax1 = 1 - top_margin - btm_margin
    w_ax1 = 1 - left_margin - right_margin

    ax1 = ax0.inset_axes([left_margin, btm_margin, w_ax1 * point2, h_ax1 * point1])
    ax0.add_child_axes(ax1)
    ax1.set_label("Axes 1 [n1, n2]")
    im0 = ax1.imshow(slice1, aspect='auto', cmap=current_map
                     , vmin=vmin, vmax=vmax
                     , extent=[sfaxes[1][0], sfaxes[1][-1], sfaxes[0][-1], sfaxes[0][0]])
    ax1.set_xlabel(fhead0["label2"])
    ax1.set_ylabel(fhead0["label1"])
    # ax1.set_yticklabels(ax1.get_yticklabels(), va='top')

    ax2 = ax0.inset_axes([left_margin + w_ax1 * point2, btm_margin, w_ax1 * (1 - point2), h_ax1 * point1], sharey=ax1)
    ax0.add_child_axes(ax2)
    ax2.set_label("Axes 2 [n2, n3]")
    im1=ax2.imshow(slice2, aspect='auto', cmap=current_map
               , vmin=vmin, vmax=vmax
               , extent=[sfaxes[2][0], sfaxes[2][-1], sfaxes[0][-1], sfaxes[0][0]])
    ax2.set_xlabel(fhead0["label3"])
    # ax2.set_xticklabels(ax2.get_xticklabels(), va='top', ha='center')

    ax2.yaxis.set_visible(False)
    ax2.xaxis.set_label_position('bottom')
    ax2.xaxis.set_ticks_position('bottom')

    ax3 = ax0.inset_axes([left_margin, btm_margin + h_ax1 * point1, w_ax1 * point2, h_ax1 * (1 - point1)], sharex=ax1)
    ax0.add_child_axes(ax3)
    ax3.set_label("Axes 3 [n1, n3]")
    im2=ax3.imshow(slice3.T, aspect='auto', cmap=current_map, origin="lower"
               , vmin=vmin, vmax=vmax
               , extent=[sfaxes[1][0], sfaxes[1][-1], sfaxes[2][0], sfaxes[2][-1]])
    ax3.xaxis.set_visible(False)
    # ax3.invert_yaxis()
    ax3.set_ylabel(fhead0["label3"])
    # ax3.set_yticklabels(ax3.get_yticklabels(), va='bottom')
    ax3.yaxis.set_label_position('left')


    if fhead0["colorbar"]:
        cbar = mbar.Colorbar(ax_bar, im1, cmap=current_map, format=fmt, label=fhead0["bartitle"])
        cbar.ax.get_yaxis().set_label_position('right')

    ax_01 = ax0.inset_axes([left_margin, btm_margin, w_ax1, h_ax1])
    ax_01.set_frame_on(False)
    ax_01.xaxis.set_visible(False)
    ax_01.yaxis.set_visible(False)
    ax_01.set_xlim([0, 1])
    ax_01.set_ylim([0, 1])

    # hline for frame1
    lframe1 = (len(sfaxes[0]) - frame1) / len(sfaxes[0]) * point1
    # lframe1 = sfaxes[0][frame1]
    ax1.axhline(sfaxes[0][frame1], linestyle='--', linewidth=1, color=linecol )
    ax2.axhline(sfaxes[0][frame1], linestyle='--', linewidth=1, color=linecol )
    ax2.text(1.01, lframe1, "%3.2f" % (sfaxes[0][frame1]), va='center', ha='left', color=textcol,transform=ax_01.transAxes, rotation=-90)
    # vline for frame2
    lframe2 = frame2 / len(sfaxes[1]) * point2
    ax1.axvline(sfaxes[1][frame2], linestyle='--', linewidth=1, color=linecol)
    ax3.axvline(sfaxes[1][frame2], linestyle='--', linewidth=1, color=linecol)
    ax3.text(lframe2, 1.01, "%3.2f" % (sfaxes[1][frame2]), va='bottom', ha='center', color=textcol,transform=ax_01.transAxes)
    # hline for frame3
    lframe3 = point1 + frame3 / len(sfaxes[2]) * (1 - point1)
    ax3.axhline(sfaxes[2][frame3], linestyle='--', linewidth=1, color=linecol )
    ax3.text(point2 + 0.01, lframe3, "%3.2f" % (sfaxes[2][frame3]), va='bottom', ha='left', color=textcol,transform=ax_01.transAxes)
    # vline for frame3
    lframe3 = point2 + frame3 / len(sfaxes[2]) * (1 - point2)
    # ax_01.plot([lframe3, lframe3], [0, point1], linestyle='--', linewidth=1, color=linecol)
    ax2.axvline(sfaxes[2][frame3], linestyle='--', linewidth=1, color=linecol)
    if frame3 / len(sfaxes[2]) >= 0.3:
        ax2.text(lframe3, point1 + 0.01, "%3.2f" % (sfaxes[2][frame3]), va='bottom', ha='left', color=textcol,transform=ax_01.transAxes)
    
    # fig.add_axes(ax_01)
    fig.add_axes(ax1)
    fig.add_axes(ax2)
    fig.add_axes(ax3)

    return ([ax1, ax2, ax3],
            [im0, im1, im2])

def grey3noflat(rsfdata, 
              frame1=0, frame2=0, frame3=0, point1=0.75, point2=0.5, 
              vmin=None, vmax=None, clip=None, pclip=100., allpos=False,
              color="seismic", 
              figure=None,
              **kargs):
    frame1 = int(frame1)
    frame2 = int(frame2)
    frame3 = int(frame3)
    point1 = float(point1)
    point2 = float(point2)
    if vmin is not None: vmin = float(vmin)
    if vmax is not None: vmax = float(vmax)
    if clip is not None: clip = float(clip)
    pclip = float(pclip)
    # if colorbar in ["True", "true", "T", "t", "1", "yes", "Yes", "YES"]: colorbar = True
    # else: colorbar = False
    # if wanttitle in ["True", "true", "T", "t", "1", "yes", "Yes", "YES"]: wanttitle = True
    # else: wanttitle = False


    # if rsfdata.header["dim"] < 3:
    #     raise IndexError("Need at least 3 dimensions. Input: %d"%rsfdata.header["dim"])

    inputshape = rsfdata.shape
    if len(inputshape) != 3: raise IndexError("Need 3 dimensions. Input data size: %s"%str(inputshape))
    if not isinstance(rsfdata, Rsfdata):
        rsfdata = Rsfdata(rsfdata, header=kargs)
    sfaxes = rsfdata.axis()
    frame1 = min(int(rsfdata.header["n1"])-1, frame1)
    frame2 = min(int(rsfdata.header["n2"])-1, frame2)
    frame3 = min(int(rsfdata.header["n3"])-1, frame3)
    slice1 = rsfdata[:,:,frame3]
    slice2 = rsfdata[:,frame2,:]
    slice3 = rsfdata[frame1,:,:].T

    allslice = np.concatenate([slice1.ravel(), slice2.ravel(), slice3.ravel()],axis=0)
    # maxelems = [np.max(np.abs(slice1)),
    #             np.max(np.abs(slice2)),
    #             np.max(np.abs(slice3))]
    maxelem = np.max(np.abs(allslice))

    if clip is not None:
        vmin = -clip
        vmax = clip
    if vmin == None and vmax==None:
        vmin = -maxelem
        vmax = maxelem
    elif vmin==None: vmin = -vmax
    elif vmax==None: vmax = -vmin
    if pclip>100 or pclip<=0: pclip=100.
    if pclip < 100.:
        vmax = np.percentile(np.abs(allslice),pclip)
        # vmin = np.percentile(allslice,100.-pclip)
        vmin = -vmax
    if type(color) == str:
        try: usecmap = mpl.colormaps[color]
        except: usecmap = mpl.colormaps["gray"]
    elif isinstance(color,mcolors.Colormap):
        usecmap = color.copy()
    else: usecmap = mpl.colormaps["gray"]
    
    # if allpos:
    #     usecmap = mcolors.LinearSegmentedColormap.from_list(
    #              'allpos', usecmap(np.linspace(0.5, 1, 256))
    #               )
    #     vmin=0.
    if allpos:
        if vmin < 0: vmin = 0.

    fhead0 = {
            "label1":"Time","label2":"Distance X","label3":"Distance X",
            "unit1":"s","unit2":"km","unit3":"km",
            "title" : "",
            "n1":inputshape[0],
            "n2":inputshape[1],
            "n3":inputshape[2],
            "d1":0.002,"o1":0.,
            "d2":1,"o2":0.,
            "d3":1,"o3":0.,
            "n1tic":5,
            "n2tic":5,
            "n3tic":5,
            "fontsize":15,
            "colorbar": False,
            "bartitle": "",
            "fullsample": False
        }
    
    if rsfdata.header is not None:
        fhead0.update(rsfdata.header)
    fhead0.update(kargs)

    nn1, nn2, nn3 = inputshape[0], inputshape[1], inputshape[2]
    j1, j2, j3 = 1, 1, 1
    if not fhead0["fullsample"]:
        while nn1*nn2*nn3 > DATATOLSIZE*np.sqrt(DATATOLSIZE):
            if nn1 > np.sqrt(DATATOLSIZE):
                nn1 = nn1//2
                j1 *= 2
            if nn2 > np.sqrt(DATATOLSIZE):
                nn2 = nn2//2
                j2 *= 2
            if nn3 > np.sqrt(DATATOLSIZE):
                nn3 = nn3//2
                j3 *= 2
    slice1 = slice1[::j1,::j2]
    slice2 = slice2[::j1,::j3]
    slice3 = slice3[::j2,::j3]
    sfaxes[0] = sfaxes[0][::j1]
    sfaxes[1] = sfaxes[1][::j2]
    sfaxes[2] = sfaxes[2][::j3]
    fhead0["d1"] = float(fhead0["d1"])*j1
    fhead0["d2"] = float(fhead0["d2"])*j2
    fhead0["d3"] = float(fhead0["d3"])*j3
    frame1 //= j1
    frame2 //= j2
    frame3 //= j3

    
    fhead0["label1"] = fhead0["label1"] + " (%s)" % fhead0["unit1"]
    fhead0["label2"] = fhead0["label2"] + " (%s)" % fhead0["unit2"]
    fhead0["label3"] = fhead0["label3"] + " (%s)" % fhead0["unit3"]
    fhead0.update(kargs)
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams["font.size"] = fhead0["fontsize"]

    n1tic, n2tic, n3tic = fhead0["n1tic"], fhead0["n2tic"], fhead0["n3tic"]
    # if label1.endswith("()"): label1 = label1[:-2]
    # if label2.endswith("()"): label2 = label2[:-2]
    # if label3.endswith("()"): label3 = label3[:-2]

    
    if color == 'gray':
        linecol = 'white'
        textcol = 'black'
    elif color == 'jet':
        linecol = 'black'
        textcol = 'black'
    else:
        linecol = kargs.get("linecol", "blue")
        textcol = kargs.get("textcol", "blue")


    default_linefat = 1
    default_linefat1 = 2.5
    default_linefat2 = 4.5
    default_linecol = linecol

    
    amax1, amax2, amax3 = sfaxes[0][-1], sfaxes[1][-1], sfaxes[2][-1]
    o1, o2, o3 = sfaxes[0][0], sfaxes[1][0], sfaxes[2][0]
    len1, len2, len3 = amax1 - o1, amax2 - o2, amax3 - o3
    if len1==0:
        len1 = 1
        o1 = 0
        amax1 = 1
    if len2==0:
        len2 = 1
        o2 = 0
        amax2 = 1
    if len3==0:
        len3 = 1
        o3 = 0
        amax3 = 1
    extents = [
        [o2, amax2, amax1, o1],
        [0, len2, len3, 0],
        [0, len3, 0, len1],
    ]

    axis_xlims = [
        [o2, amax2],
        [o2, amax2 + len2 / point2 * (1 - point2)],
        [o3, amax3]
    ]
    axis_ylims = [
        [amax1, o1],
        [o3, amax3],
        [o1, amax1 + len1 / point1 * (1 - point1)]
    ]


    axshear = (1 - point2) / len3 * len2 / point2
    axshift = -(1 - point2) / len3 * len2 / point2 * o3
    dxshear = (1 - point2) / len3
    dxshift = 0
    ayshear = (1 - point1) / len3 * len1 / point1
    ayshift = -(1 - point1) / len3 * len1 / point1 * o3
    dyshear = (1 - point1) / len3
    dyshift = 0

    l11 = (amax1 - sfaxes[0][frame1]) / len1 * point1
    l12 = (sfaxes[1][frame2] - o2) / len2 * point2

    loff1 = (sfaxes[2][frame3] - o3) / len3 * (1 - point1)
    loff2 = (sfaxes[2][frame3] - o3) / len3 * (1 - point2)

    # if reverse1:
    #     axis_ylims[0] = axis_ylims[0][::-1]
    #     extents[2][2], extents[2][3] = (len1, 0)
    #     l11 = (sfaxes[0][frame1] - o1) / len1 * point1

    # if reverse2:
    #     axis_xlims[0] = axis_xlims[0][::-1]
    #     extents[1][0], extents[1][1] = (len2, 0)
    #     l12 = (amax2 - sfaxes[1][frame2]) / len2 * point2

    # if reverse3:
    #     axis_ylims[1] = axis_ylims[1][::-1]
    #     axis_xlims[1] = axis_xlims[1][::-1]
    #     extents[1][2], extents[1][3] = (0, len3)
    #     extents[2][0],extents[2][1] = (len3, 0)
    #     loff1 = (amax3 - sfaxes[2][frame3]) / len3 * (1 - point1)
    #     loff2 = (amax3 - sfaxes[2][frame3]) / len3 * (1 - point2)




    dtrans1 = transforms.Affine2D(
        np.array([[point2 / len2, dxshear, 0],
         [0, 1 / len3, 0],
         [0, 0, 1]])
    )
    dtrans2 = transforms.Affine2D(
        np.array([[1 / len3, 0, 0],
         [dyshear, point1 / len1, 0],
         [0, 0, 1]])
    )
    atrans1 = transforms.Affine2D(
        np.array([[1, axshear, axshift],
         [0, 1, 0],
         [0, 0, 1]])
    )
    atrans2 = transforms.Affine2D(
        np.array([[1, 0, 0],
         [ayshear, 1, ayshift],
         [0, 0, 1]])
    )
    grid_helper = GridHelperCurveLinear(atrans1,
                                        extremes=[o2, amax2, o3, amax3],
                                        grid_locator1=MaxNLocator(nbins=n2tic),
                                        grid_locator2=MaxNLocator(nbins=n3tic))
    grid_helper1 = GridHelperCurveLinear(atrans2,
                                         extremes=[o3, amax3, o1, amax1],
                                         grid_locator1=MaxNLocator(nbins=n1tic),
                                         grid_locator2=MaxNLocator(nbins=n1tic)
                                         )

    topmargin = 0.05

    if figure is not None:
        fig = figure
        axbase0 = fig.gca()
    else:
        axbase0 = plt.gca()
        fig = plt.gcf()
    axbase0.axis('off')
    axbase0.set_label("Base axes")

    if fhead0["title"]!="":
        title_pos = [0, 1 - topmargin, 1, topmargin]
        ax_title = axbase0.inset_axes(title_pos)
        ax_title.set_frame_on(False)
        ax_title.get_xaxis().set_visible(False)
        ax_title.get_yaxis().set_visible(False)
        ax_title.set_xlim([0, 1])
        ax_title.set_ylim([0, 1])
        ax_title.text(0.5, 1, fhead0["title"], va='bottom', ha='center', fontsize=int(fhead0["fontsize"])+5)
        axbase0.add_child_axes(ax_title)
    else:
        topmargin = 0.
    hei_ax = 1 - topmargin
    

    if fhead0["colorbar"]:
        axbase = axbase0.inset_axes(bounds=[0, 0, 0.85, 1 - topmargin], facecolor="none")
        axbar = axbase0.inset_axes(bounds=[0.90, 0, 0.05, 1 - topmargin], facecolor="none")
        axbase.set_label("Base axes with bar")
    else:
        axbase = axbase0
    axbase.axis('off')

    ax1 = axbase.inset_axes(bounds=[0, 0, point2, point1*hei_ax], facecolor="none")
    ax1.set_label("Axes 1 [n1, n2]")

    im0 = ax1.imshow(slice1, aspect="auto", extent=extents[0],
                     cmap=usecmap, vmin=vmin, vmax=vmax)

    ax1.yaxis.set_major_locator(MaxNLocator1(nbins=n1tic))
    ax1.set_xlim(axis_xlims[0])
    ax1.set_ylim(axis_ylims[0])

    ax2 = axbase.inset_axes(bounds=[0, point1*hei_ax, 1, (1 - point1)*hei_ax],
                            axes_class=Axes, grid_helper=grid_helper,
                            facecolor="none")
    ax2.set_label("Axes 3 [n1, n3]")

    im1 = ax2.imshow(slice3, aspect="auto", extent=extents[1],
                     cmap=usecmap, vmin=vmin, vmax=vmax)
 
    im1.set_transform(dtrans1 + ax2.transAxes)

    #
    ax2.set_ylim(axis_ylims[1])
    ax2.set_xlim(axis_xlims[1])
    # ax2.axis["t1"] = ax2.new_floating_axis(0,-5)

    ax3 = axbase.inset_axes(bounds=[point2, 0, 1 - point2, hei_ax],
                            axes_class=Axes, grid_helper=grid_helper1,
                            facecolor="none")
    ax3.set_label("Axes 2 [n2, n3]")

    im2 = ax3.imshow(slice2, aspect="auto", extent=extents[2],
                     cmap=usecmap, vmin=vmin, vmax=vmax)
    im2.set_transform(dtrans2 + ax3.transAxes)
    #
    ax3.set_xlim(axis_xlims[2])
    ax3.set_ylim(axis_ylims[2])

    #  Details
    # Ticks
    # ### Calculate deformation for ticklabels rotation
    # box_width = axbase.bbox.width
    # box_height = axbase.bbox.height
    # rot_angel = np.arctan(point1/point2*box_height/box_width)/np.pi*180
    ax2.axis[:].major_ticklabels.set(visible=False)
    ax2.axis[:].major_ticks.set(visible=False)

    # if reverse3:
    #     ax2.axis["right"].major_ticks.set(tick_out=True, visible=True)
    #     ax2.axis["right"].major_ticklabels.set(visible=True,rotation=90, ha="center")
    #     ax2.axis["right"].label.set(visible=True,rotation=180, ha="right",va="bottom")
    #     ax2.axis["left"].set(visible=False)
    #     ax2.axis["right"].major_ticklabels.set(visible=True)

    
    ax2.axis["left"].major_ticks.set(tick_out=True, visible=True)
    ax2.axis["right"].set(visible=False)
    ax2.axis["left"].major_ticklabels.set(visible=True)



    ax3.axis[:].major_ticklabels.set(visible=False)
    ax3.axis[:].major_ticks.set(visible=False)
    ax3.axis["right"].major_ticks.set(visible=False)
    # bottom axis unused
    # ax3.axis["bottom"].major_ticklabels.set(rotation=-rot_angel,visible=True)
    # ax3.axis["bottom"].major_ticks.set(tick_out=True,visible=True)

    # Labels
    ax1.set_xlabel(fhead0["label2"])
    ax1.set_ylabel(fhead0["label1"])

    # if reverse3:
    #     ax2.axis["right"].label.set(text=label3)
    # else:
    ax2.axis["left"].label.set(text=fhead0["label3"])
    # ax3.axis["bottom"].label.set(text=label3)

    # Title
    # if fhead0["title"] != "":
    #     axbase0.set_title(fhead0["title"])

    # Colorbar
    if fhead0["colorbar"]:
        # cbar = fig.colorbar(im0, cax=axbar, label=bartitle, orientation=barpos)
        cbar = mbar.Colorbar(axbar, im0, cmap=usecmap, label=fhead0["bartitle"])

    # Indicating lines

    l21 = point1 + loff1
    l22 = point2 + loff2
    l221 = l12 +  (1-point2)

    l31 = l11 + 1 - point1
    l32 = point1 + loff1

    ax1.hlines(l11*hei_ax,0,point2, color=default_linecol, linewidth=default_linefat1,transform=axbase.transAxes)
    ax1.vlines(l12,0,point1*hei_ax, color=default_linecol, linewidth=default_linefat1,transform=axbase.transAxes)
    ax2.plot([loff2,l22],[l21*hei_ax,l21*hei_ax], color=default_linecol, linewidth=default_linefat1,transform=axbase.transAxes)
    ax2.plot([l12,l221],[point1*hei_ax,1*hei_ax], color=default_linecol, linewidth=default_linefat1,transform=axbase.transAxes)
    ax3.vlines(1, (1-point1)*hei_ax, 1*hei_ax, color=default_linecol, linewidth=default_linefat, transform=axbase.transAxes)
    ax3.vlines(l22, loff1*hei_ax, l32*hei_ax, color=default_linecol, linewidth=default_linefat1, transform=axbase.transAxes)
    ax3.plot([point2,1],[l11*hei_ax,l31*hei_ax], color=default_linecol, linewidth=default_linefat1,transform=axbase.transAxes)

    # Indicating labels
    lab1 = str(sfaxes[1][frame2])
    lab11 = "%.2f" % sfaxes[1][frame2]
    lab1 = lab11 if len(lab11) < len(lab1) else lab1
    lab2 = str(sfaxes[0][frame1])
    lab21 = "%.2f" % sfaxes[0][frame1]
    lab2 = lab21 if len(lab21) < len(lab2) else lab2
    lab3 = str(sfaxes[2][frame3])
    lab31 = "%.2f" % sfaxes[2][frame3]
    lab3 = lab31 if len(lab31) < len(lab3) else lab3
    axbase.text(l221, 1*hei_ax, "%s" % lab1, va='bottom',
                ha='center', color=default_linecol)
    axbase.text(1, l31*hei_ax, "%s" % lab2, va='center',
                ha='left', color=default_linecol, rotation=-90)
    axbase.text(l22, loff1*hei_ax,
                "%s" % lab3, va='top', ha='left', color=default_linecol)

    # return ([ax1, ax2, ax3,  # Axes: child ax 0-2, base ax and cbar ax
    #          axbase0,axbase, 
    #          axbar if colorbar else None], 
    #         [im0, im1, im2], # ims
    #         cbar if colorbar else None #colorbar
    #         )
    fig.add_axes(ax1)
    fig.add_axes(ax3)
    fig.add_axes(ax2)


def grey3(*args,**kwargs):
    if "flat" not in kwargs.keys():
        return grey3flat(*args,**kwargs)
    elif kwargs["flat"]== False: return grey3noflat(*args,**kwargs)
    else: return grey3flat(*args,**kwargs)


def wiggle(rsfdata, figure=None, **kargs):
    '''
    Plot wiggle traces.

    Parameters
    ----------
    rsfdata : rsfpy.Rsf
        Input data.
    figure : matplotlib.figure.Figure, optional
        Figure to plot on. If None, a new figure will be created.
    **kargs : dict
        Parameters for plotting. Details are as follows.
    pclip : float, optional
        Clip data at pclip percentile. Default is 99.
    vmin : float, optional
        Minimum value for color scale. Default is None.
    vmax : float, optional
        Maximum value for color scale. Default is None.
    zplot : float, optional
        Vertical spacing between traces. Default is 1.
    fill : bool, optional
        Fill positive and negative values with different colors. Default is True.
    interpolate : bool, optional
        Interpolate between points. Default is False.
    maxtrc : int, optional
        Maximum number of traces to plot. Default is 100.
    min1 : float, optional
        Minimum value for axis 1. Default is None.
    max1 : float, optional
        Maximum value for axis 1. Default is None.
    min2 : float, optional
        Minimum value for axis 2. Default is None.
    max2 : float, optional
        Maximum value for axis 2. Default is None.
    colorbar : bool, optional
        Show colorbar. Default is False.
    poscolor : str, optional
        Color for positive values. Default is "red".
    negcolor : str, optional
        Color for negative values. Default is "blue".
    color : str, optional
        Color for traces. Default is "black".
    linewidth : float, optional
        Line width for traces. Default is 1.
    label1 : str, optional
        Label for axis 1. Default is "Time".
    unit1 : str, optional
        Unit for axis 1. Default is "".
    label2 : str, optional
        Label for axis 2. Default is "Distance".
    unit2 : str, optional
        Unit for axis 2. Default is "".
    share_ax : matplotlib.axes.Axes, optional
        Share axes with another axes. Default is None.
    '''
    def float1(x):
        try:return float(x)
        except:return None
    def int1(x):
        try:return int(x)
        except:return 0
    
    
    if hasattr(rsfdata, "header"):
        input_args = rsfdata.header
    else:
        input_args = {}
    input_args.update(kargs)
    pclip = float1(input_args.get("pclip", 99))
    vmin = float1(input_args.get("vmin", None))
    vmax = float1(input_args.get("vmax", None))
    zplot = float1(input_args.get("zplot", 1))
    fill = input_args.get("fill", True)
    interpolate = input_args.get("interpolate", False)
    # interval2 = float1(kargs.get("interval2", 1))
    maxtrc = int1(input_args.get("maxtrc", 100))
    min1 = float1(input_args.get("min1", None))
    max1 = float1(input_args.get("max1", None))
    min2 = float1(input_args.get("min2", None))
    max2 = float1(input_args.get("max2", None))
    colorbar = input_args.get("colorbar", False)
    poscolor = input_args.get("poscolor", "red")
    negcolor = input_args.get("negcolor", "blue")
    color = input_args.get("color", "black")
    linewidth = float1(input_args.get("linewidth", 1))
    label1 = input_args.get("label1", "Time")
    unit1 = input_args.get("unit1", "")
    label2 = input_args.get("label2", "Distance")
    unit2 = input_args.get("unit2", "")
    share_ax = input_args.get("share_ax", None)
    use_ax = input_args.get("use_ax", None)
    plotstyle = input_args.get("plotstyle", {})
    colorthres = input_args.get("colorthres", 0.1) # Threshold value for filling color

    if pclip >=0 and pclip < 100:
        vmax = np.percentile(np.abs(rsfdata.view(np.ndarray)), pclip)
        vmin = -vmax

    elif vmin is None and vmax is None:
        vmax = np.max(rsfdata)
        vmin = np.min(rsfdata)
    elif vmin is None:
        vmin = -vmax
    elif vmax is None:
        vmax = -vmin
    if vmax < vmin:
        vmin, vmax = vmax, vmin
    clipbias = (vmax + vmin)/2

    if linewidth < 0:
        linewidth = 1

    datashape = rsfdata.shape
    if len(datashape) < 2:
        datashape = (datashape[0], 1)
    if len(datashape) > 2:
        datashape = (datashape[0], datashape[1])
        data = rsfdata[:,:,0]
    else:
        data = rsfdata

    # Downsample if needed
    # if interval2 >= 2:
    #     data = rsfdata[:,::interval2]
    # else:
    #     data = rsfdata
    interval2 = 1
    while datashape[1] > maxtrc:
        data = data[:,::2]
        datashape = data.shape
        interval2 *= 2

    # clip data
    data = np.clip(data.view(np.ndarray), vmin, vmax)
    # data = data.view(np.ndarray)

    # Get axis parameters
    if hasattr(rsfdata, "axis"):
        axis2 = rsfdata.axis()
        axist = axis2[0]
        axisx = axis2[1]
        if interval2 >= 2:
            axisx = axisx[::interval2]
    else:
        axist = np.arange(data.shape[0])
        axisx = np.arange(data.shape[1])
    
    # Normalize data
    xspace = abs(axisx[-1]  - axisx[0]) / len(axisx) * zplot
    # data = data / max(abs(vmin), abs(vmax)) * xspace * zplot + axisx[np.newaxis, :]
    data = (data- clipbias) / (vmax - clipbias) * xspace + axisx[np.newaxis, :]

    # Plot data
    if figure is None:
        fig = plt.figure()
        axbase = fig.add_subplot(111)
    else:
        fig = figure
        axbase = fig.gca()
    axbase.axis('off')
    axbase.set_label("Base axes")
    if use_ax is not None:
        ax = use_ax
        if ax == axbase:
            axbase.axis('on')
    else:
        if colorbar:
            if share_ax is not None:
                ax = axbase.inset_axes([0,0,1-0.15,1],facecolor='none', sharex=share_ax, sharey=share_ax)
            else:ax = axbase.inset_axes([0,0,1-0.15,1],facecolor='none')
        else:
            if share_ax is not None:
                ax = axbase.inset_axes([0,0,1,1],facecolor='none', sharex=share_ax, sharey=share_ax)
            else:ax = axbase.inset_axes([0,0,1,1],facecolor='none')

    for i in range(axisx.size):
        if fill:
            ax.fill_betweenx(axist, axisx[i], data[:, i], where=data[:, i] > axisx[i] + colorthres*xspace*zplot, color=poscolor, interpolate=interpolate)
            ax.fill_betweenx(axist, axisx[i], data[:, i], where=data[:, i] < axisx[i] - colorthres*xspace*zplot, color=negcolor, interpolate=interpolate)
        icurve = ax.plot(data[:, i], axist, color=color, linewidth=linewidth,**plotstyle)
        icurve[0].set_label(f"Trace{i}")

    
    ax.set_xlim(axisx[0]-xspace * zplot, axisx[-1]+ xspace * zplot)
    ax.set_ylim(axist[-1], axist[0])
    if min1 is not None and max1 is not None:
        ax.set_ylim(max1, min1)
    if min2 is not None and max2 is not None:
        ax.set_xlim(min2, max2)

    if unit1 != "":
        label1 = label1 + " (%s)" % unit1
    if unit2 != "":
        label2 = label2 + " (%s)" % unit2
    ax.set_xlabel(label2)
    ax.set_ylabel(label1)

    if use_ax is None: fig.add_axes(ax)
    if share_ax is not None:
        fig.sca(share_ax)
    else:fig.sca(ax)
    ax.set_label("Wiggle plot")