def readslice(inputfilename,data0):
#def readslice(inputfilename,ndim,data):
 #   shape = (ndim,ndim,ndim)
    fd = open(fname, 'rb')
   # data = np.fromfile(file=fd, dtype=np.double).reshape(shape)
    data0 = np.fromfile(file=fd, dtype=np.single)
    fd.close()
    return data0


