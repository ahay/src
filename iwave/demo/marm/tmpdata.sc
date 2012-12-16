import os
Import('MADAGASCAR', 'TMPDATAPATH', 'fetches')

try: 
    os.mkdir(TMPDATAPATH)
except OSError:
    pass
try:
    os.mkdir('./tmpdata')
except OSError:
    pass




