def Bay(name,f2=500,n2=1600,f1=400,n1=700):
    '''Extracts the Bay Area dataset (digital topography) and windows it'''
    
    datdir  = 'ftp://begpc132.beg.utexas.edu/data/bay'
    Fetch('bay.H',datdir)

    Flow(name,'bay.H',
         '''
         dd data_format=native_float |
         window f2=%d n2=%d f1=%d n1=%d |
         reverse which=1 |
         costaper nw=50
         ''' % (f2,n2,f1,n1))
