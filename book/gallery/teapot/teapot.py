from rsf.proj import *
import rsf.gallery 

Flow('vel.asc',None,
     '''
     echo 
                           0.00        9132.86
                          .47950       10553.29
                          .63794       10921.60
                          .79222       10791.97
                          .87769       11074.19
                         1.00903       11649.54
                         1.10493       11807.96
                         1.19458       12325.03
                         1.61570       14410.47
                         3.01000       17216.64
     in=$TARGET data_format=ascii_float n1=2 n2=10 
     ''')
Flow('vel','vel.asc','dd form=native | transp')

def get_vrms1(vrms):
    Flow(vrms,'vel','linear n1=2049 o1=0 d1=0.002')

def get_vint1(vint):
    Flow(vint,'vel','linear n1=2049 o1=0 d1=0.002 | dix niter=100 rect1=15')
