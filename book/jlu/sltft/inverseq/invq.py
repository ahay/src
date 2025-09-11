from rsf.proj import Flow,Plot,Result

def invq(data,name,tfxmap,dw=1):
    pf = name + "_pf"         # peak frequency
    lcf = name + "_lcf"       # local frequency
    eqvq = name + "_eqvq"     # equivalent Q value
    intq = name + "_intq"     # interval Q value
    result = name + "result"  # filtered
    spec = name + "spec"      # spectrum

    Flow(pf,tfxmap,
        '''
        cabs |
        transp |
        findmax1 verb=n|
        dd type=float | math output="(input-1)*%g*2.5"|
        dd type=float | smooth rect1=100
        ''' % dw)

    Flow([lcf, name+'_cf', name+'_var2'],[tfxmap, pf],
         '''
         window max2=100|
         lcf trange=${SOURCES[1]} avef=${TARGETS[1]}
         var2=${TARGETS[2]} rect1=100
         ''')
    
    Flow(name+'_repos',None,
         'math n1=1 o1=0. d1=1 n2=247 output="0" | dd type=int')
    Flow(name+'_qts',[ name+'_cf', name+'_var2', name+'_repos'],
         '''
         lcfseq var2=${SOURCES[1]} repos=${SOURCES[2]}
         ''')
    
    Flow(eqvq, name+'_qts','smooth rect1=1 rect2=10')
    
    Flow(intq,[ name+'_cf', name+'_var2'] ,
         '''
         lcfsiq var2=${SOURCES[1]} |
         clip clip=500 |
         smooth rect1=1 rect2=20
         ''')
    
    Flow(name+'_rasf',data,'fft1 | spray axis=2 d=0.004 n=1126 o=1.| transp')

    Flow(name+'_rcmtf1',[name+'_rasf',eqvq],
         '''
         invqfilt eqt=${SOURCES[1]} gim=30
         ''')

    Flow(result,name+'_rcmtf1','transp | fft1 inv=y | window n1=1 f1=0')

    ## ---------frequency spectrum
    Flow(spec,result,'window min1=2.5 max1=4.5 |spectra all=y ')

    return  pf, lcf, eqvq, intq, result, spec
    
