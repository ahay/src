Import('MADAGASCAR', 'TMPDATAPATH')

#MADAGASCAR = False
#TMPDATAPATH = '/var/tmp/iwave'

if MADAGASCAR:
    from rsf.proj import *
# plot model fields
    Result('vp2d_5m', './vp2d_5m/vp2d_5m.rsf','put label1=Depth unit1=m label2 = Distance unit2=m label=Velocity unit=m/ms | grey color=3 mean=y title="Dome Velocity Model dx = dz = 5 m" scalebar=y barreverse=y')
    Result('dn2d_5m', './dn2d_5m/dn2d_5m.rsf','put label1=Depth unit1=m label2 = Distance unit2=m label=Density unit=g/cm^3 | grey color=3 mean=y title="Dome Velocity Model dx = dz = 5 m" scalebar=y barreverse=y')
# 2-4 results
    traces = []
    gathers =['demo20m', 'demo10m', 'demo5m', 'demo2.5m']
    for r in range(4):
        data = './' + gathers[r] + '.rsf'
        Flow(data,'./' + gathers[r] + '/data.su','suread < $SOURCE read=data endian=0',stdin=0)
        trace = 'trace' + str(r)
        print trace
        Flow(trace,data,'window n2=1 f2=100')
        traces.append(trace)
    Result('data2.5m','./demo2.5m.rsf','grey title="Dome Data, 2-4 scheme, dx = dz = 2.5 m"')
    Result('trace',traces,
       '''
       cat axis=2 ${SOURCES[1:4]} | 
       graph plotcol=7,6,2,4 wanttitle=n label2=Pressure unit2=MPa
       ''')
    Result('wtrace',traces,
       '''
       cat axis=2 ${SOURCES[1:4]} | window min1=1.8 max1=2.5 |
       graph plotcol=7,6,2,4 wanttitle=n label2=Pressure unit2=MPa
       ''')
# 2-8 results
    traces8k = []
    gathers8k =['demo20m8k', 'demo10m8k', 'demo5m8k', 'demo2.5m8k']
    for r in range(4):
        data = './' + gathers8k[r] + '.rsf'
        Flow(data,'./' + gathers8k[r] + '/data.su','suread < $SOURCE read=data endian=0',stdin=0)
        trace = 'trace8k' + str(r)
        print trace
        Flow(trace,data,'window n2=1 f2=100')
        traces8k.append(trace)
    Result('data2.5m8k','./demo2.5m8k.rsf','grey title="Dome Data, 2-8 scheme, dx = dz = 2.5 m"')
    Result('trace8k',traces8k,
       '''
       cat axis=2 ${SOURCES[1:4]} | 
       graph plotcol=7,6,2,4 wanttitle=n label2=Pressure unit2=MPa
       ''')
    Result('wtrace8k',traces8k,
       '''
       cat axis=2 ${SOURCES[1:4]} | window min1=1.8 max1=2.5 |
       graph plotcol=7,6,2,4 wanttitle=n label2=Pressure unit2=MPa
       ''')
    End()
else:
    import os
    os.system(
        '''suwind key=tracl min=100 max=100 < ./demo20m/data.su >  ./Fig/trace100.su;
           suwind key=tracl min=100 max=100 < ./demo10m/data.su >> ./Fig/trace100.su;
           suwind key=tracl min=100 max=100 < ./demo5m/data.su  >> ./Fig/trace100.su;
           suwind key=tracl min=100 max=100 < ./demo2.5m/data.su  >> ./Fig/trace100.su;
           suwind key=tracl min=100 max=100 < ./demo20m8k/data.su >  ./Fig/trace1008k.su;
           suwind key=tracl min=100 max=100 < ./demo10m8k/data.su >> ./Fig/trace1008k.su;
           suwind key=tracl min=100 max=100 < ./demo5m8k/data.su  >> ./Fig/trace1008k.su;
           suwind key=tracl min=100 max=100 < ./demo2.5m8k/data.su  >> ./Fig/trace1008k.su;
        ''')
    os.system('psimage n1=361 f2=0 d2=0.005 f1=0 d1=0.005 label1="depth (km)" label2="offset (km)" legend=1 lstyle=horibottom units="velocity (km/s)" d1num=0.5 width=6 height=4 wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0 < ' + TMPDATAPATH + '/vp2d_5m.rsf@ > ./Fig/fig1.ps;')
    os.system('psimage n1=361 f2=0 d2=0.005 f1=0 d1=0.005 label1="depth (km)" label2="offset (km)" legend=1 lstyle=horibottom units="density (kg/m^3)" d1num=0.5 width=6 height=4 wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0 < ' + TMPDATAPATH + '/dn2d_5m.rsf@ > ./Fig/fig2.ps;')
    os.system(
        '''supsimage f2=0.1 d2=0.02 width=6 height=6 label1="time (s)" label2="gx (km)" < ./demo2.5m/data.su > ./Fig/fig3.ps;
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (ms)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 < ./Fig/trace100.su > ./Fig/fig4.ps;
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (ms)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 x1beg=1.8 x1end=2.5 x2beg=-1e+2 x2end=1e+2 < ./Fig/trace100.su > ./Fig/fig5.ps;
           supsimage f2=0.1 d2=0.02 width=6 height=6 label1="time (s)" label2="gx (km)" < ./demo5m8k/data.su > ./Fig/fig6.ps;
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (ms)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 < ./Fig/trace1008k.su > ./Fig/fig7.ps;   
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (ms)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 x1beg=1.8 x1end=2.5 x2beg=-1e+2 x2end=1e+2 < ./Fig/trace1008k.su > ./Fig/fig8.ps;
    ''')

