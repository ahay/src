Import('MADAGASCAR TMPDATAPATH ROOT')

if MADAGASCAR:
    from rsf.proj import *
# plot model fields
    os.system('< ./vp2d_5m/vp2d_5m.rsf ' + ROOT + 'sfput label1=Depth unit1=m label2 = Distance unit2=m label=Velocity unit=m/ms | ' + ROOT + 'sfgrey color=3 mean=y title="Dome Velocity Model dx = dz = 5 m" scalebar=y barreverse=y > Figs/vp2d_5m.vpl')
    os.system('< ./dn2d_5m/dn2d_5m.rsf ' + ROOT + 'sfput label1=Depth unit1=m label2 = Distance unit2=m label=Density unit=g/cm^3 | ' + ROOT + 'sfgrey color=3 mean=y title="Dome Velocity Model dx = dz = 5 m" scalebar=y barreverse=y > Figs/dn2d_5m.vpl')
# 2-4 results
    traces = ''
    gathers =['demo20m', 'demo10m', 'demo5m', 'demo2.5m']
    for r in range(4):
        data = './' + gathers[r] + '.rsf'
	os.system('< ' + gathers[r] + '/data.su ' + ROOT + 'sfsuread read=data endian=0 > ' + data)
        trace = 'trace' + str(r)
        os.system('< ' + data + ' ' + ROOT + 'sfwindow n2=1 f2=100 > ' + trace)
        traces += trace + ' '
    os.system('< ./demo2.5m.rsf ' + ROOT + 'sfgrey title="Dome Data, 2-4 scheme, dx = dz = 2.5 m" > Figs/data2.5m.vpl')
    os.system(ROOT + 'sfcat axis=2 ' + traces + ' | ' + ROOT + 'sfgraph plotcol=7,6,2,4 wanttitle=n label2=Pressure unit2=MPa > Figs/trace.vpl')
    os.system(ROOT + 'sfcat axis=2 ' + traces + ' | ' + ROOT + 'sfwindow min1=1.8 max1=2.5 | ' + ROOT + 'sfgraph plotcol=7,6,2,4 wanttitle=n label2=Pressure unit2=MPa > Figs/trace.vpl')
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

