Import('MADAGASCAR TMPDATAPATH ROOT')
import os

if MADAGASCAR:
    os.system('scons -f madfig.sc')
else:
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
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (s)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 < ./Fig/trace100.su > ./Fig/fig4.ps;
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (s)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 x1beg=1.8 x1end=2.5 x2beg=-1e+2 x2end=1e+2 < ./Fig/trace100.su > ./Fig/fig5.ps;
           supsimage f2=0.1 d2=0.02 width=6 height=6 label1="time (s)" label2="gx (km)" < ./demo5m8k/data.su > ./Fig/fig6.ps;
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (s)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 < ./Fig/trace1008k.su > ./Fig/fig7.ps;   
           supsgraph nplot=4 linecolor=red,green,blue,black label1="time (s)" label2="pressure (MPa)" style=normal hbox=6 wbox=6 x1beg=1.8 x1end=2.5 x2beg=-1e+2 x2end=1e+2 < ./Fig/trace1008k.su > ./Fig/fig8.ps;
    ''')

