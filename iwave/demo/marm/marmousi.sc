from rsf.proj import *
Import('MADAGASCAR TMPDATAPATH ROOT')

os.system('cd marmousi; DATAPATH=' + TMPDATAPATH + '/; export DATAPATH; $RSFROOT/bin/sfdd form=native < ' + TMPDATAPATH + '/velocity.HH > velocity.rsf; $RSFROOT/bin/sfdd form=native < ' + TMPDATAPATH + '/density.HH > density.rsf')

os.system('cd marmousi; DATAPATH=' + TMPDATAPATH + '/; export DATAPATH; $RSFROOT/bin/sfwindow j1=6 j2=6 < ./velocity.rsf > ./vel24.rsf; $RSFROOT/bin/sfwindow j1=6 j2=6 < ./density.rsf > ./den24.rsf')

os.system('suwind key=sx j=100 < ' + TMPDATAPATH + '/line_fix.su > ' + TMPDATAPATH + '/line100.su')
