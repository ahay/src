from rsf.proj import *
Import('MADAGASCAR TMPDATAPATH ROOT')

os.system('cd marmousi; DATAPATH=' + TMPDATAPATH + '/; export DATAPATH; $RSFROOT/bin/sfdd form=native < ' + TMPDATAPATH + '/velocity.HH > velocity.rsf; $RSFROOT/bin/sfdd form=native < ' + TMPDATAPATH + '/density.HH > density.rsf')