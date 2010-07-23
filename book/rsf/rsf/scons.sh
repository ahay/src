bash$ scons 
scons: Building targets ...
sfspike n1=1000 n2=100 > spike.rsf
< spike.rsf sfclip clip=0.5 > cliped.rsf
scons: Done building targets.
bash$ sed s/0.5/0.25/ < SConstruct > SConstruct2
bash$ mv SConstruct2 SConstruct
bash$ scons
scons: Building targets ...
< spike.rsf sfclip clip=0.25 > cliped.rsf
scons: Done building targets.
