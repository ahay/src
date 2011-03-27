bash$ scons 
scons: Building targets ...
/usr/bin/sfspike n1=1000 n2=100 | /usr/bin/sfbandpass fhi=10 > spike.rsf
< spike.rsf /usr/bin/sfclip clip=0.5 > cliped.rsf
scons: Done building targets.
bash$ sed -i'' -e's/0.5/0.25/' SConstruct
bash$ scons -Q
< spike.rsf /usr/bin/sfclip clip=0.25 > cliped.rsf

