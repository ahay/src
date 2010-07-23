bash$ sfspike n1=1000 n2=100 > spike.rsf
bash$ < spike.rsf sfclip clip=0.5 > cliped.rsf
bash$ sfspike n1=1000 n2=100 | sfclip clip=0.5 > cliped2.rsf
