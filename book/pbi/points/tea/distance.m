function distance(in,out,level)
%DISTANCE plot isosurface

dims = rsf_dim(in);
dist = zeros(dims');
rsf_read(dist,in);
isosurface(dist,level);
view(75,30);
print('-depsc',out);


