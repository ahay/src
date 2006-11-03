function points(in,out)
%POINTS plot points

dims = rsf_dim(in);
pnts = zeros(dims');
rsf_read(pnts,in);
x = pnts(1,:);
y = pnts(2,:);
z = pnts(3,:);
plot3(x,y,z,'.');
view(15,30);
print('-depsc',out);

