function eicvect(aa,c,l,w)

switch nargin
 case 3
  w=3;
 case 2
  w=3;
 case 1 
  w=3;
end

%aa=eicnorm(aa);
ax=aa(1);
ay=aa(2);
az=aa(3);

asc=1.0;
%quiver3([-asc*ax],[-asc*ay],[-asc*az],...
%        [+asc*ax],[+asc*ay],[+asc*az],...
%        asc*2.0,'Color',c,'LineWidth',2);

quiver3([0],[0],[0],...
        [+asc*ax],[+asc*ay],[+asc*az],...
        asc*1.0,'Color',c,'LineWidth',w);

lsc=1.10;
h=text(lsc*ax+0.05,lsc*ay+0.05,lsc*az+0.05,l); 
set(h,'FontSize',15,'Color',c);