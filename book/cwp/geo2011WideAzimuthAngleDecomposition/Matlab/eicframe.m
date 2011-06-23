%function eicframe()
eicdims
hold on

red  =[1 0 0];
green=[0.25 0.65 0.1];
blue =[0 0 1];
cyan =[0 1 1];
yellow=[1 1 0];
black=[0 0 0];
brown=[210 105 30]/255;
white=[1 1 1];
magenta=[255 0 255]/255;
%brown=[139 69 19]/255.;

cmap=[ red ; green ; blue ; yellow ; magenta ; brown ; white ; black ];
colormap(cmap)
caxis([0 7])

% corner dots
bco=2.0;
xco=[-bco -bco -bco -bco +bco +bco +bco +bco];
yco=[-bco -bco +bco +bco -bco -bco +bco +bco];
zco=[-bco +bco -bco +bco -bco +bco -bco +bco];
plot3(xco,yco,zco,'w.');

asc=1.05;
orx=[-asc*bco  0  0];
ory=[0  -asc*bco  0];
orz=[0  0  -asc*bco];
nnx=[+asc*bco  0  0];
nny=[0  +asc*bco  0];
nnz=[0  0  +asc*bco];

axis([-bco +bco -bco +bco -bco +bco]);
daspect([1 1 1])
grid on
axis off

% coordinate system axes
cor=1.99;
h=text(-cor+0.6,-cor,-cor,'x'); set(h,'FontSize',10,'Color','k');
h=text(-cor,-cor+0.6,-cor,'y'); set(h,'FontSize',10,'Color','k');
h=text(-cor,-cor,-cor+0.6,'z'); set(h,'FontSize',10,'Color','k');

quiver3([-cor],[-cor],[-cor],[+1],[0],[0],...
        0.5,'Color','k','LineWidth',2);
quiver3([-cor],[-cor],[-cor],[0],[+1],[0],...
        0.5,'Color','k','LineWidth',2);
quiver3([-cor],[-cor],[-cor],[0],[0],[+1],...
        0.5,'Color','k','LineWidth',2);


