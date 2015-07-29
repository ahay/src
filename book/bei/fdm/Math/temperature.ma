xM=41;xC=20;
yM=41;yC=20;
u=Table[20 Exp[-0.05(x-xC)^2-0.05(y-yC)^2],{x,xM},{y,yM}];
p1=ListPlot3D[u,PlotRange -> {0,20},ViewPoint->{2,-2,1},Ticks->None, 
    Mesh->False, Boxed->False, Axes->False];

v=Table[0,{x,xM},{y,yM}];
v[[xC]]=u[[xC]];
p2=ListPlot3D[v, PlotRange -> {0,20},ViewPoint->{2,-2,1},Ticks->None, 
    Mesh->False,Boxed->False, Axes->False];
		
w=Table[0,{x,xM},{y,yM}];
w[[xC,yC]]=u[[xC,yC]];
p3=ListPlot3D[w,PlotRange -> {0,20},ViewPoint->{2,-2,1},Ticks->None, 
    Mesh->False,Boxed->False, Axes->False];

p=Show[
	Graphics[
		{
			Rectangle[{0,0},{1,1},p3],
			Rectangle[{1,0},{2,1},p2],
			Rectangle[{2,0},{3,1},p1]
	}
], ImageSize->500
];

Display["junk_ma.eps",p,"EPS"];

