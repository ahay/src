n=1;
d[1]=1;

A = Sum[ (d[i] - m)^2,{i,1,n}];
B =  Sum[ Abs[d[i] - m],{i,1,n}];
CCC= Sum[ (Abs[d[i] - m])^(0.1) ,{i,1,n}];

a = Plot [A,{m,0,2}];
b = Plot [B,{m,0,2},PlotRange->{0.,1.2}];
c = Plot [CCC,{m,0,2},PlotRange->{0.,1.2}];

n=5;
d[1]=1 ; d[2]=2 ; d[3]=3 ; d[4]=6 ; d[5]=6;
d[1]=1 ; d[2]=2 ; d[3]=3 ; d[4]=5 ; d[5]=5;
d[1]=1 ; d[2]=1 ; d[3]=1 ; d[4]=2 ; d[5]=5;
d[1]=1 ; d[2]=1 ; d[3]=2 ; d[4]=3 ; d[5]=5;

A = Sum[ (d[i] - m)^2,{i,1,n}];
B =  Sum[ Abs[d[i] - m],{i,1,n}];
CCC= Sum[ (Abs[d[i] - m])^(0.1) ,{i,1,n}];

aa=Plot [A,  {m,0,6},PlotRange->{0.0,60}];
bb=Plot [B,  {m,0,6},PlotRange->{0.0,15}];
cc=Plot [CCC,{m,0,6},PlotRange->{0.0,7},PlotPoints->100];

aaa=Show[GraphicsArray[{{a,b,c},{aa,bb,cc}}]];

Display["junk_ma.eps", aaa, "EPS", ImageSize->1000];



