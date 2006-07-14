spline[s_, x_, y_] := 
	1/10 (114 - 73 s) + 
	8/15 (-11 + 8 s) (Cos[x] + Cos[y]) - 
	7/15 (-1 + s) (Cos[2 x] + Cos[2 y]) + 
      	28/45 (-3 + 2 s) (Cos[x - y] + Cos[x + y]) + 
      	1/180 (-6 + 5 s) (Cos[2 (x - y)] + Cos[2 (x + y)]) - 
      	4/45 (-9 + 8 s) (Cos[x-2 y]+Cos[2 x-y]+Cos[2 x+y]+Cos[x+2 y]);
cp[s_]:=ContourPlot[spline[s,x,y],{x,-Pi,Pi},{y,-Pi,Pi},
    	DisplayFunction->Identity,PlotRange->{0,16}];
Show[GraphicsArray[{
	Show[cp[#],PlotLabel->SequenceForm["tension=",#],TextStyle->{FontSize->18}]& /@ {0,0.3},
	Show[cp[#],PlotLabel->SequenceForm["tension=",#],TextStyle->{FontSize->18}]& /@ {0.7,1}}]];
Display["junk_ma.eps",%,"EPS",ImageSize->{576,576}];

