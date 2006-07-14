spline[s_, x_, y_] := 
	1/10 (114 - 73 s) + 
	8/15 (-11 + 8 s) (Cos[x] + Cos[y]) - 
	7/15 (-1 + s) (Cos[2 x] + Cos[2 y]) + 
      	28/45 (-3 + 2 s) (Cos[x - y] + Cos[x + y]) + 
      	1/180 (-6 + 5 s) (Cos[2 (x - y)] + Cos[2 (x + y)]) - 
      	4/45 (-9 + 8 s) (Cos[x-2 y]+Cos[2 x-y]+Cos[2 x+y]+Cos[x+2 y]);
exact[s_,x_,y_] := (1-s)(x^2+y^2)^2 + s (x^2+y^2);
pl[s_]:=Plot[{spline[s,x,0],exact[s,x,0]},{x,-Pi,Pi},
    	DisplayFunction->Identity,PlotRange->{0,16},
	Frame->True,AspectRatio->1,TextStyle->{FontSize->18},
	PlotLabel->SequenceForm["tension=",s],
	PlotStyle->{{},{GrayLevel[0.5],Dashing[{0.02,0.02}]}}];
Show[GraphicsArray[{Show[pl[#]]& /@ {0,0.3},Show[pl[#]]& /@ {0.7,1}}]];
Display["junk_ma.eps",%,"EPS",ImageSize->{576,864}];

