h[x_]:=-Sqrt[1-0.5 x^2];
Plot[h[x],{x,-0.999999 Sqrt[2],0.999999 Sqrt[2]},
	AspectRatio->Automatic,PlotStyle->Thickness[0.01]];
r[h_]:=(0.01+1-h+Sign[h-1]*
	Sqrt[(0.01-1-h)^2-4 h])/0.1;
pr[s_]:=ParametricPlot[{
	0.1+Sqrt[s]+(r[s]-0.1-Sqrt[s])*t,
	t*h[r[s]]},{t,0,1}];
ps[s_]:=ParametricPlot[{
	0.1-Sqrt[s]+(r[s]-0.1+Sqrt[s])*t,
	t*h[r[s]]},{t,0,1}];
Show[pr[0.16],ps[0.16]];
Show[pr[0.64],ps[0.64]];
Show[pr[1.44],ps[1.44]];
Show[%6,%7,%8,%2,
	AspectRatio->1/(2 Sqrt[2]),Axes->{True,False}];
Display["junk_ma.eps", %, "EPS"];
