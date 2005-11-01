1+1;
t[x_,v_]:=Sqrt[1 + x^2/(1-v^2)];
xr[y_,v_]:=y (1-v^2);
tr[y_,v_]:=Sqrt[1 + y^2 (1-v^2)];
vr[y_]:=ParametricPlot[{xr[y,v],-tr[y,v]},{v,0,1.4},
	PlotStyle->Dashing[{0.01,0.01}]];
Show[vr[-1],vr[-0.66],vr[0.66],vr[1]];
tt[y_,v_]:=Plot[-t[x,v],{x,-y*(1-v^2)*1.1,y*(1-v^2)*1.1},
	PlotStyle->Thickness[0.01]];
trt[v_]:=tt[1.,v];
Show[trt[0],trt[0.6],trt[0.8],trt[1.2],trt[1.4]];
Show[%6,%9,PlotRange->{{-1.2,1.2},{0,-1.6}},
	Ticks->{{-1.2,-0.8,-0.4,0.4,0.8,1.2},
	{{-0.4,"0.4"},{-0.8,"0.8"},{-1.2,"1.2"},{-1.6,"1.6"}}},
	AxesLabel->{"x",""}];
pt[x_,v_]:=0.5*x/Sqrt[1 - 0.25 v^2];
pxr[y_,v_]:=y (1 - 0.25 v^2);
ptr[y_,v_]:=0.5 y Sqrt[1 - 0.25 v^2];
pvr[y_]:=ParametricPlot[{pxr[y,v],-ptr[y,v]},{v,0,2},
	PlotStyle->Dashing[{0.01,0.01}]];
Show[pvr[0.25],pvr[0.5],pvr[0.75],pvr[1]];
ptt[y_,v_]:=Plot[-pt[x,v],{x,0,y*(1-0.25 v^2)*1.1},
	PlotStyle->Thickness[0.01]];
prt[v_]:=ptt[1.,v];
Show[prt[0],prt[0.8],prt[1.2],prt[1.6]];
Show[%15,%18,PlotRange->{{0,1.2},{0,-0.6}},
	Ticks->{{0.4,0.8,1.2},
	{{-0.2,"0.2"},{-0.4,"0.4"},{-0.6,"0.6"}}},
	AxesLabel->{"x",""}];
Display["junk.eps",%10,"EPS"];

