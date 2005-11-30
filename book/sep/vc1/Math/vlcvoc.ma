f[t_,p_,h_]:=(t + Sqrt[t^2 + 4 h^2 p^2])/2;
t2[t_,p_,g_,h_]:=t f[t,p,h] + (1-1/g^2) h^2;
t0[t_,p_,g_,h_]:=f[t2[t,p,g,h],p f[t,p,h],h];
x0[t_,p_,h_]:=h^2 p/ f[t,p,h];
xr[z_,t_,p_,g_,h_]:=x0[t,p,h] - t0[t,p,g,h]/(p f[t,p,h]) * (1 - z/t2[t,p,g,h]);
tp[t_,p_,s_,h_]:=t2[t,p,s,h]/Sqrt[f[t2[t,p,s,h],p f[t,p,h],h]];
xp[t_,p_,g_,h_]:=h^2 p/f[t,p,h] * (1 - f[t,p,h]^2/f[t2[t,p,g,h],p f[t,p,h],h]);
smile[t_,g_,h_]:=ParametricPlot[{x0[t,p,h],-Sqrt[t2[t,p,g,h]]},{p,-2.1,2.1}];
smile[1,1.2,1];
smile[1,1,1];
smile[1,0.8,1];
triangle[t_,g_,h_]:=ParametricPlot[{xp[t,p,g,h],-tp[t,p,g,h]},{p,-10,10},
	PlotStyle->Thickness[0.01]];
triangle[1,1.2,1];
triangle[1,0.8,1];
ray[p_,t_,g_,h_]:=ParametricPlot[{xr[z,t,p,g,h],-Sqrt[z]},
	{z,tp[t,p,g,h]^2,t2[t,p,g,h]}, PlotStyle->Dashing[{0.01,0.01}]];
Show[	ray[-2,1,1.2,1],ray[-1.75,1,1.2,1],ray[-1.5,1,1.2,1],
	ray[-1.25,1,1.2,1],ray[-1,1,1.2,1],ray[-0.75,1,1.2,1],
	ray[-0.5,1,1.2,1],ray[-0.25,1,1.2,1],ray[2,1,1.2,1],
	ray[1.75,1,1.2,1],ray[1.5,1,1.2,1],ray[1.25,1,1.2,1],
	ray[1,1,1.2,1],ray[0.75,1,1.2,1],ray[0.5,1,1.2,1],ray[0.25,1,1.2,1]];
Show[	ray[-2,1,1,1],ray[-1.75,1,1,1],ray[-1.5,1,1,1],
	ray[-1.25,1,1,1],ray[-1,1,1,1],ray[-0.75,1,1,1],
	ray[-0.5,1,1,1],ray[-0.25,1,1,1],ray[2,1,1,1],
	ray[1.75,1,1,1],ray[1.5,1,1,1],ray[1.25,1,1,1],
	ray[1,1,1,1],ray[0.75,1,1,1],ray[0.5,1,1,1],ray[0.25,1,1,1]];
Show[ray[-2,1,0.8,1],ray[-1.75,1,0.8,1],ray[-1.5,1,0.8,1],
	ray[-1.25,1,0.8,1],ray[-1,1,0.8,1],ray[-0.75,1,0.8,1],
	ray[-0.5,1,0.8,1],ray[-0.25,1,0.8,1],ray[2,1,0.8,1],
	ray[1.75,1,0.8,1],ray[1.5,1,0.8,1],ray[1.25,1,0.8,1],
	ray[1,1,0.8,1],ray[0.75,1,0.8,1],ray[0.5,1,0.8,1],ray[0.25,1,0.8,1]];
Show[%9,%13,%16,PlotRange->{{-0.85,0.85},{-0.95,-1.9}},
	Frame->True,Axes->False,AspectRatio->1/2,
	FrameTicks->{{-0.5,0,0.5},{{-1,"1"},{-1.4,"1.4"},{-1.8,"1.8"}}},
	FrameLabel->{"midpoint (km)","","",""}];
Show[%10,%17,PlotRange->{{-0.85,0.85},{-0.75,-1.7}},
	Frame->True,Axes->False,AspectRatio->1/2,
	FrameTicks->{{-0.5,0,0.5},{{-0.8,"0.8"},{-1.2,"1.2"},{-1.6,"1.6"}}},
	FrameLabel->{"midpoint (km)","","",""}];
Show[%11,%14,%18,PlotRange->{{-0.85,0.85},{-0.55,-1.5}},
	Frame->True,Axes->False,AspectRatio->1/2,
	FrameTicks->{{-0.5,0,0.5},{{-0.6,"0.6"},{-1,"1"},{-1.4,"1.4"}}},
	FrameLabel->{"midpoint (km)","","",""}];
Show[GraphicsArray[{{%19},{%20},{%21}}]];
Display["junk_ma.eps",%,"EPS",ImageSize->{594,783}];