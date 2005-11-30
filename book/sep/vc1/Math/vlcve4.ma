tzo[t_,a_,g_]:= Sqrt[t^2 +a^2 (1-g^2)];
xzo[a_,g_]:=a (1-g^2);
f[t_,p_,h_]:=(t + Sqrt[t^2 + 4 h^2 p^2])/2;
t2[t_,p_,g_,h_]:=t f[t,p,h] + (1-1/g^2) h^2;
tp[t_,p_,s_,h_]:=t2[t,p,s,h]/Sqrt[f[t2[t,p,s,h],p f[t,p,h],h]];
xp[t_,p_,g_,h_]:=h^2 p/f[t,p,h] * (1 - f[t,p,h]^2/f[t2[t,p,g,h],p f[t,p,h],h]);
tco[t_,a_,g_,h_]:=tp[tzo[t,a,g],a/tzo[t,a,g],g,h];
xco[t_,a_,g_,h_]:=xzo[a,g] + xp[tzo[t,a,g],a/tzo[t,a,g],g,h];
Push[t_,g_]:=ParametricPlot[{xco[t,a,g,1], g tco[t,a,g,1]}, {a,-3.5*t,3.5*t}];
Show[Push[0.8,0.8],Push[1.2,0.8],Push[1.6,0.8],Push[2.,0.8],
	Push[2.4,0.8],Push[2.8,0.8],Push[3.2,0.8],Push[3.6,0.8],
	DefaultFont->{"Times",12},PlotRange->{{-4.75,4.75},{0.,7.}},
	AspectRatio->1/GoldenRatio,Frame->True,Axes->False,
	FrameTicks->{{-4,0,4},{0,4},None,None},
	FrameLabel->{"x (km)","","",""}];
Display["junk_ma.eps", %,"EPS"];
