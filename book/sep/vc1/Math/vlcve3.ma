tzo[t_,a_,g_]:= Sqrt[t^2 +a^2 (1-g^2)];
xzo[a_,g_]:=a (1-g^2);
f[t_,p_,h_]:=(t + Sqrt[t^2 + 4 h^2 p^2])/2;
t2[t_,p_,g_,h_]:=t f[t,p,h] + (1-1/g^2) h^2;
tp[t_,p_,s_,h_]:=t2[t,p,s,h]/Sqrt[f[t2[t,p,s,h],p f[t,p,h],h]];
xp[t_,p_,g_,h_]:=h^2 p/f[t,p,h] * (1 - f[t,p,h]^2/f[t2[t,p,g,h],p f[t,p,h],h]);
tco[t_,a_,g_,h_]:=tp[tzo[t,a,g],a/tzo[t,a,g],g,h];
xco[t_,a_,g_,h_]:=xzo[a,g] + xp[tzo[t,a,g],a/tzo[t,a,g],g,h];
Pop[t_,g_]:=ParametricPlot[{xco[t,a,g,1],g tco[t,a,g,1]}, {a,-g^2*t,g^2*t}];
Show[Pop[0.4,1.2],Pop[0.8,1.2],Pop[1.2,1.2],Pop[1.6,1.2],Pop[2,1.2],
	Pop[2.4,1.2],Pop[2.8,1.2],Pop[3.2,1.2],Pop[3.6,1.2],
	DefaultFont->{"Times",12},PlotRange->{{-2.5,2.5},{0.,4.5}},
	AspectRatio->1/GoldenRatio,Frame->True,Axes->False,
	FrameTicks->{{-2,0,2},{0,2,4},None,None},
	FrameLabel->{"x (km)","","",""}];
Display["junk_ma.eps", %,"EPS"];
