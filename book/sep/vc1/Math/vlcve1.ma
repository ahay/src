tn2[t_,h_,g_]:=g^2 t^2 - h^2 (1-g^2);
zo[t_,h_,g_,x_]:=Sqrt[tn2[t,h,g]+g^2*x^2/(1-g^2)];
zop[t_,h_,g_,a_]:=Sqrt[tn2[t,h,g]+a^2 (1-g^2) g^2];
xop[a_,g_]:=a (1-g^2);
Pop[t_,g_]:=ParametricPlot[{xop[a,g],-zop[t,1.,g,a]},
	{a,-g^2*t,g^2*t}];
Show[Pop[0.4,1.2],Pop[0.8,1.2],Pop[1.2,1.2],Pop[1.6,1.2],Pop[2.,1.2],
	Pop[2.4,1.2],Pop[2.8,1.2],Pop[3.2,1.2],Pop[3.6,1.2],
	DefaultFont->{"Times",12},PlotRange->{{-2.5,2.5},{0.,-4.5}},
	AspectRatio->1/GoldenRatio,Frame->True,Axes->False,
	FrameTicks->{{-2,0,2},{0,{-2,"2"},{-4,"4"}},None,None},
	FrameLabel->{"Midpoint (km)","Pseudo-depth (km)","",""}];
Display["junk_ma.eps", %, "EPS"];
