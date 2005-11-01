tn2[t_,h_,g_]:=g^2 t^2 - h^2 (1-g^2);
zo[t_,h_,g_,x_]:=Sqrt[tn2[t,h,g]+g^2*x^2/(1-g^2)];
zop[t_,h_,g_,a_]:=Sqrt[tn2[t,h,g]+a^2 (1-g^2) g^2];
xop[a_,g_]:=a (1-g^2);
Push[t_,g_]:=ParametricPlot[{xop[a,g],-zop[t,1.,g,a]},
	{a,-3.5*t,3.5*t}];
Show[Push[0.8,0.8],Push[1.2,0.8],Push[1.6,0.8],Push[2.,0.8],
	Push[2.4,0.8],Push[2.8,0.8],Push[3.2,0.8],Push[3.6,0.8],
	DefaultFont->{"Times",12},PlotRange->{{-4.75,4.75},{0.,-7.}},
	AspectRatio->1/GoldenRatio,Frame->True,Axes->False,
	FrameTicks->{{-4,0,4},{0,{-4,"4"}},None,None},
	FrameLabel->{"Midpoint (km)","Pseudo-depth (km)","",""}];
Display["junk_ma.eps",%,"EPS"];
