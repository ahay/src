major = {{-0.4,{RGBColor[1,0,0]}},{-0.8,{RGBColor[1,0,0]}},
	{-1.2,{RGBColor[1,0,0]}},{-1.6,{RGBColor[1,0,0]}}};
f[t_,p_,h_]:=(t + Sqrt[t^2 + 4 h^2 p^2])/2;
t2[t_,p_,g_,h_]:=t f[t,p,h] + (1-1/g^2) h^2;
tp[t_,p_,s_,h_]:=t2[t,p,s,h]/Sqrt[f[t2[t,p,s,h],p f[t,p,h],h]];
xp[t_,p_,g_,h_]:=h^2 p/f[t,p,h] * (1 - f[t,p,h]^2/f[t2[t,p,g,h],p f[t,p,h],h]);
triangle[t_,g_,h_]:=ParametricPlot[{xp[t,p,g,h],-tp[t,p,g,h]},{p,-10,10},
	PlotStyle->{RGBColor[1,130/255,3/255],Thickness[0.01]}];
Show[triangle[0.5,1.333,1],triangle[1,1.333,1],triangle[1.5,1.333,1]];
Show[%,PlotRange->{{-0.5,0.5},{0,-2}},Frame->True,Axes->False,
	FrameTicks->{{-0.4,-0.2,0,0.2,0.4},
	{0,{-0.4,"0.4"},{-0.8,"0.8"},{-1.2,"1.2"},{-1.6,"1.6"}},None,None},
	GridLines->{None,major},	
	AspectRatio->1,DefaultColor->GrayLevel[1]];
Show[triangle[0.75,0.833,1],triangle[1.25,0.833,1],triangle[1.75,0.833,1]];
Show[%,PlotRange->{{-0.5,0.5},{0,-2}},Frame->True,Axes->False,
	FrameTicks->{{-0.4,-0.2,0,0.2,0.4},
	{0,{-0.4,"0.4"},{-0.8,"0.8"},{-1.2,"1.2"},{-1.6,"1.6"}},None,None},
	GridLines->{None,major},
	AspectRatio->1,DefaultColor->GrayLevel[1]];
Show[GraphicsArray[{%8,%10}]];
Display["junk_ma.eps",%,"EPS"];
