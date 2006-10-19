<<Graphics`Legend`;
p[y_, a_] := 180/Pi*(2 ArcCot[(1 - 4 a^2) Csc[(Pi y)/180]^2/a^3]/a + 
	(1 - Cos[(Pi y)/180]));
Plot[{p[y,2/3],p[y,4/3], p[y,1]},{y,0,90},
	FrameTicks->{Automatic,Automatic, None, None}, PlotRange->All, 
  FrameLabel->{"dip angle",None,None,None} ,GridLines->Automatic, Frame->True,
  PlotStyle->{{AbsoluteThickness[3],Dashing[{0.03}]},{AbsoluteThickness[3],
        Dashing[{0.06}]},AbsoluteThickness[3]},
  PlotLegend->{"a=4/3","a=2/3","a=1"},LegendPosition->{-0.5,0.1}, 
  LegendSize->{0.6,0.4}, PlotLabel->"Phase Error (in degrees)"];
Display["junk_ma.eps",%,"EPS"];
