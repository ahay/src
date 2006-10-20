<<Graphics`Legend`;
r[k_,b_]:=-Sin[k]^2/(1-4 b Sin[k]^2);
Plot[{	4 r[k/2,0],
	4 r[k/2,1/12], 
	4 r[k/2,0.122996],	
	-k^2, 
	4 r[k/2,0.148679]}, 
	{k,0,Pi},PlotRange->All,
	FrameTicks->{Automatic,Automatic, None, None}, 
	PlotRange->All, 
  	FrameLabel->{"wavenumber",None,None,None},
	GridLines->Automatic, 
	Frame->True,
  	PlotStyle->{
	{AbsoluteThickness[2],Dashing[{0.01}]},
	{AbsoluteThickness[2],Dashing[{0.03,0.03}]},
	{AbsoluteThickness[2],Dashing[{0.03,0.03,0.01}]},
	AbsoluteThickness[2],
	{AbsoluteThickness[2],Dashing[{0.06}]}},
  	PlotLegend->{"beta=0","beta=1/12","beta=1/8.13","exact","beta=1/6.73"},
	LegendPosition->{-0.7,-0.4}, 
  	LegendSize->{0.9,0.5}, 
	PlotLabel->"Second Derivative Approximations"];
Display["junk_ma.eps",%,"EPS"];
