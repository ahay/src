ph[o_]:=Sqrt[2] Gamma[1/2-I o/2] /
		Gamma[-I o/2];
ap[o_]:=(-I o)^(1/2);
Plot[Abs[ph[o]],{o,0,10}];
Plot[Abs[ap[o]],{o,0,10},
	PlotStyle->Dashing[{0.01,0.01}]];
Show[%3,%4,
	AxesLabel->{"omega","|F|"}];
Plot[Arg[ph[o]]/Pi,{o,0,10}];
Plot[Arg[ap[o]]/Pi,{o,0,10},
	PlotStyle->Dashing[{0.01,0.01}]];
Show[%7,%6,AxesLabel->{"omega","Arg(F)/Pi"},
	PlotRange->{-0.5,-0.2}];
Show[GraphicsArray[{%5,%8}],AspectRatio->1/2];
Display["junk_ma.eps", %, "EPS", ImageSize->432];
