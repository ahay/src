sinc[x_]:= Sin[Pi x]/(Pi x);
nn[x_ /; Abs[x] <= Pi]:= 1;
nn[x_ /; Abs[x] > Pi] := 0;
Plot[sinc[x], {x, -6, 6}, PlotStyle -> {Thickness[0.01]}, 
    	Frame -> True, PlotRange->All,
   	FrameLabel -> {None, None, 
        "Sinc Interpolation", None}];
Plot[nn[w],{w,-2Pi,2Pi},PlotStyle->{Thickness[0.01]},
	Frame->True,
	FrameLabel->{None,None,"Spectrum",None}];
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];
