lin[x_ /; Abs[x] <= 1]:= 1-Abs[x];
lin[x_ /; Abs[x] > 1] := 0;
linf[w_]:=2 Integrate[Cos[w x] (1-x),{x,0,1}];
Plot[lin[x], {x, -2, 2}, PlotStyle -> {Thickness[0.01]}, 
    	Frame -> True, 
    	FrameLabel -> {None, None, 
        "Linear Interpolation: B-1", None}];
Plot[linf[w],{w,-2Pi,2Pi},PlotStyle->{Thickness[0.01]},
	Frame->True,
 	FrameLabel->{None,None,"Spectrum",None}];
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];
