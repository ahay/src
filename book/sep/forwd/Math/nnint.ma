nn[x_ /; Abs[x] <= 1/2]:= 1;
nn[x_ /; Abs[x] > 1/2] := 0;
nnf[w_]:=Integrate[Cos[w x],{x,-1/2,1/2}];
Plot[nn[x], {x, -2, 2}, PlotStyle -> {Thickness[0.01]}, 
    	Frame -> True, 
   	FrameLabel -> {None, None, 
        "Nearest Neighbor Interpolation: B-0", None}];
Plot[nnf[w],{w,-2Pi,2Pi},PlotStyle->{Thickness[0.01]},
	Frame->True,
	FrameLabel->{None,None,"Spectrum",None}];
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];
