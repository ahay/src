ccin[x_ /; Abs[x] <= 1]:= 3/2 Abs[x]^3 - 5/2 Abs[x]^2 + 1;
ccin[x_ /; Abs[x] > 1 && Abs[x] <= 2]:= 
	-1/2 Abs[x]^3 + 5/2 Abs[x]^2 - 4 Abs[x] + 2;
ccin[x_ /; Abs[x] > 2] := 0;
ccinf[w_]:= 2 Integrate[Cos[w x] (3/2 x^3 - 5/2 x^2 + 1),{x,0,1}] +
	    2 Integrate[Cos[w x] (-1/2 x^3 + 5/2 x^2 - 4 x + 2),{x,1,2}];
Plot[ccin[x], {x, -2, 2}, PlotStyle -> {Thickness[0.01]}, 
    	Frame -> True, 
    	FrameLabel -> {None, None, 
        "Cubic Convolution Interpolation", None}];
Plot[ccinf[w],{w,-2Pi,2Pi},PlotStyle->{Thickness[0.01]},
	Frame->True,
 	FrameLabel->{None,None,"Spectrum",None}];
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];
