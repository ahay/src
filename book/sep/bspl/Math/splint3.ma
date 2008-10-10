z1[x_]:=1/6(4-6 x^2+3 x^3);
z2[x_]:=1/6(2-x)^3;
spl[x_ /; Abs[x] <= 1]:= z1[Abs[x]];
spl[x_ /; Abs[x] > 1 && Abs[x] <= 2] := z2[Abs[x]];
spl[x_ /; Abs[x] > 2] := 0;
splf[w_]:=
  2(Integrate[Cos[w x] z1[x],{x,0,1}]+Integrate[Cos[w x] z2[x],{x,1,2}]);
Plot[spl[x], {x, -6, 6}, PlotStyle -> {Thickness[0.01]}, 
    	Frame -> True, 
    	FrameLabel -> {None, None, "B-spline Representation: B-3",None}];
Plot[splf[w],{w,-2Pi,2Pi},PlotStyle->{Thickness[0.01]},
	Frame->True,
  	FrameLabel->{None,None,"Spectrum",None}];
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];
