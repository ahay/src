z1[x_]:=(2416 - 1680 x^2 + 560 x^4 - 140 x^6 + 35 x^7)/5040;
z2[x_]:=(2472 - 392 x - 504 x^2 - 1960 x^3 + 2520 x^4 - 1176 x^5 + 252 x^6 - 21 x^7)/5040;
z3[x_]:=(-1112 + 12152 x - 19320 x^2 + 13720 x^3 - 5320 x^4 + 1176 x^5 -
 140 x^6 + 7 x^7)/5040;
z4[x_]:=(4-x)^7/5040;
spl[x_ /; Abs[x] <= 1]:= z1[Abs[x]];
spl[x_ /; Abs[x] > 1 && Abs[x] <= 2] := z2[Abs[x]];
spl[x_ /; Abs[x] > 2 && Abs[x] <= 3] := z3[Abs[x]];
spl[x_ /; Abs[x] > 3 && Abs[x] <= 4] := z4[Abs[x]];
spl[x_ /; Abs[x] > 4] := 0;
 Simplify[2(Integrate[Cos[w x] z1[x],{x,0,1}]+
            Integrate[Cos[w x] z2[x],{x,1,2}]+
            Integrate[Cos[w x] z3[x],{x,2,3}]+
	    Integrate[Cos[w x] z4[x],{x,3,4}])];
Plot[spl[x], {x, -6, 6}, PlotStyle -> {Thickness[0.01]}, 
    	Frame -> True, PlotRange->All,
    	FrameLabel -> {None, None, "B-spline Representation: B-7",None}];
Plot[%%,{w,-2Pi,2Pi},PlotStyle->{Thickness[0.01]},
	Frame->True, PlotRange->All,
  	FrameLabel->{None,None,"Spectrum",None}];
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];
