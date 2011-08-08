m[p_]:=1/12(1+p)(2+p);
m1[p_]:=1/1680(1+p)(2+p)(3+p)(4+p);
m2[p_]:=1/420 (4-p)(2+p)(3+p)(4+p);
f3[p_,z_]:= m[p]  z + (1-m[p]-m[-p]) + m[-p] /z;
f5[p_,z_]:= m1[p]  z^2 + m2[p] z + 
	   (1-m1[p]-m1[-p]-m2[p]-m2[-p]) + 
	   m2[-p] /z + m1[-p]/z^2;
a3[p_,z_]:=f3[p,z] /f3[p,1/z];
a5[p_,z_]:=f5[p,z] /f5[p,1/z];
<< PlotLegends`;
phase[p_]:= Plot[{w p,Arg[a3[p,Exp[I w]]],Arg[a5[p,Exp[I w]]]},{w,0,Pi},
    PlotStyle->{{GrayLevel[0]},{Dashing[{0.03}]},{Dashing[{0.01}]}},
    Frame->True,FrameLabel->{"Frequency",None,None,None},
    PlotLabel->"Phase",
    PlotLegend->{"Exact","3-point","5-point"},LegendPosition->{-0.85,0.1},
    LegendSize->{0.65,0.4}];
phase[1/2];
phase[4/5];
ga = Show[GraphicsArray[{%%,%}]];
Export["junk_ma.eps",ga,"EPS",ImageSize->864];

