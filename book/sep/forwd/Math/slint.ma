slf[w_,t_]:=NIntegrate[
    Exp[I w x] (1 - Abs[x - t]), {x, t - 1, t + 1}]/((1 - t) + t Exp[I w]);
a = 2 - Sqrt[3];
c = (2+ Sqrt[3])/6;
f[n_]:=Sum[(-1)^k (a z)^k,{k,0,n}] Sum[(-1)^k (a / z)^k,{k,0,n}] /c;
Series[Expand[N[f[100]]],{z,0,100}];
Plot[Sum[SeriesCoefficient[%,k] spl[x-k], {k,-10,10}],{x,-6,6},PlotStyle->{Thickness[0.01]},
	Frame->True,PlotRange->All,
  	FrameLabel->{None,None,"B-spline interpolator: B-3",None}];
Plot[slf[w,1/2 (1 - Sqrt[3]/3)],{w,-2Pi,2Pi},PlotStyle->{Thickness[0.01]},
	Frame->True,
  	FrameLabel->{None,None,"Spectrum",None}]; 	
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];
