<<Calculus`Pade`;
<<Graphics`Legend`;
f[a_,k_]:=Exp[- a k^2];
pd[m_,n_,a_,k_]:=Pade[f[a,k],{k,0,m,n}];
pd[4,0,2/3,k]-f[2/3,k];
pd[4,0,4/3,k]-f[4/3,k];
Plot[{%%,%},{k,0,Pi},FrameTicks->{Automatic,Automatic, None, None}, 
  PlotRange->All, FrameLabel->{"wavenumber",None,None,None} ,
  GridLines->Automatic, Frame->True,
  PlotStyle->{{AbsoluteThickness[3],Dashing[{0.05}]},AbsoluteThickness[3]},
  PlotLegend->{"a=2/3","a=4/3"},LegendPosition->{-0.5,0}, 
  LegendSize->{0.4,0.3}, PlotLabel->"Explicit Error"];
pd[2,2,2/3,k]-f[2/3,k];
pd[2,2,4/3,k]-f[4/3,k];
Plot[{%%,%},{k,0,Pi},FrameTicks->{Automatic,Automatic, None, None}, 
  PlotRange->All, FrameLabel->{"wavenumber",None,None,None} ,
  GridLines->Automatic, Frame->True,
  PlotStyle->{{AbsoluteThickness[3],Dashing[{0.05}]},AbsoluteThickness[3]},
  PlotLegend->{"a=2/3","a=4/3"},LegendPosition->{-0.5,0}, 
  LegendSize->{0.4,0.3}, PlotLabel->"Implicit Error"];
Show[GraphicsArray[{%7,%10}]];
Display["junk_ma.eps",%,"EPS",ImageSize->864];

