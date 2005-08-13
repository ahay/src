err1[a_,p_]:=((-1 - p)*ArcCoth[Cos[a]] + 
  (1 + p)*ArcTanh[(1 + p)/Sqrt[(1 + p)^2 - Sin[a]^2]] + 
  Cos[a] - Sqrt[(1 + p)^2 - Sin[a]^2])/
 ((1 + p)*ArcCoth[Cos[a]] - 
  (1 + p)*ArcTanh[(1 + p)/Sqrt[(1 + p)^2 - Sin[a]^2]] - 
  (1 + p)*Cos[a] + Sqrt[(1 + p)^2 - Sin[a]^2]);
err2[a_,p_]:=(3*(1 + p)*ArcCoth[Cos[a]]*(Cos[a] + 
    Sqrt[(1 + p)^2 - Sin[a]^2]) - 
  3*(1 + p)*ArcCoth[Sqrt[(1 + p)^2 - Sin[a]^2]/(1 + p)]*
   (Cos[a] + Sqrt[(1 + p)^2 - Sin[a]^2]) - 
  (Cos[a] - Sqrt[(1 + p)^2 - Sin[a]^2])*
   (2*Cos[a] + (1 + p)*Cos[a] + 
    Sqrt[(1 + p)^2 - Sin[a]^2] + 
    2*(1 + p)*Sqrt[(1 + p)^2 - Sin[a]^2]))/
 (3*((-1 - p)*ArcCoth[Cos[a]] + 
   (1 + p)*ArcCoth[Sqrt[(1 + p)^2 - Sin[a]^2]/(1 + p)] + 
   (1 + p)*Cos[a] - Sqrt[(1 + p)^2 - Sin[a]^2])*
  (Cos[a] + Sqrt[(1 + p)^2 - Sin[a]^2]));
errplot[p_,label_]:=Plot[{100 Abs[err2[a Pi/180, p/100]], 
		          100 Abs[err1[a Pi/180, p/100]]}, 
			  {a, 0.01,80}, 
  Frame -> True, 
  PlotStyle -> {{Thickness[0.01]}, {Thickness[0.01],Dashing[{0.03}]}}, 
  GridLines -> None, PlotRange -> All, PlotLabel->label,
  FrameLabel -> {"Angle (degrees)","Error (%)"}];
errplot[0.1,"a"];
errplot[1,"b"];
errplot[10,"c"];
Show[GraphicsArray[{{%%%},{%%},{%}}],AspectRatio->3/GoldenRatio];
Display["junk_ma.eps", %, "EPS",ImageSize->300];
