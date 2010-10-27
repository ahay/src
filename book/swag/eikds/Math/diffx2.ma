<< PlotLegends`;
LinVeikdsF[v00_, dvdx_, x_, z_, dy_] := Block[{b, tt},
  a =  -(2 dvdx (2 v00 + dvdx x) (x^2 + z^2) (2 v00^2 (2 x^2 + z^2) + 
         dvdx x^2 (dvdx (x^2 + z^2) + 
            Sqrt[(x^2 + z^2) ((2 v00 + dvdx x)^2 + dvdx^2 z^2)]) + 
         2 v00 x (dvdx (2 x^2 + z^2) + 
            Sqrt[(x^2 + z^2) ((2 v00 + dvdx x)^2 + 
               dvdx^2 z^2)])))/(v00 (v00 + 
        dvdx x) Sqrt[(x^2 + z^2) ((2 v00 + dvdx x)^2 + 
         dvdx^2 z^2)] (2 v00 x + dvdx (x - z) (x + z) + 
        Sqrt[(x^2 + z^2) ((2 v00 + dvdx x)^2 + 
           dvdx^2 z^2)]) (2 v00 x + dvdx (x^2 + z^2) + 
        Sqrt[(x^2 + z^2) ((2 v00 + dvdx x)^2 + dvdx^2 z^2)]));
  zc = x;
       xc = -z;
  v = v00 + dvdx zc;
  
  z0 = v00/dvdx;
       zc1 = zc + z0;
       x0 = -(xc^2 + zc1^2 - z0^2)/(2 xc);
       xc1 = xc + x0;
  
  r = Sqrt[xc1^2 + zc1^2];
  
  c0 = x0/r;
       c = xc1/r;
  b = Simplify[Log[(v (1 + c0))/(v00 (1 + c))]/dvdx];
  tt = b + dy  a;
  Return[{b, tt, a}]
  ];
aa = ShowLegend[
  ContourPlot[
   Abs[LinVeikdsF[2000, 0.3, xx, zz, 100][[2]] - 
     LinVeikdsF[2030, 0.3, xx, zz, 100][[1]]], {xx, 0, 5000}, {zz, 0, 
    5000}, Contours -> 20, FrameLabel -> {{z[m],}, {x[m], "a"}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow"], {ColorData["Rainbow"][1 - #1] &, 13, 
   "0.0005", "0", LegendPosition -> {1.1, -0.7}, 
   LegendSize -> {0.4, 1.5}, LegendLabel -> Style["\[Tau] (s)", 13]}];
bb = ShowLegend[
  ContourPlot[
   Abs[LinVeikdsF[2000, 0.3, xx, zz, 200][[2]] - 
     LinVeikdsF[2060, 0.3, xx, zz, 200][[1]]], {xx, 0, 5000}, {zz, 0, 
    5000}, Contours -> 20, FrameLabel -> {{z[m],}, {x[m], "b"}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow"], {ColorData["Rainbow"][1 - #1] &, 13, 
   "0.002", "0", LegendPosition -> {1.1, -0.7}, 
   LegendSize -> {0.4, 1.5}, LegendLabel -> Style["\[Tau] (s)", 13]}];
ff = GraphicsRow[{aa, bb}, ImageSize -> {900, 400}, AspectRatio -> 1, 
  Spacings -> 0];
Display["junk_ma.eps", ff, "EPS"];
