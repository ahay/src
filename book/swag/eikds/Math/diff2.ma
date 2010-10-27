LinVeikdsZF[v00_, dvdz_, x_, z_, dy_] := Block[{b, tt},
  a = -(((2 v00 + dvdz z) Sqrt[(dvdz^2 (x^2 + z^2))/(
     dvdz^2 x^2 + (2 v00 + dvdz z)^2)])/(v00 (v00 + dvdz z)));
  b = 1/dvdz  ArcCosh[
     1 + (dvdz z)^2 /(2 v00 (v00 + dvdz z)) (1 + x^2 /z^2)];
  tt = b +   dy  a;
  Return[{b, tt, a}]
  ];
<< PlotLegends`;
aa = ShowLegend[
  ContourPlot[
   Abs[LinVeikdsZF[2000, 0.5, xx, zz, 200][[2]] - 
     LinVeikdsZF[2100, 0.5, xx, zz, 200][[1]]], {xx, 0, 5000}, {zz, 0,
     4000}, Contours -> 20, FrameLabel -> {{z[m],}, {x[m], "a"}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow"], {ColorData["Rainbow"][1 - #1] &, 13, 
   "0.004", "0", LegendPosition -> {1.1, -0.7}, 
   LegendSize -> {0.4, 1.5}, LegendLabel -> Style["\[Tau] (s)", 13]}];
bb = ShowLegend[
  ContourPlot[
   Abs[LinVeikdsZF[2000, 0.7, xx, zz, 200][[2]] - 
     LinVeikdsZF[2140, 0.7, xx, zz, 200][[1]]], {xx, 0, 5000}, {zz, 0,
     4000}, Contours -> 20, FrameLabel -> {{z[m],}, {x[m], "b"}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow"], {ColorData["Rainbow"][1 - #1] &, 13, 
   "0.007", "0", LegendPosition -> {1.1, -0.7}, 
   LegendSize -> {0.4, 1.5}, LegendLabel -> Style["\[Tau] (s)", 13]}];
ff = GraphicsRow[{aa, bb}, ImageSize -> {900, 400}, AspectRatio -> 1, 
  Spacings -> 0];
Display["junk_ma.eps", ff, "EPS"];
