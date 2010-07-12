Vp[\[Theta]_, v_, vv_, \[Eta]_] := 
 Block[{}, 
  vp = 1/4 (2 v^2 (1 + 2 \[Eta]) Sin[\[Theta]]^2 + 
      2 Cos[\[Theta]]^2 vv^2 + Sqrt[
      4 v^4 (1 + 2 \[Eta])^2 Sin[\[Theta]]^4 + 
       2 v^2 (1 - 2 \[Eta]) Sin[2 \[Theta]]^2 vv^2 + 
       4 Cos[\[Theta]]^4 vv^4]); 
  vs = 1/4 (2 v^2 (1 + 2 \[Eta]) Sin[\[Theta]]^2 + 
      2 Cos[\[Theta]]^2 vv^2 - Sqrt[
      4 v^4 (1 + 2 \[Eta])^2 Sin[\[Theta]]^4 + 
       2 v^2 (1 - 2 \[Eta]) Sin[2 \[Theta]]^2 vv^2 + 
       4 Cos[\[Theta]]^4 vv^4]);
  Return[{Sqrt[vp], Sqrt[vs]}]];
Kz[v_, vv_, \[Eta]_, km_, freq_, \[Theta]_] :=
 
 Block[{stheta, \[Omega]},
  \[Omega] = 2 Pi freq;
  stheta = 1.0/Vp[\[Theta], v, vv, \[Eta]][[1]];
  k\[Lambda] = 
   Tan[ \[Theta]] Sqrt[(2 \[Omega] stheta Cos [\[Theta]])^2 - km^2];
  kz = Sqrt[\[Omega]^2 stheta^2 - (km - k\[Lambda])^2] + 
    Sqrt[\[Omega]^2 stheta^2 - (km + k\[Lambda])^2];
  Return[kz]];
<< PlotLegends`
a = ShowLegend[
   ContourPlot[
     Kz[2, 1.8, 0.3, km, 25, theta*Pi/180.], 
    {km, -90, 90}, {theta, 0, 30}, 
     Contours -> 20, FrameLabel -> {Subscript[k, mx], \[Theta]^o}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14],
     ColorFunction -> "Rainbow"], {ColorData["Rainbow"][1 - #1] &, 13, 
     "150", "0", LegendPosition -> {1.1, -0.7}, LegendSize -> 1.5,  
     LegendLabel -> Style["\!\(\*SubscriptBox[\"k\", \"z\"]\)", 14]}];
b = ShowLegend[
   ContourPlot[
     -Kz[2, 1.8, 0.3, km, 25, theta*Pi/180.]+Kz[1.8, 1.8, 0.0, km, 25, theta*Pi/180.], 
    {km, -90, 90}, {theta, 0, 30}, 
     Contours -> 20, FrameLabel -> {Subscript[k, mx], \[Theta]^o}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14],
     ColorFunction -> "Rainbow"], {ColorData["Rainbow"][1 - #1] &, 13, 
     "150", "0", LegendPosition -> {1.1, -0.7}, LegendSize -> 1.5,  
     LegendLabel -> Style["\!\(\*SubscriptBox[\"k\", \"z\"]\)", 14]}];
mm = GraphicsRow[{a, b}, ImageSize -> {600, 320}, AspectRatio -> 1, 
    Spacings -> 0];
Display["junk_ma.eps", %, "EPS"];
