slice[t_, a_, r_] := {
    Sqrt[2]*r*Cos[a]*Sqrt[(t*(-0.04 + t^2))/
    (2*t + Sqrt[2]*Sqrt[0.04 + t^2 - (-0.04 + t^2)*Cos[4*a]])], 
    Sqrt[2]*r*Sqrt[(t*(-0.04 + t^2))/
    (2*t + Sqrt[2]*Sqrt[0.04 + t^2 - (-0.04 + t^2)*Cos[4*a]])]*Sin[a], t};
pyramid = ParametricPlot3D[-slice[t, a, 1], {t, 0.0001, 1}, {a, 0, 2 Pi}, 
    PlotPoints -> 100, ViewPoint -> {1.3, 0, 3}, Mesh -> None, Axes -> False];
cont[t_] := ParametricPlot3D[-slice[t, a, 1.001], {a, 0, 2 Pi}, 
            PlotPoints -> 100, Axes -> False, ViewPoint -> {1.3, 0, 3}, 
            PlotStyle-> {Black}];
hyp[x_, r_] := ParametricPlot3D[{r(x + h), r(x - h), 
	       -Sqrt[0.01 + (x - h)^2] -Sqrt[0.01 + (x + h)^2]}, 
               {h, -0.5, 0.5},ViewPoint -> {1.3, 0, 3},
               PlotStyle -> {Black, Thickness[0.01]}];
Show[Join[{pyramid}, Table[cont[t], {t, 0.2, 1, 0.04}], 
    Table[hyp[x, 1.001], {x, 0, 0.4, 0.1}]], 
    PlotRange -> {Automatic, Automatic, {0, -1}}];
Show[Graphics3D[{
      Line[{{-1, 1, -1}, {1, -1, -1}}], Text["y", {-1, 0.85, -1}], 
      Polygon[{{-1, 1, -1}, {-0.875, 0.925, -1}, {-0.925, 0.875, -1}}], 
      Line[{{2, 0, -1}, {1, -1, -1}}], Text["g", {0.85, 1, -1}],
      Polygon[{{2, 0, -1}, {1.925, -0.125, -1}, {1.875, -0.075, -1}}],
      Line[{{1, 1, -1}, {1, -1, -1}}], Text["s", {-1, -0.85, -1}], 
      Polygon[{{1, 1, -1}, {1.025, 0.9, -1}, {0.975, 0.9, -1}}],
      Line[{{-1, -1, -1}, {1, -1, -1}}], Text["h", {1.9, 0.05, -1}],
      Polygon[{{-1, -1, -1}, {-0.9, -1.025, -1}, {-0.9, -0.975, -1}}]}, 
    ViewPoint -> {1.3, 0, 3}]];
left = Show[{%, %%}, PlotRange -> {All, All, {0, -1}}, Boxed -> False];
Show[Graphics[{
      Arrow[{{1, -1}, {3, 1}}], Text["y", {2.85, 1}], Arrow[{{1, -1}, {2, -2}}], 
      Text["g", {3, -0.85}],
      Arrow[{{1, -1}, {3, -1}}], Text["s", {1.15, 1}], Arrow[{{1, -1}, {1, 1}}], 
      Text["h", {2.05, -1.9}], 
      Line[{{2.8944271909999157, 0}, {1.1055728090000843, 0}}], 
      Line[{{2, -0.8944271909999159}, {2, 0.8944271909999159}}]}], 
  AspectRatio -> Automatic];
ParametricPlot[{slice[1, a, 1][[2]] + 2, -slice[1, a, 1][[1]]}, {a, 0, 2 Pi}, 
  PlotPoints -> 100, Axes -> False];
right = Show[{%%, %}, AspectRatio -> Automatic];
Show[GraphicsArray[{Rasterize[left], Rasterize[right]}]];
Export["junk_ma.eps",%, "EPS"];
