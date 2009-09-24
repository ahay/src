z[x_, h_, t_] := -Sqrt[(-h + t)*(h + t)*(t - x)*(t + x)]/t;
hyper[a_, h_, t_] := Plot[-Sqrt[(-h^2 + t^2)*Cos[a]^2 + (x - t*Sin[a])^2], 
 {x, -2*t, 2*t}, PlotRange -> {{-1.5, 1.5}, {0.2, -1}}];
Table[hyper[a, 0.8, 1.2], {a, 0, 2  Pi, Pi/20}];
Plot[z[x, 0.8, 1.2], {x, -1.2, 1.2}, PlotRange -> All, 
  PlotStyle -> {Thickness[0.01]}, AspectRatio -> Automatic];
Plot[-Sqrt[1.2^2 - 0.8^2]*Sqrt[(0.8 - x)*(0.8 + x)]/0.8,{x, -0.8, 0.8}, 
  PlotRange -> All, PlotStyle -> {Thickness[0.01]}, AspectRatio -> Automatic];
Show[{%%%, %%, %, 
    Graphics[{Text["-h", {-0.8, 0.1}], Text["h", {0.8, 0.1}],
    Line[{{-1.5,0},{-1.5,-2.1}}],Line[{{0,0},{0,-2.1}}]}]}, 
    AspectRatio -> Automatic, AxesLabel -> {"y", None}, 
    Axes->{Automatic,None},
    Ticks -> None];
Export["junk_ma.eps",%, "EPS"];

