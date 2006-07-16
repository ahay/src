x[j_, q_, t_] := (j Cos[t])/Sqrt[q];
z[j_, q_, t_] := - (j Sin[t])/q - j (1/q - 1) - 1; 
g[q_, a_] := ParametricPlot[{x[a, q, t], z[a, q, t]}, {t, 0, 2 Pi},  
        DisplayFunction -> Identity];
gplot[q_] := Show[g[q, 0.001], g[q, 0.2], g[q, 0.4], g[q, 0.6], g[q, 0.8], 
        PlotRange -> {{-1.2, 1.2}, {0, -4}}, Frame -> True, 
        Axes -> False, AspectRatio -> Automatic];
Show[GraphicsArray[{gplot[1/2], gplot[3/2]}]];
Display["junk_ma.eps", %, "EPS"];

