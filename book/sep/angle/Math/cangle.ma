heops[x_, t_, h_] := 0.5*(Sqrt[t^2 + (x - h)^2] + Sqrt[t^2 + (x + h)^2]);
piramid = Plot3D[-heops[x, 0.5, h], {x, -4, 4}, {h, -4, 4}, 
        PlotRange -> {{-4, 4}, {-4, 4}, {-4, 0}}, Boxed -> False, 
        AxesLabel -> {"Midpoint", "Offset", "Time"}, 
        PlotPoints -> {100, 100}, Mesh -> False, 
        AxesLabel -> {"Offset", "Midpoint", "Time"}, 
        DefaultFont -> {"Times-Roman", 12}, 
        Ticks -> {Automatic, 
        Automatic, {0, {-2, "2"}, {-4, "4"}}}];	
tq[t_, a_, g_] := t Cos[a] Cos[g]/(Cos[a]^2 - Sin[g]^2);
xq[t_, a_, g_] := t Cos[a] Sin[a]/(Cos[a]^2 - Sin[g]^2);
hq[t_, a_, g_] := t Cos[g] Sin[g]/(Cos[a]^2 - Sin[g]^2);
cg1[g_] := ParametricPlot3D[{xq[0.5, a, Pi/2g], hq[0.5, a, Pi/2g], 
         0.1 - tq[0.5, a, Pi/2g]}, 
	{a,-Pi/2(0.9999 - Abs[g]), Pi/2(0.9999 - Abs[g])}, 
        PlotRange -> {{-4, 4}, {-4, 4}, {-4, 0}}, Boxed -> False, 
        PlotPoints -> 100, Axes -> False, DisplayFunction -> Identity];
cangle = Table[cg1[g], {g, -1, 1, 0.05}];
Show[Graphics3D[piramid], cangle];
Display["junk_ma.eps", %, "EPS"];


