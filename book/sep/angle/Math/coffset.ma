heops[x_, t_, h_] := 0.5*(Sqrt[t^2 + (x - h)^2] + Sqrt[t^2 + (x + h)^2]);
piramid = Plot3D[-heops[x, 0.5, h], {x, -4, 4}, {h, -4, 4}, 
        PlotRange -> {{-4, 4}, {-4, 4}, {-4, 0}}, Boxed -> False, 
        AxesLabel -> {"Midpoint", "Offset", "Time"}, 
        PlotPoints -> {100, 100}, Mesh -> False, 
        AxesLabel -> {"Offset", "Midpoint", "Time"}, 
        DefaultFont -> {"Times-Roman", 12}, 
        Ticks -> {Automatic, 
        Automatic, {0, {-2, "2"}, {-4, "4"}}}];	
co[h_] :=  ParametricPlot3D[{x, h, -heops[x, 0.5, h]}, {x, -4, 4}, 
        Boxed -> False, Axes -> False, 
        PlotRange -> {{-4, 4}, {-4, 4}, {-4, 0}}, 
        PlotPoints -> 100, DisplayFunction -> Identity];	
coffset = Table[co[h], {h, -4, 4, 0.25}];
Show[Graphics3D[piramid], coffset];
Display["junk_ma.eps", %, "EPS"];


