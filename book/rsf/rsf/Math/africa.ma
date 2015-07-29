<<Miscellaneous`WorldPlot`;
map = WorldPlot[Africa];
RSE=ToMinutes[{{-33,-52},{151,16}}];
RSG = ToMinutes[{{-5, -57}, {-49, -39}}];
RSF=(RSE+RSG)/2;
Show[{map, 
    WorldGraphics[{PointSize[0.02], RGBColor[1, 0, 0], 
        Point[RSF - {0, 150}]}], 
    WorldGraphics[{RGBColor[0, 1, 0], 
        Text[FontForm["RSF", {"Helvetica-Bold", 12}], RSF + {0, 120}]}]},
    WorldGridStyle -> {Thickness[0.001], RGBColor[1, 1, 1]}];
Display["junk_ma.eps",Graphics[%],"EPS"];
