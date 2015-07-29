<< Miscellaneous`WorldPlot`;
earth = WorldPlot[{World, RandomGrays}];
RSE = ToMinutes[{{-33, -52}, {151, 16}}];
RSG = ToMinutes[{{-5, -57}, {-49, -39}}];
RSF = (RSE + RSG)/2;
Show[{earth, WorldGraphics[{RGBColor[1, 1, 0], Line[{RSE, RSF, RSG}]}], 
    WorldGraphics[{PointSize[0.03], RGBColor[1, 0, 0], Point[RSE], Point[RSF],
         Point[RSG]}], 
    WorldGraphics[{RGBColor[0, 1, 0], 
        Text[FontForm["RSE", {"Helvetica-Bold", 12}], RSE + {600, 100}], 
        Text[FontForm["RSF", {"Helvetica-Bold", 12}], RSF + {600, 100}], 
        Text[FontForm["RSG", {"Helvetica-Bold", 12}], RSG + {600, 100}]}]},
    WorldGridStyle -> {Thickness[0.001], RGBColor[1, 1, 1]}];
Display["junk_ma.eps", Graphics[%], "EPS"];
