ymax[x_, h_] := (-1 + x^2 + Sqrt[1 + (2 + 4 h^2) x^2 + x^4])/(2 x);
ray[y_, x_] := Sqrt[1 + x y];
fro[y_, h_] := Sqrt[1/2 (y^2  +  1  - h^2  +  
               Sqrt[(y^2 + 1 - h^2)^2 + 4 h^2])];
plotfro[h_] := Plot[-fro[y, h], {y, -ymax[5.1, h], ymax[5.1, h]}, 
        PlotStyle -> {Thickness[0.01]},  
        DisplayFunction -> Identity];
fronts  = Show[plotfro[0], plotfro[1], plotfro[2], plotfro[3], plotfro[4], 
                plotfro[5]];
plotray[x_] := Plot[-ray[y, x], {y, -1/x, ymax[x, 5]}, 
        PlotStyle -> {Dashing[{0.02}]}];
plot1 = Show[plotray[1], plotray[2], plotray[3], plotray[4], plotray[5], 
        plotray[-1], plotray[-2], plotray[-3], plotray[-4], 
        plotray[-5], fronts, AspectRatio -> 1, PlotRange->{-7,0},
	Frame->True, FrameLabel -> {"midpoint", None},
	PlotLabel->"Diffraction Point",FrameTicks->{Automatic,
	{{-1,"1"},{-2,"2"},{-3,"3"},{-4,"4"},{-5,"5"},{-6,"6"}},
	None,None}]; 
ray[y_, x_] := Sqrt[1 - x y];
plotray1[x_] := Plot[-ray[y, x], {y, -1, 1/x}, 
      PlotStyle -> {Dashing[{0.02}]}]; 
plotray2[x_] := Plot[-ray[y, x], {y, 1/x, 1}, 
      PlotStyle -> {Dashing[{0.02}]}];
fro1[y_, h_] := Sqrt[1/2 (1 + h^2 - y^2 + 
              Sqrt[-4 h^2 + (-1 - h^2 + y^2)^2])];
fro2[y_, h_] := Sqrt[1/2 (1 + h^2 - y^2 - 
              Sqrt[-4 h^2 + (-1 - h^2 + y^2)^2])];
plotfro1[h_] := Plot[-fro1[y, h], {y, - (1 - h), 1 - h}, 
        PlotStyle -> {Thickness[0.01]},
        DisplayFunction -> Identity];
fronts1 = Show[plotfro1[0.8], plotfro1[0.6], plotfro1[0.4], plotfro1[0.2]];
plotfro2[h_] := Plot[-fro2[y, h], {y, - (h - 1), h - 1}, 
        PlotStyle -> {Thickness[0.01]},
        DisplayFunction -> Identity];
fronts2 = Show[plotfro2[1.2], plotfro2[1.4], plotfro2[1.6], plotfro2[1.8]];
plot2  =  Show[plotray1[1], plotray1[0.4], plotray1[0.6], plotray1[0.8], 
        plotray2[-1], plotray2[-0.4], plotray2[-0.6], plotray2[-0.8], 
	fronts1, fronts2, 
	Frame->True, PlotLabel->"Elliptic Reflector",
	FrameLabel -> {"midpoint", None},	
        PlotRange -> {{-1.1, 1.1}, {-1.4,0}},
	FrameTicks -> {{-1,-0.5,0,0.5,1},
	{{-0.2,"0.2"},{-0.4,"0.4"},{-0.6,"0.6"},
	 {-0.8,"0.8"},{-1,"1.0"},{-1.2,"1.2"}},
	None,None}, AspectRatio -> 1];
Show[GraphicsArray[{plot1, plot2}], AspectRatio -> 1/2];
Display["junk_ma.eps", %, "EPS", ImageSize -> 432];


