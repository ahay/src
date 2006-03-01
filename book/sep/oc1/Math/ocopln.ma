$DefaultFont={"Times",12};
ymax[x_, h_]:= 1/2 (x + Sqrt[4 h^2 + x^2]);		
f2[h_] := Plot[-Sqrt[y^2-h^2], {y, h, ymax[5.1, h]},
        PlotStyle -> {Thickness[0.01]},  
        DisplayFunction -> Identity];	
fronts2  = Show[f2[0], f2[1], f2[2], f2[3], f2[4], f2[5]];		
r2[x_] := Plot[-Sqrt [y x], {y, 0, ymax[x, 5]}, 
        PlotStyle -> {Dashing[{0.02}]}, 
        DefaultFont -> {"Times", 12}];
plot2  =  Show[r2[1], r2[2], r2[3], r2[4], r2[5], fronts2, 
	Frame->True, FrameLabel -> {"midpoint","t_n"},
	PlotLabel->"After NMO", RotateLabel->True,
        AspectRatio -> 1, PlotRange->{-12,0},FrameTicks->{Automatic,
	{{-2,"2"},{-4,"4"},{-6,"6"},{-8,"8"},{-10,"10"},{-12,"12"}},
	None,None}];
f1[h_] := Plot[-Sqrt[x^2 + 3 h^2], {x, ymax[0, h], ymax[5.1, h]}, 
        PlotStyle -> {Thickness[0.01]}, 
	DisplayFunction -> Identity]; 
fronts1 = Show[f1[0], f1[1], f1[2], f1[3], f1[4], f1[5]];
r1[y_] := Plot[-Sqrt[4 x^2 - 3 x y], {x, 3/4 y, ymax[y, 5]},
        PlotStyle -> {Dashing[{0.02}]}, 
        DefaultFont -> {"Times", 12}];
plot1  =  Show[r1[0], r1[1], r1[2], r1[3], r1[4], r1[5], fronts1, 
	Frame->True, FrameLabel -> {"midpoint", "t"},
	PlotLabel->"Before NMO", 
        AspectRatio -> 1, PlotRange->{-12,0},FrameTicks->{Automatic,
	{{-2,"2"},{-4,"4"},{-6,"6"},{-8,"8"},{-10,"10"},{-12,"12"}},
	None,None}];
Show[GraphicsArray[{plot1, plot2}],AspectRatio->1/2];
Display["junk_ma.eps", %, "EPS", ImageSize -> 432];

