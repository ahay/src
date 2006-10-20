r[x_,y_,g_]:=2(Cos[Sqrt[x^2+y^2]] - 1)-
	2(g (Cos[x] + Cos[y] - 2) + (1-g) (Cos[x] Cos[y] -1));
p[g_,a_]:=
  ContourPlot[Abs[r[x,y,g]],{x,0,Pi},{y,0,Pi},
    Contours->{0.2,0.4,0.6,0.8,1,1.2, 1.4, 1.5,1.6,1.8,2,2.2,2.4,2.6,2.8,3,
        3.2,3.4,3.5}, PlotPoints->30, ColorFunction->(GrayLevel[1-#]&),
    PlotLabel->FontForm[a,{"Helvetica",14}]];
Show[GraphicsArray[{p[1,"a"],p[0,"b"],p[0.5,"c"],p[0.65,"d"]}], 
  AspectRatio->1/4];
Display["junk_ma.eps",%,"EPS",ImageSize->864];
