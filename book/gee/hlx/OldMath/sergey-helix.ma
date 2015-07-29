helix[z_,a_,c_]:= ParametricPlot3D[{(t+b)/c,Cos[t],Sin[t]},{t,2,z+2},{b,-a,a},
	PlotPoints->{101,2}];
patch[z1_,z2_,a_,c_]:=
  ParametricPlot3D[{(t+b)/c,Cos[t],Sin[t]},{t,2*Pi*(z1-1)/25,2*Pi*(z2-1)/25},{
      b,-a,a},PlotPoints->{10*(z2-z1+1),20}];
h1=helix[8*Pi,Pi,2];
h2=helix [8*Pi,2,2];
p21=patch[35,36, 2,2];
p11 = patch[35,36, Pi,2];
p22 = patch[59,62,2,2];
p12 =patch[59,62,Pi,2];
p23 = patch[85,86, 2,2];
p13 = patch[85,86, Pi,2];
p3 = Show[h1, p11, p12, p13, BoxRatios->{1,1,1},Boxed->False,Axes->False, 
    PlotRange-> All, PlotLabel->FontForm[b,{"Helvetica",14}]];
p2 = Show[h2, p22, p21, p23, BoxRatios->{1,1,1},Boxed->False,Axes->False, 
    PlotRange-> All, PlotLabel->FontForm[c,{"Helvetica",14}]];
t={{0,0,0,0},
   {0,0,0,0},
   {0,0,0,0},
   {0,0,0,0},
   {0,0,0,0},
   {0,0,0,0},
   {0,0,0,0},
   {0,0,1,0},
   {0,1,1,1},
   {0,0,1,0}};
color[f_]:=If[f==0,GrayLevel[0.6],GrayLevel[0.2]];
p1 = ListDensityPlot[t , ColorFunction->(color[#]&), AspectRatio->Automatic, 
    Axes->False, Frame->False, PlotLabel->FontForm[a,{"Helvetica",14}]];
p4=Show[Graphics[Raster[{{0.6,0.6,0.6,0.6,0.2,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.2,0.2,0.2,0.6,
    	0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.2,0.6,0.6,0.6,0.6}}]],
    	AspectRatio->1/31,PlotLabel->FontForm[d,{"Helvetica",14}]];
Show[Graphics[{Rectangle[{-6,2},{25,3},p4],Rectangle[{0,3},{5,8},p1],Rectangle[{4,-2},{11,12},p3], 
      Rectangle[{10,-2},{17,12},p2]}]];
Display["junk_ma.eps",%, "EPS", ImageSize->864];

