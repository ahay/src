helix[z_, a_, c_] := 
 ParametricPlot3D[{(t + b)/c, Cos[t], Sin[t]}, {t, 2, z + 2.002}, {b, -a, a}, 
  Mesh -> {Table[x, {x, 2.001, z + 2.001, Pi/12}], {-a + 0.001, 
     a - 0.001}}, PlotPoints -> 101, MeshStyle -> Thick];
patch[z1_, z2_, a_, c_] := 
 ParametricPlot3D[{(t + b)/c, Cos[t], Sin[t]}, {t, 2*Pi*(z1 - 1)/25, 2*Pi*(z2 - 1)/25}, {b, -a, a}, 
  PlotPoints -> {10*(z2 - z1 + 1), 20}, Mesh -> 50];
h1=helix[8*Pi,Pi,2];
h2=helix [8*Pi,2,2];
p21=patch[35,36, 2,2];
p11 = patch[35,36, Pi,2];
p22 = patch[59,62,2,2];
p12 =patch[59,62,Pi,2];
p23 = patch[85,86, 2,2];
p13 = patch[85,86, Pi,2];
p3 = Show[h1, p11, p12, p13, BoxRatios -> {1, 1, 1}, Boxed -> False, 
  Axes -> False, PlotRange -> All, PlotLabel -> b];
p2 = Show[h2, p22, p21, p23, BoxRatios -> {1, 1, 1}, Boxed -> False, 
  Axes -> False, PlotRange -> All, PlotLabel -> c];
t = {{0, 0, 0, 0, 0}, 
     {0, 0, 0, 0, 0}, 
     {0, 0, 0, 0, 0}, 
     {0, 0, 0, 0, 0}, 
     {0, 0, 0, 0, 0}, 
     {0, 0, 0, 0, 0}, 
     {0, 0, 0, 0, 0}, 
     {0, 0, 1, 0, 0}, 
     {0, 1, 1, 1, 0}, 
     {0, 0, 1, 0, 0}, 
     {0, 0, 0, 0, 0}};
color[f_] := If[f == 0, GrayLevel[0.6], GrayLevel[0.2]];
p1 = ListDensityPlot[t, ColorFunction -> (color[#] &), 
  AspectRatio -> Automatic, Axes -> False, Frame -> False, 
  PlotLabel -> a, InterpolationOrder -> 0, Mesh -> All];
p4 = Show[
  Graphics[Raster[{{0.6, 0.6, 0.6, 0.6, 0.2, 0.6, 0.6, 0.6, 0.6, 0.6, 
      0.6, 0.6, 0.6, 0.6, 0.2, 0.2, 0.2, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
       0.6, 0.6, 0.6, 0.2, 0.6, 0.6, 0.6, 0.6}}]], 
  AspectRatio -> 1/31, PlotLabel -> d];
Show[Graphics[{Rectangle[{0, 5}, {40, 10}, p4], 
   Rectangle[{0, 10}, {25, 30}, p1], 
   Rectangle[{10, 10}, {30, 30}, p3], 
   Rectangle[{25, 10}, {45, 30}, p2]}]];
Display["junk_ma.eps",Rasterize[%,ImageResolution->600],"EPS"];

