pcrs[r_, hmax_, tmax_] := 
 ParametricPlot3D[{Sqrt[((1 + d - n) (1 + (1 + d) n + 2 r + 
      r (r - Sqrt[n^2 + (1 + r)^2])))/n], d, 
   100 Abs[1 - 
      Sqrt[-(n Sqrt[
        n^2 + (1 + r)^2] ((
          d^2 (1 + r)^2 (-r + Sqrt[1 + (1 + r)^2]))/(1 + (1 + r)^2)^(
          3/2) + (r - (1 + d + (1 + r)^2)/Sqrt[
            1 + (1 + r)^2])^2 + ((1 + d - n) (1 + 
             r)^2 (1 + (1 + d) n + 2 r + 
             r (r - Sqrt[n^2 + (1 + r)^2])))/(
          n (1 + (1 + r)^2))))/((n r - (1 + d) Sqrt[
           n^2 + (1 + r)^2]) (1 + (1 + d) n + 2 r + 
          r (r - Sqrt[n^2 + (1 + r)^2])))]]}, {d, 0, 1}, {n, 0.01, 
   10}, PlotRange -> {{0, hmax}, {0, 1}, {0, tmax}}, 
  PlotLabel -> "(a) CRS",
  AxesLabel -> {"h/D", "d/D", "%"}, Ticks->{{0,1},{0,1},Automatic},
  Mesh -> Full, NormalsFunction -> None,
  BoxRatios -> {hmax, 1, 0.4}, PlotPoints -> {50, 50}];
pnh[r_, hmax_, tmax_] := 
 ParametricPlot3D[{Sqrt[((1 + d - n) (1 + (1 + d) n + 2 r + 
      r (r - Sqrt[n^2 + (1 + r)^2])))/n], d, 
   100 Abs[1 - 
      1/Sqrt[2] (\[Sqrt](-(n Sqrt[
              n^2 + (1 + r)^2] ((
                d^2 (1 + r)^2 (-r + Sqrt[
                   1 + (1 + r)^2]))/(1 + (1 + r)^2)^(
                3/2) + (-r + d/Sqrt[1 + (1 + r)^2] + Sqrt[
                  1 + (1 + r)^2])^2 + ((1 + d - n) (1/(
                   1 + (1 + r)^2) + (2 (1 + r)^2)/(
                   1 + (1 + r)^2) - ((1 + r)^2 (-r + Sqrt[
                    1 + (1 + r)^2]))/(1 + (1 + r)^2)^(
                   3/2)) (1 + (1 + d) n + 2 r + 
                   r (r - Sqrt[n^2 + (1 + r)^2])))/
                n + \[Sqrt]((((1 + r)^2 (-r + Sqrt[
                    1 + (1 + r)^2]) (d - 
                    Sqrt[((1 + d - n) (1 + (1 + d) n + 2 r + 
                    r (r - Sqrt[n^2 + (1 + r)^2])))/
                    n])^2)/(1 + (1 + r)^2)^(
                    3/2) + (-r + Sqrt[1 + (1 + r)^2] + (
                    d - Sqrt[((1 + d - n) (1 + (1 + d) n + 2 r + 
                    r (r - Sqrt[n^2 + (1 + r)^2])))/n])/Sqrt[
                    1 + (1 + r)^2])^2) (((1 + r)^2 (-r + Sqrt[
                    1 + (1 + r)^2]) (d + 
                    Sqrt[((1 + d - n) (1 + (1 + d) n + 2 r + 
                    r (r - Sqrt[n^2 + (1 + r)^2])))/
                    n])^2)/(1 + (1 + r)^2)^(
                    3/2) + (-r + Sqrt[1 + (1 + r)^2] + (
                    d + Sqrt[((1 + d - n) (1 + (1 + d) n + 2 r + 
                    r (r - Sqrt[n^2 + (1 + r)^2])))/n])/Sqrt[
                    1 + (1 + r)^2])^2))))/((n r - (1 + d) Sqrt[
                 n^2 + (1 + r)^2]) (1 + (1 + d) n + 2 r + 
                r (r - Sqrt[n^2 + (1 + r)^2])))))]}, {d, 0, 1}, {n, 
   0.01, 10}, PlotRange -> {{0, hmax}, {0, 1}, {0, tmax}}, 
  PlotLabel -> "(b) Non-hyperbolic CRS",
  AxesLabel -> {"h/D", "d/D", "%"}, Ticks->{{0,1},{0,1},Automatic},
  Mesh -> Full, NormalsFunction -> None,
  BoxRatios -> {hmax, 1, 0.4}, PlotPoints -> {50, 50}];
crs=pcrs[1, 1, 3];
nh=pnh[1, 1, 3];
GraphicsGrid[{{crs},{nh}}];
Export["junk_ma.eps", %, ImageSize-> {200,Automatic}];
