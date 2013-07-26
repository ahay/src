t[x_] := Sqrt[t0^2 + x^2/w - ((-1 + g^2)^2 x^4 (1 + 2 Sqrt[1 - ((-1 + g^2)^2 x^2)/(16 g^2 H^2)]))/
                               (144 g^3 H^2 v^2 (1 +   Sqrt[1 - ((-1 + g^2)^2 x^2)/(16 g^2 H^2)])^2)];
tf[x_] := t0^2 + x^2/w + A x^4/(w^2 (t0^2 + B x^2/w + Sqrt[t0^4 + 2 B t0^2 x^2/w + c x^4/w^2]));
x[g_,t_] := 4 t/Sqrt[-1 + g^2];
Abs[(Sqrt[tf[x]] - t[x])/t[x]];
% /. {B -> (-1 + g + g^3 - g^4)/(18 g^2),
      c -> -((-1 + g)^4 (1 + g + g^2)^2)/(81 g^4),
      t0 -> (4 (1 + g + g^2) H)/(3 (g + g^2) v), 
      w -> 3 g^2 v^2/(1 + g + g^2), A -> -(-1 + g)^2/(6 g)} /. {H -> 1, v -> 1};
err = 100 % /. {x -> x[g,t]};
pp = ParametricPlot3D[{x[g,t], g, err}, {g, 1.001, 4}, {t, 0, 1}, 
 BoxRatios -> {1, 1, 0.4}, PlotLabel -> "(d) Generalized", 
 AxesLabel -> {"x/H", "v(H)/v0", "%      "}, 
 PlotRange -> {{0, 4}, All, All}];
Export["junk_ma.eps", Rasterize[pp, ImageResolution->200], "EPS"];
