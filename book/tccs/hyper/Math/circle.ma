X[n_, m_, r_, h_] := (
 4 (m - n) (h^2 + m n + 2 h r + r (r - Sqrt[n^2 + (h + r)^2])))/n; 
T[n_, m_, r_, h_] := -(
  4 (n r - m Sqrt[n^2 + (h + r)^2]) (h^2 + m n + 2 h r + 
     r (r - Sqrt[n^2 + (h + r)^2])))/(n Sqrt[n^2 + (h + r)^2]);
Hyp[n_, m_, r_, h_] := 
 4 ((r - Sqrt[
      m^2 + (h + r)^2])^2 + ((m - n) (h + r)^2 (h^2 + m n + 2 h r + 
       r (r - Sqrt[n^2 + (h + r)^2])))/(n (m^2 + (h + r)^2)));
Alk[n_, m_, r_, h_] := 
 4 ((r - Sqrt[
      m^2 + (h + r)^2])^2 + ((m - n) (h + r)^2 (h^2 + m n + 2 h r + 
       r^2 - r Sqrt[n^2 + (h + r)^2]))/(
    n (m^2 + (h + r)^2)) + (m^2 (m - n)^2 (h + r)^2 (-r + Sqrt[
         m^2 + (h + r)^2]) (h^2 + m n + 2 h r + r^2 - 
         r Sqrt[n^2 + (h + r)^2])^2)/(n^2 (m^2 + (h + r)^2)^(
       5/2) ((r - Sqrt[m^2 + (h + r)^2])^2 + 
         1/(n (m^2 + (h + r)^2)) (m - n) (h + r)^2 (h^2 + m n + 
             2 h r + r^2 - r Sqrt[n^2 + (h + r)^2]) (1 + (
             m^2 (-1 + r/Sqrt[m^2 + (h + r)^2]))/(h + r)^2))));
Mal[n_, m_, r_, 
  h_] := (2 (-r + Sqrt[m^2 + (h + r)^2]) (1 - 1/(
      1 + (4 m^2 (-1 + r/Sqrt[m^2 + (h + r)^2]))/(h + r)^2)) + (
   2 Sqrt[(r - Sqrt[
       m^2 + (h + r)^2])^2 + ((m - n) (h + r)^2 (h^2 + m n + 2 h r + 
        r^2 - r Sqrt[n^2 + (h + r)^2]) (1 + (
        4 m^2 (-1 + r/Sqrt[m^2 + (h + r)^2]))/(h + r)^2))/(
     n (m^2 + (h + r)^2))])/(
   1 + (4 m^2 (-1 + r/Sqrt[m^2 + (h + r)^2]))/(h + r)^2))^2;
Gen[n_, m_, r_, h_] := 
 4 (r - Sqrt[m^2 + (h + r)^2])^2 + (
  4 (m - n) (h + r)^2 (h^2 + m n + 2 h r + 
     r (r - Sqrt[n^2 + (h + r)^2])))/(
  n (m^2 + (h + r)^2)) + (32 m^2 (m - n)^2 (h + r)^2 (-r + Sqrt[
       m^2 + (h + r)^2]) (h^2 + m n + 2 h r + 
       r (r - Sqrt[n^2 + (h + r)^2]))^2)/(n^2 (m^2 + (h + r)^2)^(
     5/2) (4 (r - Sqrt[m^2 + (h + r)^2])^2 + (
       4 (m - n) (h + r)^2 (h^2 + m n + 2 h r + 
          r (r - Sqrt[n^2 + (h + r)^2])) (2 - (2 r)/Sqrt[
          m^2 + (h + r)^2] - (
          m^2 (r - Sqrt[m^2 + (h + r)^2])^2)/((h + r)^2 (m^2 + 
             2 r (h + r - Sqrt[m^2 + (h + r)^2])))))/(
       n (m^2 + (h + 
            r)^2)) + \[Sqrt](16 (r - Sqrt[m^2 + (h + r)^2])^4 + (
          16 m^4 (m - n)^2 (r - Sqrt[m^2 + (h + r)^2])^4 (h^2 + m n + 
             2 h r + r (r - Sqrt[n^2 + (h + r)^2]))^2)/(
          n^2 (m^2 + (h + r)^2)^2 (m^2 + 
             2 r (h + r - Sqrt[m^2 + (h + r)^2]))^2) + (
          32 (m - n) (h + r)^2 (r - Sqrt[m^2 + (h + r)^2])^2 (h^2 + 
             m n + 2 h r + r (r - Sqrt[n^2 + (h + r)^2])) (2 - (2 r)/
             Sqrt[m^2 + (h + r)^2] - (
             m^2 (r - Sqrt[m^2 + (h + r)^2])^2)/((h + r)^2 (m^2 + 
                2 r (h + r - Sqrt[m^2 + (h + r)^2])))))/(
          n (m^2 + (h + r)^2)))));
HypErr[n_, m_, r_] := Abs[Sqrt[Hyp[n, m, r, 1]/T[n, m, r, 1]] - 1];
AlkErr[n_, m_, r_] := Abs[Sqrt[Alk[n, m, r, 1]/T[n, m, r, 1]] - 1];
MalErr[n_, m_, r_] := Abs[Sqrt[Mal[n, m, r, 1]/T[n, m, r, 1]] - 1];
GenErr[n_, m_, r_] := Abs[Sqrt[Gen[n, m, r, 1]/T[n, m, r, 1]] - 1];
h = ParametricPlot3D[{Sqrt[X[n, 1, r, 1]], r, 100 HypErr[n, 1, r]}, 
 {r, 0, 4}, {n, 0.1, 1}, PlotRange -> {{0, 4}, All, {0,12}}, 
 BoxRatios -> {1, 1, 0.4}, PlotLabel -> "(a) Hyperbolic", 
 AxesLabel -> {"x/H", "R/H", "%"}];
m = ParametricPlot3D[{Sqrt[X[n, 1, r, 1]], r, 100 MalErr[n, 1, r]}, 
 {r, 0, 4}, {n, 0.1, 1}, PlotRange -> {{0, 4}, All, {0,40}}, 
 BoxRatios -> {1, 1, 0.4}, PlotLabel -> "(b) Shifted hyperbola", 
 AxesLabel -> {"x/H", "R/H", "%"}]; 
a = ParametricPlot3D[{Sqrt[X[n, 1, r, 1]], r, 100 AlkErr[n, 1, r]}, 
 {r, 0, 4}, {n, 0.1, 1}, PlotRange -> {{0, 4}, All, {0,7}}, 
 BoxRatios -> {1, 1, 0.4}, PlotLabel -> "(c) Alkhalifah-Tsvankin", 
 AxesLabel -> {"x/H", "R/H", "%"}]; 
g = ParametricPlot3D[{Sqrt[X[n, 1, r, 1]], r, 100 GenErr[n, 1, r]}, 
 {r, 0, 4}, {n, 0.1, 1}, PlotRange -> {{0, 4}, All, All}, 
 BoxRatios -> {1, 1, 0.4}, PlotLabel -> "(d) Generalized", 
 AxesLabel -> {"x/H", "R/H", "%"}]; 
ga = GraphicsArray[{{Rasterize[h], Rasterize[m]}, 
                    {Rasterize[a], Rasterize[g]}}];
Export["junk_ma.eps", ga, "EPS"];
