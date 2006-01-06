GS = {A -> 3.41, B -> 0.54, C -> 2.27, D -> 1.07};
GS = {B -> 7.19312, C -> 15.2178, A -> 19.3885, D -> 0.646813};
PhaP[a_, b_] := 
  1/2((A + B) b^2 + (C + B) a^2 + 
        Sqrt[((A - B) b^2 - (C - B) a^2)^2 + 4 (D + B)^2 a^2 b^2]);
Muir[a_, b_] := (a^4 v^4 + (1 + q)a^2 b^2 v^2 h^2 + b^4 h^4)/(a^2 v^2 + 
        b^2 h^2);
MPP = {v -> Sqrt[C], h -> Sqrt[A], q -> ((B + D)^2 + B(C - B))/(A(C - B))};
WeP[a_, b_] := vp^2 (1 + 2 e b^4 + 2 d a^2 b^2);
Th = {vp -> Sqrt[C], vs -> Sqrt[B], e -> (A - C)/(2 C), 
    d -> ((B + D)^2 - (C - B)^2)/(2 C (C - B))};
Fow[a_, b_, x_] := (v^2 a^2 + h^2 b^2)(1 - x) + 
    x Sqrt[(v^2 a^2 + h^2 b^2)^2 + 2(q - 1) a^2 b^2 v^2 h^2/x];
ParametricPlot[{a 180/Pi, 
        Abs[Sqrt[WeP[Cos[a], Sin[a]]/PhaP[Cos[a], Sin[a]]] - 1]} /. Th /. 
    GS, {a, 0, Pi/2}, PlotStyle -> AbsoluteDashing[{2}]];
ParametricPlot[{a 180/Pi, 
          Abs[Sqrt[Muir[Cos[a], Sin[a]]]/Sqrt[PhaP[Cos[a], Sin[a]]] - 
              1]} /. MPP /. GS, {a, 0, Pi/2}];
ParametricPlot[{a 180/Pi, 
          Abs[Sqrt[Fow[Cos[a], Sin[a],1/2]]/Sqrt[PhaP[Cos[a], Sin[a]]] - 
              1]} /. MPP /. GS, {a, 0, Pi/2}, PlotStyle->AbsoluteDashing[{4}]];
Show[{%,%%,%%%}];
Display["junk_ma.ps",%,"EPS"];
