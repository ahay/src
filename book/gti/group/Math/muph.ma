GS = {A -> 3.41, B -> 0.54, C -> 2.27, D -> 1.07};
GS = {A -> 14.47, C -> 9.57, B -> 2.28, D -> 4.51}; 
PhaSV[a_, b_] := 
  1/2((A + B) b^2 + (C + B) a^2 - 
        Sqrt[((A - B) b^2 - (C - B) a^2)^2 + 4 (D + B)^2 a^2 b^2]);
PhaP[a_, b_] := 
  1/2((A + B) b^2 + (C + B) a^2 + 
        Sqrt[((A - B) b^2 - (C - B) a^2)^2 + 4 (D + B)^2 a^2 b^2]);
Muir[a_, b_] := (a^4 v^4 + (1 + q)a^2 b^2 v^2 h^2 + b^4 h^4)/(a^2 v^2 + 
        b^2 h^2);
Muir2[a_, b_] := (a^4 v^4 + 
        (1 + (q1 a^2 v^2 + q2 b^2 h^2)/(a^2 v^2 + b^2 h^2))a^2 b^2 v^2 h^2 + 
                  b^4 h^4)/(a^2 v^2 + b^2 h^2);
MPP = {v -> Sqrt[C], h -> Sqrt[A], q -> ((B + D)^2 + B(C - B))/(A(C - B))};
MPSV = {v -> Sqrt[B], h -> Sqrt[B], 
        q1 -> (A(C - B) - (B + D)^2)/(B(C - B)), 
        q2 -> (C(A - B) - (B + D)^2)/(B(A - B))};
ParametricPlot[{Sqrt[PhaSV[Cos[a], Sin[a]]] Sin[a], 
      Sqrt[PhaSV[Cos[a], Sin[a]]]Cos[a]} /. GS, {a, 0, 2 Pi}, 
  AspectRatio -> Automatic, Frame -> True];
ParametricPlot[{Sqrt[PhaP[Cos[a], Sin[a]]] Sin[a], 
      Sqrt[PhaP[Cos[a], Sin[a]]]Cos[a]} /. GS, {a, 0, 2 Pi}, 
  AspectRatio -> Automatic, Frame -> True];
ParametricPlot[{Sqrt[Muir[Cos[a], Sin[a]]] Sin[a], 
        Sqrt[Muir[Cos[a], Sin[a]]] Cos[a]} /. MPP /. GS, {a, 0, 2Pi}, 
  PlotStyle -> Dashing[{0.01}], Frame -> True, AspectRatio -> Automatic];
ParametricPlot[{Sqrt[Muir2[Cos[a], Sin[a]]] Sin[a], 
        Sqrt[Muir2[Cos[a], Sin[a]]] Cos[a]} /. MPSV /. GS, {a, 0, 2Pi}, 
  PlotStyle -> Dashing[{0.01}], Frame -> True, AspectRatio -> Automatic];
Show[{%,%%,%%%,%%%%},Axes->False];
Display["junk_ma.ps",%,"EPS"];
