GS = {A -> 3.41, B -> 0.54, C -> 2.27, D -> 1.07};
GS = {A -> 14.47, C -> 9.57, B -> 2.28, D -> 4.51}; 
PhaSV[a_, b_] := 
  1/2((A + B) b^2 + (C + B) a^2 - 
        Sqrt[((A - B) b^2 - (C - B) a^2)^2 + 4 (D + B)^2 a^2 b^2]);
Muir2[a_, b_] := (a^4 v^4 + 
        (1 + (q1 a^2 v^2 + q2 b^2 h^2)/(a^2 v^2 + b^2 h^2))a^2 b^2 v^2 h^2 + 
                  b^4 h^4)/(a^2 v^2 + b^2 h^2);
MPSV = {v -> Sqrt[B], h -> Sqrt[B], 
        q1 -> (A(C - B) - (B + D)^2)/(B(C - B)), 
        q2 -> (C(A - B) - (B + D)^2)/(B(A - B))};
WeSV[a_, b_] := vs^2 + 2 (e - d) a^2 b^2 vp^2
Th = {vp -> Sqrt[C], vs -> Sqrt[B], e -> (A - C)/(2 C), 
    d -> ((B + D)^2 - (C - B)^2)/(2 C (C - B))};
Fow2[a_, b_, x_] := B ((1 - x) + 
    x Sqrt[1 + 2 (q1 a^2 + q2 b^2 - 1) a^2 b^2 /x]);
New[a_,b_] := (B*(3 + q1) + Sqrt[B^2*(-1 + q1)^2 + 4*(B+D)^2] - 
  2*Sqrt[(B+D)^2 + ((1 - 2*a2)^2*B*(-1 + q1)*(B*(-1 + q1) + 
        Sqrt[B^2*(-1 + q1)^2 + 4*(B+D)^2]))/2])/4;
New[a_,b_] := (A*b^2 + 2*B + (a^2*(A - B)*(D^2 + B*(C + 2*D)))/
   (B + D)^2 - 
  Sqrt[((-4*a^2*b^2*(A - B)*(B + D)^2*(D^2 + B*(C + 2*D))*
       (A*(B - C) + D^2 + B*(C + 2*D)))/(B - C) + 
     (a^2*B*(D^2 + B*(C + 2*D)) - 
       A*(b^2*(B + D)^2 + a^2*(B*C + 2*B*D + D^2)))^2)/
    (B + D)^4])/2;
Stopin[a_,b_] := (A*b^2 + 2*B + a^2*C - Sqrt[(A*b^2 + a^2*C)^2 + 
    4*a^2*A*b^2*C*(-1 + (B*C + 2*B*D + D^2)/
       (-(A*B) + A*C))])/2;
Stopin2[a_,b_]:= (a^2*B*C + A*b^2*(B - a^2*C*(-1 + (B*C + 2*B*D + D^2)/
       (-(A*B) + A*C))))/(A*b^2 + a^2*C);
ParametricPlot[{a 180/Pi, 
        Abs[Sqrt[WeSV[Cos[a], Sin[a]]/PhaSV[Cos[a], Sin[a]]] - 1]} /. Th /. 
    GS, {a, 0, Pi/2}, PlotStyle -> AbsoluteDashing[{2}]];
ParametricPlot[{a 180/Pi, 
          Abs[Sqrt[Muir2[Cos[a], Sin[a]]]/Sqrt[PhaSV[Cos[a], Sin[a]]] - 
              1]} /. MPSV /. GS, {a, 0, Pi/2}];
ParametricPlot[{a 180/Pi, 
          Abs[Sqrt[Stopin[Cos[a], Sin[a]]]/Sqrt[PhaSV[Cos[a], Sin[a]]] - 
              1]} /. MPSV /. GS, {a, 0, Pi/2},
	      PlotStyle -> AbsoluteDashing[{6}]];
ParametricPlot[{a 180/Pi, 
          Abs[Sqrt[New[Cos[a], Sin[a]]]/Sqrt[PhaSV[Cos[a], Sin[a]]] - 
              1]} /. MPSV /. GS, {a, 0, Pi/2},
	      PlotStyle -> AbsoluteDashing[{4}]];
Show[{%,%%,%%%,%%%%},PlotRange->All];
Display["junk_ma.ps",%,"EPS"];
