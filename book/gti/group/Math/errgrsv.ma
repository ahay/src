GS = {A -> 3.41, B -> 0.54, C -> 2.27, D -> 1.07};
GS = {A -> 14.47, C -> 9.57, B -> 2.28, D -> 4.51}; 
AngSV[a_,b_]:=-(a^2*(A*b^2*(B - C) + a^2*(B - C)^2 + 
     b^2*(B*(B + C) + 4*B*D + 2*D^2) - 
     (B + C)*Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2])^2)/
 (2*(-(A^4*b^6) - a^6*(B - C)^2*(B^2 + C^2) + 
   b^4*B^3*(-(b^2*B) + Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
       4*a^2*b^2*(B + D)^2]) + 
   A^3*b^4*(2*b^2*B - 2*a^2*(B - C) + 
     Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
       4*a^2*b^2*(B + D)^2]) + 
   a^2*b^2*(-(b^2*(B^2*(B + C)*(3*B + C) + 
        4*B^2*(3*B + C)*D + 2*B*(7*B + C)*D^2 + 8*B*D^3 + 
        2*D^4)) + (2*B + C)*(B*(B + C) + 4*B*D + 2*D^2)*
      Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2]) + 
   A^2*b^2*(-2*b^4*B^2 - a^4*(B - C)^2 - 
     (b^2*B + a^2*(-B + C))*
      Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2] - a^2*b^2*(3*B^2 + C^2 + 
       4*D^2 + 2*B*(C + 4*D))) + 
   a^4*((B - C)^2*(B + C)*Sqrt[(b^2*(A - B) + a^2*(B - C))^
         2 + 4*a^2*b^2*(B + D)^2] - 
     b^2*(3*B^4 + 2*B^3*(C + 6*D) + 2*B*(C + 4*D)*
        (C^2 + D^2) + 2*D^2*(2*C^2 + D^2) + 
       B^2*(3*C^2 + 4*C*D + 14*D^2))) + 
   A*b^2*(Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2]*(-(b^2*B^2) + 
       a^2*(3*B^2 - C^2 + 4*B*D + 2*D^2)) + 
     2*(b^4*B^3 - a^4*(B - C)*(2*B^2 + C^2 + 2*B*D + D^2) + 
       a^2*b^2*(-B^3 + 2*B^2*(C - D) + C*D^2 + 
         B*(C^2 + 2*C*D - D^2))))));
GruSV[a_,b_]:=(b^2*(A^2*b^2 + a^2*A*B - 2*A*b^2*B + a^2*B^2 + b^2*B^2 - 
     a^2*A*C + a^2*B*C + 4*a^2*B*D + 2*a^2*D^2 - 
     (A + B)*Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2])^2 + 
  a^2*(A*b^2*(B - C) + a^2*(B - C)^2 + 
     b^2*(B*(B + C) + 4*B*D + 2*D^2) - 
     (B + C)*Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2])^2)/
 (2*(A^2*b^4 + b^4*B^2 + a^4*(B - C)^2 - 
   2*A*b^2*(b^2*B + a^2*(-B + C)) + 
   2*a^2*b^2*(B*(B + C) + 4*B*D + 2*D^2))*
  (b^2*(A + B) + a^2*(B + C) - 
   Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
     4*a^2*b^2*(B + D)^2]));
Muir2[a_, b_] := (a^4 v^4 + 
        (1 + (q1 a^2 v^2 + q2 b^2 h^2)/(a^2 v^2 + b^2 h^2))a^2 b^2 v^2 h^2 + 
                  b^4 h^4)/(a^2 v^2 + b^2 h^2);
MGSV = {v -> 1/Sqrt[B], h -> 1/Sqrt[B], p->1/Sqrt[A],
        q1 -> 1/((A(C - B) - (B + D)^2)/(B(C - B))), 
        q2 -> 1/((C(A - B) - (B + D)^2)/(B(A - B)))};
WeSV[a_, b_] := vs^2 + 2 (e - d) a^2 b^2 vp^2;
Th = {vp -> Sqrt[C], vs -> Sqrt[B], e -> (A - C)/(2 C), 
    d -> ((B + D)^2 - (C - B)^2)/(2 C (C - B))};
Fow2[a_, b_, x_] := (v^2 a^2 + h^2 b^2)(1 - x) + 
    x Sqrt[(v^2 a^2 + h^2 b^2)^2 + 2(q1 a^2 + q2 b^2 - 1) a^2 b^2 v^2 h^2/x];
Stopin[a_,b_,x_]:= B^(-1) + (b^2/A + a^2/C)*x - 
 Sqrt[(b^2/A + a^2/C)^2 + 
    (2*a^2*b^2*(-1 + (A*(-B + C))/(D^2 + B*(C + 2*D))))/
     (A*C*x)]*x;
ParametricPlot[{ArcCos[Sqrt[AngSV[Cos[a], Sin[a]]]] 180/Pi, 
        Abs[Sqrt[
              WeSV[Sqrt[AngSV[Cos[a], Sin[a]]], 
                      Sqrt[1 - AngSV[Cos[a], Sin[a]]]]/ GruSV[Cos[a], 
                      Sin[a]]] - 1]} /. Th /. GS, {a, 0, Pi/2},
		PlotStyle->AbsoluteDashing[{2}]];
ParametricPlot[{ArcCos[Sqrt[AngSV[Cos[a], Sin[a]]]] 180/Pi, 
        Abs[Sqrt[
              1/(Muir2[Sqrt[AngSV[Cos[a], Sin[a]]], 
                      Sqrt[1 - AngSV[Cos[a], Sin[a]]] ] GruSV[Cos[a], 
                      Sin[a]])] - 1]} /. MGSV /. GS, {a, 0, Pi/2}];
ParametricPlot[{ArcCos[Sqrt[AngSV[Cos[a], Sin[a]]]] 180/Pi, 
        Abs[Sqrt[
              1/(Stopin[Sqrt[AngSV[Cos[a], Sin[a]]], 
                      Sqrt[1 - AngSV[Cos[a], Sin[a]]],
        (D^2 + B*(C + 2*D))/(2*(A*(-B + C) + D^2 + B*(C + 2*D)))] 
	GruSV[Cos[a], Sin[a]])] - 1]} /. MGSV /. GS, {a, 0, Pi/2},
		PlotStyle->AbsoluteDashing[{6}]];
Show[{%,%%,%%%}];
Display["junk_ma.ps",%,"EPS"];
