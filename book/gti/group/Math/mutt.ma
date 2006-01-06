GS = {A -> 3.41, B -> 0.54, C -> 2.27, D -> 1.07};
GS = {A -> 14.47, C -> 9.57, B -> 2.28, D -> 4.51}; 
AngP[a_,b_]:=(a^2*(A*b^2*(B - C) + a^2*(B - C)^2 + 
    b^2*(B*(B + C) + 4*B*D + 2*D^2) + 
    (B + C)*Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
       4*a^2*b^2*(B + D)^2])^2)/
 (2*(A^4*b^6 + a^6*(B - C)^2*(B^2 + C^2) + 
   b^4*B^3*(b^2*B + Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
       4*a^2*b^2*(B + D)^2]) + 
   A^3*b^4*(-2*b^2*B + 2*a^2*(B - C) + 
     Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
       4*a^2*b^2*(B + D)^2]) + 
   A^2*b^2*(a^4*B^2 + 3*a^2*b^2*B^2 + 2*b^4*B^2 - 
     2*a^4*B*C + 2*a^2*b^2*B*C + a^4*C^2 + a^2*b^2*C^2 + 
     8*a^2*b^2*B*D + 4*a^2*b^2*D^2 + 
     (-(b^2*B) + a^2*(B - C))*
      Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2]) + 
   a^2*b^2*(b^2*(B^2*(B + C)*(3*B + C) + 4*B^2*(3*B + C)*
        D + 2*B*(7*B + C)*D^2 + 8*B*D^3 + 2*D^4) + 
     (2*B + C)*(B*(B + C) + 4*B*D + 2*D^2)*
      Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2]) + 
   a^4*((B - C)^2*(B + C)*Sqrt[(b^2*(A - B) + a^2*(B - C))^
         2 + 4*a^2*b^2*(B + D)^2] + 
     b^2*(3*B^4 + 2*B^3*(C + 6*D) + 2*B*(C + 4*D)*
        (C^2 + D^2) + 2*D^2*(2*C^2 + D^2) + 
       B^2*(3*C^2 + 4*C*D + 14*D^2))) + 
   A*b^2*(Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2]*(-(b^2*B^2) + 
       a^2*(3*B^2 - C^2 + 4*B*D + 2*D^2)) + 
     2*(-(b^4*B^3) + a^4*(B - C)*(2*B^2 + C^2 + 2*B*D + 
         D^2) + a^2*b^2*(B*(B^2 - 2*B*C - C^2) + 
         2*B*(B - C)*D + (B - C)*D^2)))));
GruP[a_,b_]:=(b^2*(A^2*b^2 + a^2*A*B - 2*A*b^2*B + a^2*B^2 + b^2*B^2 - 
     a^2*A*C + a^2*B*C + 4*a^2*B*D + 2*a^2*D^2 + 
     (A + B)*Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2])^2 + 
  a^2*(A*b^2*(B - C) + a^2*(B - C)^2 + 
     b^2*(B*(B + C) + 4*B*D + 2*D^2) + 
     (B + C)*Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
        4*a^2*b^2*(B + D)^2])^2)/
 (2*(A^2*b^4 + b^4*B^2 + a^4*(B - C)^2 - 
   2*A*b^2*(b^2*B + a^2*(-B + C)) + 
   2*a^2*b^2*(B*(B + C) + 4*B*D + 2*D^2))*
  (b^2*(A + B) + a^2*(B + C) + 
   Sqrt[(b^2*(A - B) + a^2*(B - C))^2 + 
     4*a^2*b^2*(B + D)^2]));
Muir[a_, b_] := (a^4 v^4 + (1 + q)a^2 b^2 v^2 h^2 + b^4 h^4)/(a^2 v^2 + 
        b^2 h^2);
MGP = {v -> 1/Sqrt[C], h -> 1/Sqrt[A], 
    q -> 1/(((B + D)^2 + B(C - B))/(A(C - B)))};
AlT[a_, b_] := (v^2 a^2 + h^2 b^2) + 
	(q-1) a^2 b^2 v^2 h^2/(v^2 a^2 + q^2 h^2 b^2)
ParametricPlot[{ArcCos[Sqrt[AngP[Cos[a], Sin[a]]]] 180/Pi, 
        Abs[Sqrt[
              1/(Muir[Sqrt[AngP[Cos[a], Sin[a]]], 
                      Sqrt[1 - AngP[Cos[a], Sin[a]]] ] GruP[Cos[a], 
                      Sin[a]])] - 1]} /. MGP /. GS, {a, 0, Pi/2}];
ParametricPlot[{ArcCos[Sqrt[AngP[Cos[a], Sin[a]]]] 180/Pi, 
        Abs[Sqrt[
              1/(AlT[Sqrt[AngP[Cos[a], Sin[a]]], 
                      Sqrt[1 - AngP[Cos[a], Sin[a]]]] GruP[Cos[a], 
                      Sin[a]])] - 1]} /. MGP /. GS, {a, 0, Pi/2},
		PlotStyle->AbsoluteDashing[{4}]];
Show[{%,%%}];
Display["junk_ma.ps",%,"EPS"];
