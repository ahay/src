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
MGP = {v -> 1/Sqrt[C], h -> 1/Sqrt[A], 
    q -> 1/(((B + D)^2 + B(C - B))/(A(C - B)))}
Fow[a_, b_, x_] := (v^2 a^2 + h^2 b^2)(1 - x) + 
    x Sqrt[(v^2 a^2 + h^2 b^2)^2 + 2(q - 1) a^2 b^2 v^2 h^2/x];
Tariq[a_,b_]:= (b*Sqrt[(b^2*h^4*q^2*(b^6*h^6*q^3 - 3*a^2*b^4*h^4*(-3 + q)*
        q^2*v^2 + 3*a^4*b^2*h^2*q*(1 + 2*q)*v^4 + 
       4*a^6*v^6))/((b^2*h^2*q + a^2*v^2)*
      (b^6*h^6*q^4 + a^2*b^4*h^4*q^2*(1 + 5*q)*v^2 + 
       a^4*b^2*h^2*q*(-13 + 22*q)*v^4 + 4*a^6*v^6))] + 
  a*v*Sqrt[(a^2*b^6*h^6*q^3*(1 + 3*(-1 + q)*q)*v^2 - 
      6*a^4*b^4*h^4*q^2*(2 + (-4 + q)*q)*v^4 + 
      9*a^6*b^2*h^2*q*(-1 + 2*q)*v^6 + 4*a^8*v^8)/
     (b^8*h^8*q^4 + a^2*b^6*h^6*q^3*(10 + 3*(-2 + q)*q)*
       v^2 - 3*a^4*b^4*h^4*q^2*(3 + 2*(-5 + q)*q)*v^4 + 
      a^6*b^2*h^2*q*(-5 + 18*q)*v^6 + 4*a^8*v^8)])^2;
Tsvan[a_,b_] := a^2/C + (b^2*(-B + C)*(b^2*C*(-B + C) + 
    a^2*(D^2 + B*(C + 2*D))))/(A*b^2*(B - C)^2*C + 
   a^2*(D^2 + B*(C + 2*D))^2);
ParametricPlot[{ArcCos[Sqrt[AngP[Cos[a], Sin[a]]]] 180/Pi, 100
        (Sqrt[
              1/(Fow[Sqrt[AngP[Cos[a], Sin[a]]], 
                      Sqrt[1 - AngP[Cos[a], Sin[a]]],1/(2 (1+q))] 
	GruP[Cos[a], Sin[a]])] - 1)} /. MGP /. GS, {a, 0, Pi/2}]
ParametricPlot[{ArcCos[Sqrt[AngP[Cos[a], Sin[a]]]] 180/Pi, 100
        (Sqrt[
              1/(Tariq[Sqrt[AngP[Cos[a], Sin[a]]], 
                      Sqrt[1 - AngP[Cos[a], Sin[a]]]] 
	GruP[Cos[a], Sin[a]])] - 1)} /. MGP /. GS, {a, 0, Pi/2},
		PlotStyle->AbsoluteDashing[{2}]];
Show[{%,%%},PlotRange->All,
Frame->True,FrameLabel->{"Group Angle (degrees)",
"Relative Error (%)",None,None},PlotLabel->"Group Velocity Error"];
Display["junk_ma.ps",%,"EPS"];
