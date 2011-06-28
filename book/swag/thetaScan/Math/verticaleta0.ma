VnmTTI[angle_, delta1_, epsilon1_, v1_, vs1_, theta1_, delta2_, 
  epsilon2_, v2_, vs2_, theta2_, delta3_, epsilon3_, v3_, vs3_, 
  theta3_, lyz1_, lyz2_, lyz3_] := Block[{c, a, l, f}, cc1 = 1;
  aa1 = 2*epsilon1 + 1;
  ll1 = vs1*vs1/(v1 v1);
  ff1 = Sqrt[2*delta1*(1 - ll1) + (1 - ll1)^2] - ll1;
  vv1 = v1 v1;
  ovv1 = 1/vv1;
  cc2 = 1;
  aa2 = 2*epsilon2 + 1;
  ll2 = vs2*vs2/(v2 v2);
  ff2 = Sqrt[2*delta2*(1 - ll2) + (1 - ll2)^2] - ll2;
  vv2 = v2 v2;
  ovv2 = 1/vv2;
  cc3 = 1;
  aa3 = 2*epsilon3 + 1;
  ll3 = vs3*vs3/(v3 v3);
  ff3 = Sqrt[2*delta3*(1 - ll3) + (1 - ll3)^2] - ll3;
  vv3 = v3 v3;
  ovv3 = 1/vv3;
  c1 = 1 + ll1;
  c2 = aa1 + ll1;
  w1 = aa1 - ll1;
  w2 = 1 - ll1;
  w3 = (ff1 + ll1)*(ff1 + ll1);
  si = Sin[angle - theta1];
  ci = Cos[angle - theta1];
  co2 = ci*ci;
  si2 = si*si;
  sc = si*ci;
  sz2 = si*si;
  cz2 = ci*ci;
  alfa = c2*sz2 + c1*cz2;
  bbb = w1*sz2 - w2*cz2;
  beta = Sqrt[bbb*bbb + 4*sz2*cz2*w3];
  vcps = N[v1*Sqrt[.5*(alfa + beta)]];
  pz = ci/vcps;
  px = si/vcps;
  px2 = px px;
  pz2 = pz pz;
  pxz = pz px;
  pzl1 = pz;
  gamma11 = aa1*px2 + ll1*pz2;
  gamma33 = pz2 + ll1*px2;
  gamma13 = (ff1 + ll1)*pxz;
  den = 1/(gamma11 + gamma33 - 2*ovv1);
  g11 = (gamma33 - ovv1)*den;
  g33 = (gamma11 - ovv1)*den;
  g13 = -gamma13*den;
  dx = N[vv1*(aa1*px*g11 + (ff1 + ll1)*pz*g13 + ll1*g33*px)];
  dz = N[vv1*(pz*g33 + (ff1 + ll1)*px*g13 + ll1*g11*pz)];
  dxT = dx Cos[theta1] + dz Sin[theta1];
  dzT = dz Cos[theta1] - dx Sin[theta1];
  
  t1 = N[lyz1/dzT];
  x1 = N[t1*dxT];
  dz1 = dz;
  g131 = g13;
  g111 = g11;
  g331 = g33;
  f1 = vv2*(aa2 + ll2);
  f2 = vv2*(1 + ll2);
  f3 = vv2*(aa2 - ll2);
  f4 = vv2*(1 - ll2);
  f5 = 2.*(ff2 + ll2)*(ff2 + ll2)*vv2*vv2;
  f6 = 2.;
  eps = .0001;
  px2 = px*px;
  alpha = f2*f2 - f4*f4;
  beta1 = 2*((f1*f2 + f3*f4 - f5)*px2 - f2*f6);
  gamma = N[f6*f6 - (2.*f1*f6 - (f1*f1 - f3*f3)*px2)*px2];
  det = beta1*beta1 - 4.*alpha*gamma;
  rad = Sqrt[det];
  signbeta = Abs[beta1]/beta1;
  q = -.5*(beta1 + signbeta*rad);
  pz2 = N[gamma/q];
  vp = 1/Sqrt[px2 + pz2];
  pz = Sqrt[pz2];
  pxz = px pz;
  pzl2 = pz;
  gamma11 = aa2*px2 + ll2*pz2;
  gamma33 = pz2 + ll2*px2;
  gamma13 = (ff2 + ll2)*pxz;
  den = 1/(gamma11 + gamma33 - 2*ovv2);
  g11 = (gamma33 - ovv2)*den;
  g33 = (gamma11 - ovv2)*den;
  g13 = -gamma13*den;
  dx = N[vv2*(aa2*px*g11 + (ff2 + ll2)*pz*g13 + ll2*g33*px)];
  dz = N[vv2*(pz*g33 + (ff2 + ll2)*px*g13 + ll2*g11*pz)];
  t2 = N[lyz2/dz];
  x2 = N[t2*dx];
  dz2 = dz;
  f1 = vv3*(aa3 + ll3);
  f2 = vv3*(1 + ll3);
  f3 = vv3*(aa3 - ll3);
  f4 = vv3*(1 - ll3);
  f5 = 2.*(ff3 + ll3)*(ff3 + ll3)*vv3*vv3;
  f6 = 2.;
  eps = .0001;
  px2 = px*px;
  alpha = f2*f2 - f4*f4;
  beta1 = 2*((f1*f2 + f3*f4 - f5)*px2 - f2*f6);
  gamma = N[f6*f6 - (2.*f1*f6 - (f1*f1 - f3*f3)*px2)*px2];
  det = beta1*beta1 - 4.*alpha*gamma;
  rad = Sqrt[det];
  signbeta = Abs[beta1]/beta1;
  q = -.5*(beta1 + signbeta*rad);
  pz2 = N[gamma/q];
  vp = 1/Sqrt[px2 + pz2];
  pz = Sqrt[pz2];
  pxz = px pz;
  gamma11 = aa3*px2 + ll3*pz2;
  gamma33 = pz2 + ll3*px2;
  gamma13 = (ff3 + ll3)*pxz;
  den = 1/(gamma11 + gamma33 - 2*ovv3);
  g11 = (gamma33 - ovv3)*den;
  g33 = (gamma11 - ovv3)*den;
  g13 = -gamma13*den;
  dx = N[vv3*(aa3*px*g11 + (ff3 + ll3)*pz*g13 + ll3*g33*px)];
  dz = N[vv3*(pz*g33 + (ff3 + ll3)*px*g13 + ll3*g11*pz)];
  t3 = N[lyz3/dz];
  x3 = N[t3*dx];
  offset = Abs[2*(x1 + x2 + x3)];
  offset = 2*(x1 + x2 + x3);
  time = N[2*(t1 + t2 + t3)];
  (*vnmo=1/Sqrt[((t1+t2)^2-t^2)/(src-rec)^2];*)
  Return[N[{offset, time}]];];
LinVeikdsEll[X_, delta1_, epsilon1_, v1_, theta1_, lyz1_] := 
 Block[{b, tt},
  z = lyz1;
  Z = 2 z;
  s = Sin[theta1];
  x = 0.5  X;
  x = X;
  eta = (epsilon1 - delta1)/(1 + 2 delta1);
  vv = v1;
  v = v1 Sqrt[1 + 2 delta1];
  z = Z;
  A = -4 eta; B = (1 + 8 eta + 8 eta^2)/(1 + 2 eta);
  CC = 1/(1 + 2 eta)^2;
  ttt = (Sqrt[
    x^2/v^2 + z^2/
     vv^2] (s vv^2 x^4 + 2 vv^2 x^3 z + 2 s (-v^2 + vv^2) x^2 z^2 + 
      2 v^2 x z^3 - s v^2 z^4))/(
   s vv^2 x^4 + 2 vv^2 x^3 z + 2 v^2 x z^3 - s v^2 z^4);
  ttt2 = Sqrt[
     x^2/v^2 + z^2/vv^2] + ((-v^2 + vv^2) x z Sqrt[
      x^2/v^2 + z^2/vv^2])/(vv^2 x^2 + v^2 z^2) s + (
     Sqrt[x^2/v^2 + z^2/
       vv^2] (-vv^4 x^4 - v^4 z^4 + v^2 vv^2 (x^4 + z^4)))/(
     2 (vv^2 x^2 + v^2 z^2)^2) s^2;
  ttt3 = (2 (vv^2 x^2 + v^2 z^2)^2 Sqrt[
      x^2/v^2 + z^2/
       vv^2])/(2 s (v - vv) (v + vv) x z (vv^2 x^2 + v^2 z^2) + 
      2 (vv^2 x^2 + v^2 z^2)^2 + 
      s^2 (v - vv) (v + vv) (-vv^2 x^4 + 2 (v - vv) (v + vv) x^2 z^2 +
          v^2 z^4));
  ttt4 = (2 (vv^2 x^2 + v^2 z^2) Sqrt[
      x^2/v^2 + z^2/
       vv^2] (vv^4 x^3 (x^3 + 2 s x^2 z + 2 x z^2 + 2 s z^3) - 
        v^4 z^3 (2 x^2 z + z^3 - 2 s x (x^2 + z^2)) - 
        v^2 vv^2 x z (x (x - z) z (x + z) + 
           2 s (x^2 + z^2)^2)))/(2 s (v - vv) (v + vv) x z (-vv x^2 + 
         v z^2) (vv x^2 + v z^2) (vv^2 x^2 + v^2 z^2) - 
      s^2 (v - vv) (v + vv) (vv^2 x^4 - v^2 z^4)^2 + 
      2 (vv^2 x^2 + v^2 z^2)^2 (-v^2 z^2 (2 x^2 + z^2) + 
         vv^2 (x^4 + 2 x^2 z^2)));
  ttt5 = Sqrt[x^2/v^2 + z^2/vv^2]/(
   1 + (s (v - vv) (v + vv) x z)/(vv^2 x^2 + v^2 z^2));
  ttT = Sqrt[
    Z^2 (1 - 2 delta1 s^2 + 2 (delta1 - epsilon1) s^4)/vv^2 + 
     x^2 (1 - 2 delta1 + 2 epsilon1 s^2 - 
         14 (epsilon1 - delta1) s^2 (1 - s^2))/vv^2];
  ttNew = 
   Sqrt[ (((-2 + s^2) v^2 - s^2 vv^2)^2 z^2)/(
     4 v^4 vv^2) - (((2 + 3 s^2) v^2 - 3 s^2 vv^2) ((-2 + s^2) v^2 - 
         s^2 vv^2))/(4 v^6) x^2 - 
     0.5 (s^2 (v^2 - vv^2)^2 ((-2 + s^2) v^2 - 3 s^2 vv^2))/(
      2 v^8 z^2) x^4];
  ttSh = Sqrt[(0.125` z^2 (-16.` v^8 vv^2 x^2 + 
         s^2 v^6 (8.` v^4 - 48.` v^2 vv^2 + 40.` vv^4) x^2 - 
         16.` v^10 z^2 + 
         s^6 v^2 ((6.` v^8 - 8.` v^6 vv^2 - 12.` v^4 vv^4 + 
               24.` v^2 vv^6 - 10.` vv^8) x^2 + 
            v^2 (-16.` v^6 + 48.` v^4 vv^2 - 48.` v^2 vv^4 + 
               16.` vv^6) z^2) + 
         s^8 ((-1.` v^10 - 2.` v^8 vv^2 + 18.` v^6 vv^4 - 
               32.` v^4 vv^6 + 23.` v^2 vv^8 - 6.` vv^10) x^2 + 
            v^2 (3.` v^8 - 12.` v^6 vv^2 + 18.` v^4 vv^4 - 
               12.` v^2 vv^6 + 3.` vv^8) z^2) + 
         s^4 v^4 (28.` vv^6 x^2 + v^4 vv^2 (52.` x^2 - 48.` z^2) + 
            v^2 vv^4 (-68.` x^2 + 24.` z^2) + 
            v^6 (-12.` x^2 + 24.` z^2))))/(v^4 vv^2 (s^2 ((1.` - 
               0.5` s^2) v^6 + (-2.` + 2.5` s^2) v^4 vv^2 + (1.` - 
               3.5` s^2) v^2 vv^4 + 1.5` s^2 vv^6) x^2 + 
         v^2 ((-2.` - 2.` s^2 + 1.5` s^4) v^4 + 
            s^2 (2.` - 3.` s^2) v^2 vv^2 + 1.5` s^4 vv^4) z^2))];
  Return[{ttt, ttt2, ttt3, ttt4, ttt5, ttt6, ttT, ttNew, ttSh}]];
LinVeikdsTTI[X_, delta1_, epsilon1_, v1_, theta1_, lyz1_] := 
 Block[{b, tt},
  z = lyz1;
  Z = 2 z;
  s = Sin[theta1];
  x = 0.5  X;
  x = X;
  theta = ArcTan[X/Z];
  sg = Sin[theta1 - theta];
  eta = (epsilon1 - delta1)/(1 + 2 delta1);
  vv = v1;
  v = v1 Sqrt[1 + 2 delta1];
  z = Z;
  
  ttt = (Sqrt[
      x^2/v^2 + z^2/
       vv^2] (2 (vv^2 x^2 + v^2 z^2)^2 ((2 + eta) vv^4 x^4 + 
           4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 
           2 v^4 z^4) ((2 + 3 eta) vv^4 x^4 + 
           4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 2 v^4 z^4)^2 + 
        s^2 ((2 + 3 eta)^3 (v - vv) vv^14 (v + vv) x^16 - 
           12 vv^12 (-(4 + 3 eta (8 + 15 eta)) v^4 + (4 + 
                 3 eta (8 + 15 eta (1 + eta))) v^2 vv^2 + 
              3 eta^3 vv^4) x^14 z^2 + 
           v^2 vv^10 (6 (20 + 
                 3 eta (44 + 9 (11 - 8 eta) eta)) v^4 - (128 + 
                 9 eta (92 + 3 eta (68 + 145 eta))) v^2 vv^2 + (8 + 
                 9 eta (4 + (6 - 61 eta) eta)) vv^4) x^12 z^4 - 
           4 v^4 vv^8 (4 (-10 + 
                 9 eta (-7 + 3 eta (-5 + 8 eta))) v^4 + (52 + 
                 27 eta (12 + 5 eta (5 + 23 eta))) v^2 vv^2 + 
              
              3 (-4 + 
                 3 eta (-8 + 3 eta (-5 + 29 eta))) vv^4) x^10 z^6 + 
           6 v^6 vv^6 (2 (10 + 51 eta + 72 eta^2) v^4 - (40 + 
                 9 eta (26 + 7 eta (7 + 40 eta))) v^2 vv^2 + (20 - 
                 33 eta (-4 + eta (-9 + 40 eta))) vv^4) x^8 z^8 + 
           16 v^8 vv^4 ((3 + 9 eta) v^4 - (1 + 3 eta) (13 + 33 eta + 
                 36 eta^2) v^2 vv^2 + (10 + 
                 9 eta (7 + (15 - 52 eta) eta)) vv^4) x^6 z^10 + 
           4 v^10 (v - vv) vv^2 (v + vv) (2 v^2 - 
              3 (10 + 51 eta + 72 eta^2) vv^2) x^4 z^12 - 
           48 (1 + 3 eta) v^12 (v - vv) vv^2 (v + vv) x^2 z^14 + 
           8 v^14 (-v^2 + vv^2) z^16) - 
        2 s x z (vv^2 x^2 + v^2 z^2) ((2 + 3 eta) vv^4 x^4 + 
           4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 
           2 v^4 z^4) ((-4 + eta (-8 + 3 eta)) vv^10 x^8 + 
           4 v^10 z^8 + 
           8 (1 + 2 eta) v^6 vv^4 x^2 z^4 ((3 + 9 eta) x^2 - 2 z^2) - 
           4 v^8 vv^2 z^6 (-4 (1 + 3 eta) x^2 + z^2) + 
           v^2 vv^8 x^6 ((4 + 3 eta (8 + 15 eta)) x^2 + 
              8 (-2 + 3 (-2 + eta) eta) z^2) + 
           8 v^4 vv^6 x^4 z^2 ((2 + 3 eta (4 + 9 eta)) x^2 + 
              3 (-1 + eta (-3 + 2 eta)) z^2))))/(2 (vv^2 x^2 + 
        v^2 z^2)^2 ((2 + 3 eta) vv^4 x^4 + 
        4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 2 v^4 z^4)^3);
  ttt2 = (2 (vv^2 x^2 + v^2 z^2)^2 Sqrt[
      x^2/v^2 + z^2/
       vv^2] ((2 + eta) vv^4 x^4 + 4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 
        2 v^4 z^4)^3)/(2 (vv^2 x^2 + v^2 z^2)^2 ((2 + eta) vv^4 x^4 + 
         4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 
         2 v^4 z^4)^2 ((2 + 3 eta) vv^4 x^4 + 
         4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 2 v^4 z^4) + 
      s^2 ((2 + eta) (2 + 3 eta)^2 vv^14 (-v^2 + vv^2) x^16 + 
         2 vv^12 ((-16 + eta (-44 + 9 eta (16 + 71 eta))) v^4 + 
            2 (2 + 9 eta) (2 + 3 eta (-1 + 4 eta)) v^2 vv^2 + (8 + 
               eta (20 + eta (-10 + 9 eta))) vv^4) x^14 z^2 + 
         v^2 vv^10 (6 (-4 + 
               3 eta (12 + eta (133 + 460 eta))) v^4 + (-64 + 
               eta (-508 + 3 eta (-380 + 1347 eta))) v^2 vv^2 + (88 + 
               eta (292 + eta (-230 + 351 eta))) vv^4) x^12 z^4 + 
         4 v^4 vv^8 ((20 + 
               eta (296 + 3 eta (461 + 1080 eta))) v^4 + (-68 + 
               eta (-452 - 735 eta + 3501 eta^2)) v^2 vv^2 + 
            3 (16 + eta (52 + eta (-88 + 201 eta))) vv^4) x^10 z^6 + 
         2 v^6 vv^6 (2 (50 + 
               eta (449 + 432 eta (3 + 2 eta))) v^4 + (-200 + 
               3 eta (-362 + eta (-299 + 2856 eta))) v^2 vv^2 + (100 +
                eta (188 + 7 eta (-169 + 504 eta))) vv^4) x^8 z^8 + 
         8 v^8 vv^4 (3 (8 + eta (49 + 72 eta)) v^4 + 
            2 (-17 + eta (-55 + 3 eta (19 + 36 eta))) v^2 vv^2 + (10 +
                eta (-37 - 298 eta + 936 eta^2)) vv^4) x^6 z^10 + 
         4 v^10 (v - vv) vv^2 (v + vv) ((22 + 72 eta) v^2 + 
            3 (2 + eta (37 + 72 eta)) vv^2) x^4 z^12 + 
         16 v^12 (v - vv) (v + 
            vv) (v^2 + (2 + 9 eta) vv^2) x^2 z^14 + 
         8 v^14 (v - vv) (v + vv) z^16) + 
      2 s x z (vv^2 x^2 + v^2 z^2) ((2 + eta) vv^4 x^4 + 
         4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 
         2 v^4 z^4) ((-4 + eta (-8 + 3 eta)) vv^10 x^8 + 4 v^10 z^8 + 
         8 (1 + 2 eta) v^6 vv^4 x^2 z^4 ((3 + 9 eta) x^2 - 2 z^2) - 
         4 v^8 vv^2 z^6 (-4 (1 + 3 eta) x^2 + z^2) + 
         v^2 vv^8 x^6 ((4 + 3 eta (8 + 15 eta)) x^2 + 
            8 (-2 + 3 (-2 + eta) eta) z^2) + 
         8 v^4 vv^6 x^4 z^2 ((2 + 3 eta (4 + 9 eta)) x^2 + 
            3 (-1 + eta (-3 + 2 eta)) z^2)));
  ttt3 = (Sqrt[
      x^2/v^2 + z^2/vv^2] + ((-v^2 + vv^2) x z Sqrt[
       x^2/v^2 + z^2/vv^2])/(vv^2 x^2 + v^2 z^2) s + (-vv^4 x^4 Sqrt[
       x^2/v^2 + z^2/vv^2])/(vv^2 x^2 + v^2 z^2)^2 eta + (
      Sqrt[x^2/v^2 + z^2/
        vv^2] (-vv^4 x^4 - v^4 z^4 + v^2 vv^2 (x^4 + z^4)))/(
      2 (vv^2 x^2 + v^2 z^2)^2) s^2 + (
      3 vv^6 x^6 (vv^2 x^2 + 4 v^2 z^2) Sqrt[x^2/v^2 + z^2/vv^2])/(
      2 (vv^2 x^2 + v^2 z^2)^4) eta^2 + (-vv^4 x^3 z ((3 v^2 + vv^2) x^2 + 4 v^2 z^2) Sqrt[
       x^2/v^2 + z^2/vv^2])/(vv^2 x^2 + v^2 z^2)^3 eta s);
  ttt4 = Sqrt[
   x^2/v^2 + z^2/vv^2 + (A X^4)/(
    v^2 (Z^2 + B X^2 + Sqrt[Z^4 + 2 B Z^2  X^2 + CC X^4]))];
  A0 = Z^2 (1 - 2 delta1 s^2 + 2 (delta1 - epsilon1) s^4)/vv^2;
  A2 = (1 - 2 delta1 + 2 epsilon1 s^2 - 
      14 (epsilon1 - delta1) s^2 (1 - s^2))/vv^2;
  A4 = -2 eta (1 - s^2)^2/(z^2 vv^2);
  A = A4/(1/(v^2 (1 + 2 eta)) - 1/v^2 );
  ttT = Sqrt[
    Z^2 (1 - 2 delta1 s^2 + 2 (delta1 - epsilon1) s^4)/vv^2 + 
     x^2 (1 - 2 delta1 + 2 epsilon1 s^2 - 
         14 (epsilon1 - delta1) s^2 (1 - s^2))/vv^2 - 
     2 eta (1 - s^2)^2 x^4/(z^2 vv^2)];
  ttT = Sqrt[A0 + A2 x^2 + A4 x^4];
  ttTF = Sqrt[A0 + A2 x^2 + A4 x^4/(1 + A x^2)];
  ttT9 = Sqrt[ (z^2)/ (vv^2) + 
     x^2 (1 - 2 delta1 + 2 epsilon1 s^2 - 
         14 (epsilon1 - delta1) s^2 (1 - s^2))/vv^2 - 
     2 eta (1 - s^2)^2 x^4/(z^2 vv^2)];
  ttT2 = Sqrt[-(z^2 (-((4 - 
                8 (delta1 + (-7 delta1 + 6 epsilon1) s^2 + 
                   7 (delta1 - epsilon1) s^4)) x^2 + (((-2 + 
                   s^2) v^2 - s^2 vv^2)^2 z^2)/v^4)^2 + 
          1/v^8 ((-2 + s^2) v^2 - 
             s^2 vv^2)^2 (-8 eta (-1 + s^2)^2 v^4 x^4 + 
             4 (1 - 2 (delta1 + (-7 delta1 + 6 epsilon1) s^2 + 
                   7 (delta1 - epsilon1) s^4)) v^4 x^2 z^2 + ((-2 + 
                   s^2) v^2 - 
                s^2 vv^2)^2 z^4)))/(16 vv^2 x^2 (2 eta (-1 + 
            s^2)^2 x^2 + (1 - 
            2 (delta1 + (-7 delta1 + 6 epsilon1) s^2 + 
               7 (delta1 - epsilon1) s^4)) z^2))];
  ttTSena = 
   Sqrt[X^2 + Z^2] Sqrt[
      1 - 2 delta1 sg^2 + 2 (delta1 - epsilon1) sg^4]/vv;
  Return[{ttt, ttt2, ttt3, ttt4, ttT9, ttT2, ttTSena, ttT, ttTF}]];
LinVeikdsTTIS[X_, delta1_, epsilon1_, v1_, theta1_, lyz1_] := 
 Block[{b, tt},
  z = lyz1;
  Z = 2 z;
  s = Sin[theta1];
  x = 0.5  X;
  x = X;
  theta = ArcTan[X/Z];
  sg = Sin[theta1 - theta];
  eta = (epsilon1 - delta1)/(1 + 2 delta1);
  vv = v1;
  v = v1 Sqrt[1 + 2 delta1];
  z = Z;
  Return[Sqrt[
     X^2 + Z^2] Sqrt[1 - 2 delta1 sg^2 + 2 (delta1 - epsilon1) sg^4]/
     vv]];
fxst1 = {};
fxst2 = {};
fxst3 = {};
fxst4 = {};
fxst5 = {};
fxst6 = {};
fxst7 = {};
fxst8 = {};
fxst9 = {};
vs = 0;
angle = 5. Pi/180;
delta = 0.2;
dd = 2.0;
Do[xss = 0.5*{VnmTTI[N[0.25*i*Pi/180], delta, 0.2, 2.0, 
       vs, -angle, .2, 0.2, 2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 
       0.0, 0.0][[1]] - 
     VnmTTI[-N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[1]], 
    VnmTTI[N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[2]] + 
     VnmTTI[-N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[2]]};
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[2]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[2]];
 fxst4 = Append[fxst4, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[3]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[3]];
 fxst5 = Append[fxst5, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[1]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[1]];
 fxst8 = Append[fxst8, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[6]];
 fxst6 = Append[fxst6, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - LinVeikdsTTIS[xss[[1]], delta, 0.2, 2, angle, 2];
 fxst7 = Append[fxst7, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[7]];
 fxst9 = Append[fxst9, {xss[[1]]/2, 100. dt/xss[[2]]}], {i, 0, 60}];
bb1 = ListPlot[fxst4, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`]}];
bb2 = ListPlot[fxst5, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`], 
    Dashing[{0.004`, 0.02`, 0.03`, 0.02`}]}];
bb3 = ListPlot[fxst7, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.5`]}];
bb3 = ListPlot[fxst9, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.5`]}];
bb4 = ListPlot[fxst8, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`], 
    Dashing[{0.004`, 0.02`, 0.03`, 0.02`}]}];
a1 = Show[bb1, bb2, bb3, bb4, Frame -> True, GridLines -> Automatic, 
  LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
  FrameLabel -> {{Error["%"],}, {Offset[km],}}, 
  PlotRange -> {{0, 0.1}, {-.05, 0.005}}];
fxst1 = {};
fxst2 = {};
fxst3 = {};
fxst4 = {};
fxst5 = {};
fxst6 = {};
fxst7 = {};
fxst8 = {};
fxst9 = {};
vs = 0;
angle = 10. Pi/180;
delta = 0.2;
dd = 2.0;
Do[xss = 0.5*{VnmTTI[N[0.25*i*Pi/180], delta, 0.2, 2.0, 
       vs, -angle, .2, 0.2, 2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 
       0.0, 0.0][[1]] - 
     VnmTTI[-N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[1]], 
    VnmTTI[N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[2]] + 
     VnmTTI[-N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[2]]};
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[2]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[2]];
 fxst4 = Append[fxst4, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[3]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[3]];
 fxst5 = Append[fxst5, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[1]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[1]];
 fxst8 = Append[fxst8, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[6]];
 fxst6 = Append[fxst6, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - LinVeikdsTTIS[xss[[1]], delta, 0.2, 2, angle, 2];
 fxst7 = Append[fxst7, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[7]];
 fxst9 = Append[fxst9, {xss[[1]]/2, 100. dt/xss[[2]]}], {i, 0, 60}];
bb1 = ListPlot[fxst4, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`]}];
bb2 = ListPlot[fxst5, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`], 
    Dashing[{0.004`, 0.02`, 0.03`, 0.02`}]}];
bb3 = ListPlot[fxst7, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.5`]}];
bb3 = ListPlot[fxst9, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.5`]}];
bb4 = ListPlot[fxst8, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`], 
    Dashing[{0.004`, 0.02`, 0.03`, 0.02`}]}];
a2 = Show[bb1, bb2, bb3, bb4, Frame -> True, GridLines -> Automatic, 
  LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
  FrameLabel -> {{Error["%"],}, {Offset[km],}}, 
  PlotRange -> {{0, 0.1}, {-.2, 0.01}}];
fxst1 = {};
fxst2 = {};
fxst3 = {};
fxst4 = {};
fxst5 = {};
fxst6 = {};
fxst7 = {};
fxst8 = {};
fxst9 = {};
vs = 0;
angle = 20. Pi/180;
delta = 0.2;
dd = 2.0;
Do[xss = 0.5*{VnmTTI[N[0.25*i*Pi/180], delta, 0.2, 2.0, 
       vs, -angle, .2, 0.2, 2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 
       0.0, 0.0][[1]] - 
     VnmTTI[-N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[1]], 
    VnmTTI[N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[2]] + 
     VnmTTI[-N[0.25*i*Pi/180], delta, 0.2, 2.0, vs, -angle, .2, 0.2, 
       2.0, vs, .0, .2, 0.2, 2.0, vs, 0.0, dd, 0.0, 0.0][[2]]};
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[2]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[2]];
 fxst4 = Append[fxst4, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[3]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[3]];
 fxst5 = Append[fxst5, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   0.5*LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[1]] - 
   0.5*LinVeikdsTTI[-xss[[1]], delta, 0.2, 2, angle, 2][[1]];
 fxst8 = Append[fxst8, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[6]];
 fxst6 = Append[fxst6, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - LinVeikdsTTIS[xss[[1]], delta, 0.2, 2, angle, 2];
 fxst7 = Append[fxst7, {xss[[1]]/2, 100. dt/xss[[2]]}];
 dt = xss[[2]] - 
   LinVeikdsTTI[xss[[1]], delta, 0.2, 2, angle, 2][[7]];
 fxst9 = Append[fxst9, {xss[[1]]/2, 100. dt/xss[[2]]}], {i, 0, 60}];
bb1 = ListPlot[fxst4, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`]}];
bb2 = ListPlot[fxst5, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`], 
    Dashing[{0.004`, 0.02`, 0.03`, 0.02`}]}];
bb3 = ListPlot[fxst7, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.5`]}];
bb3 = ListPlot[fxst9, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.5`]}];
bb4 = ListPlot[fxst8, Joined -> True, 
  PlotStyle -> {Thickness[0.006`], GrayLevel[0.0`], 
    Dashing[{0.004`, 0.02`, 0.03`, 0.02`}]}];
a3 = Show[bb1, bb2, bb3, bb4, Frame -> True, GridLines -> Automatic, 
  LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
  FrameLabel -> {{Error["%"],}, {Offset[km],}}, 
  PlotRange -> {{0, 0.1}, {-0.7, 0.15}}];
mm = GraphicsRow[{a1, a2, a3}, AspectRatio -> 0.4, Spacings -> 0, 
  ImageSize -> {1000, 300}];
Export["junk_ma.eps", mm, "EPS"];
