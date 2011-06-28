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
  A = -4 eta; B = (1 + 8 eta + 8 eta^2)/(1 + 2 eta);
  CC = 1/(1 + 2 eta)^2;
  ttt = (Sqrt[
      x^2/v^2 + z^2/
       vv^2] (-s^2 (v - vv) (v + vv) (vv^2 x^2 + 
           v^2 z^2)^2 (-vv^2 x^4 + 2 (v - vv) (v + vv) x^2 z^2 + 
           v^2 z^4) + 
        eta vv^4 x^4 ((2 + eta) vv^4 x^4 + 
           4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 2 v^4 z^4) + 
        2 s x z (vv^2 x^2 + 
           v^2 z^2) (vv^4 ((1 - 5 eta) v^2 + (-1 + eta) vv^2) x^4 + 
           2 v^2 vv^2 (v^2 - (1 + 2 eta) vv^2) x^2 z^2 + 
           v^4 (v - vv) (v + vv) z^4)))/(eta vv^4 x^4 ((2 + 
            3 eta) vv^4 x^4 + 4 (1 + 3 eta) v^2 vv^2 x^2 z^2 + 
         2 v^4 z^4) - 
      2 s x z (vv^2 x^2 + v^2 z^2) ((1 + eta) vv^6 x^4 - v^6 z^4 + 
         v^4 vv^2 z^2 (-2 x^2 + z^2) + 
         v^2 vv^4 x^2 ((-1 + 3 eta) x^2 + 2 (1 + 2 eta) z^2)) - 
      s^2 (vv^2 x^2 + v^2 z^2)^2 (vv^4 x^4 + v^4 z^4 - 
         v^2 vv^2 (x^4 + z^4)));
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
  ttt3 = -(2 (vv^2 x^2 + v^2 z^2)^6 Sqrt[
       x^2/v^2 + z^2/
        vv^2])/((vv^2 x^2 + 
         v^2 z^2)^2 ((-2 + (-2 + eta) eta) vv^8 x^8 + 
         4 (-1 + eta) (2 + 3 eta) v^2 vv^6 x^6 z^2 - 
         2 (6 + eta) v^4 vv^4 x^4 z^4 - 8 v^6 vv^2 x^2 z^6 - 
         2 v^8 z^8) - 
      2 s (v - vv) (v + vv) x z (vv^2 x^2 + 
         v^2 z^2) ((1 + 2 eta) vv^8 x^8 + 
         4 (1 + eta - 3 eta^2) v^2 vv^6 x^6 z^2 + 
         2 (3 + eta) v^4 vv^4 x^4 z^4 + 4 v^6 vv^2 x^2 z^6 + 
         v^8 z^8) - 
      s^2 (v - vv) (v + vv) (v^10 z^10 (2 x^2 + z^2) + 
         v^8 vv^2 x^2 z^8 (7 x^2 + 2 z^2) + 
         2 v^6 vv^4 x^4 z^6 ((4 + 3 eta) x^2 + (-1 + eta) z^2) + 
         2 v^4 vv^6 x^6 z^4 ((1 + (5 - 18 eta) eta) x^2 - (4 + eta + 
               6 eta^2) z^2) - 
         vv^10 x^10 ((1 + 2 eta) x^2 + (2 + 3 eta (2 + eta)) z^2) + 
         v^2 vv^8 x^8 z^2 ((-2 + eta (2 + 15 eta)) x^2 + (-7 + 
               2 eta (-5 + 18 eta)) z^2)));
  ttt4 = Sqrt[
   x^2/v^2 + z^2/vv^2 + (A X^4)/(
    v^2 (Z^2 + B X^2 + Sqrt[Z^4 + 2 B Z^2  X^2 + CC X^4]))];
  ttT = Sqrt[(((-2 + s^2) v^2 - s^2 vv^2)^2 z^2)/(4 v^4 vv^2) + 
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
  ttF = (Sqrt[
      x^2/v^2 + z^2/
       vv^2] (-s^2 (v - vv) (v + vv) (x^2 + z^2) (-vv^2 x^2 + 
           v^2 z^2) (vv^2 x^2 + v^2 z^2)^2 + 
        eta vv^4 x^4 ((1 + 2 eta) vv^4 x^4 + 
           2 (1 + 6 eta) v^2 vv^2 x^2 z^2 + v^4 z^4) + 
        s x z (vv^2 x^2 + 
           v^2 z^2) (vv^4 ((1 - 5 eta) v^2 + (-1 + eta) vv^2) x^4 + 
           2 v^2 vv^2 (v^2 - (1 + 2 eta) vv^2) x^2 z^2 + 
           v^4 (v - vv) (v + vv) z^4)))/(eta vv^4 x^4 ((1 + 
            3 eta) vv^4 x^4 + 2 (1 + 6 eta) v^2 vv^2 x^2 z^2 + 
         v^4 z^4) - 
      s x z (vv^2 x^2 + v^2 z^2) ((1 + eta) vv^6 x^4 - v^6 z^4 + 
         v^4 vv^2 z^2 (-2 x^2 + z^2) + 
         v^2 vv^4 x^2 ((-1 + 3 eta) x^2 + 2 (1 + 2 eta) z^2)) - 
      s^2 (vv^2 x^2 + v^2 z^2)^2 (vv^4 x^4 + v^4 z^4 - 
         v^2 vv^2 (x^4 + z^4)));
  ttF2 = (Sqrt[
      x^2/v^2 + z^2/
       vv^2] ((1 + 2 eta) vv^6 x^7 + (3 + 14 eta) v^2 vv^4 x^5 z^2 + 
        3 (1 + 4 eta) v^4 vv^2 x^3 z^4 + v^6 x z^6 - 
        s^3 (v - vv) (v + vv) z (-vv x^2 + v z^2) (vv x^2 + 
           v z^2) ((3 v^2 + vv^2) x^2 + 4 v^2 z^2) + 
        s z (vv^4 ((2 - 9 eta) v^2 + (2 + eta) vv^2) x^6 + 
           4 v^2 vv^2 ((1 - 3 eta) v^2 + (2 + eta) vv^2) x^4 z^2 + 
           2 v^4 (v^2 + 5 vv^2) x^2 z^4 + 4 v^6 z^6) + 
        s^2 x ((1 + 3 eta) (v - vv) vv^4 (v + vv) x^6 + 
           vv^2 (-2 v^4 + (1 - 15 eta) v^2 vv^2 - (-1 + 
                 eta) vv^4) x^4 z^2 - 
           v^2 (3 v^4 + 3 (1 + 4 eta) v^2 vv^2 + 
              2 (-3 + 2 eta) vv^4) x^2 z^4 + 
           5 v^4 (-v^2 + vv^2) z^6)))/((vv^2 x^2 + 
        v^2 z^2) (vv^4 x^4 (x + 3 eta x + s z) + 
        v^4 z^3 (3 s x^2 + x z + 4 s z^2) + 
        v^2 vv^2 x^2 z (3 s x^2 + 2 (1 + 6 eta) x z + 5 s z^2)));
  ttF3 = Sqrt[
     x^2/v^2 + z^2/vv^2] + ((-v^2 + vv^2) x z Sqrt[
      x^2/v^2 + z^2/vv^2])/(vv^2 x^2 + v^2 z^2)
      s + (-vv^4 x^4 Sqrt[
      x^2/v^2 + z^2/vv^2])/(vv^2 x^2 + v^2 z^2)^2 eta + (
     Sqrt[x^2/v^2 + z^2/
       vv^2] (-vv^4 x^4 - v^4 z^4 + 
        v^2 vv^2 (x^4 + z^4)))/ ((vv^2 x^2 + v^2 z^2)^2) s^2 + (
     3 vv^6 x^6 (vv^2 x^2 + 4 v^2 z^2) Sqrt[
      x^2/v^2 + z^2/vv^2])/ ((vv^2 x^2 + v^2 z^2)^4)
      eta^2 + (-vv^4 x^3 z ((3 v^2 + vv^2) x^2 + 4 v^2 z^2) Sqrt[
      x^2/v^2 + z^2/vv^2])/(vv^2 x^2 + v^2 z^2)^3 eta s;
  Return[{ttt, ttt2, ttt3, ttt4, ttT, ttT2, ttTSena, ttF, ttF2, 
    ttF3}]];
<< PlotLegends`;
bb = ShowLegend[
  ContourPlot[
   LinVeikdsTTI[2, 0.05, eta (1 + 2*0.05) + 0.05, 2, angle, 2][[2]] - 
    LinVeikdsTTI[2, 0.05, 0.05, 2, 0, 2][[2]], {angle, 0, 
    90 Pi/180}, {eta, -0.1, 1.0}, Contours -> 20, 
   FrameLabel -> {{\[Eta],}, {\[Theta][rad], "b"}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow", 
   PlotRange -> {{0, Pi/2}, {0, 1}, {-1, 0}}], {ColorData["Rainbow"][
     1 - #1] &, 13, "0.0", "-1", LegendPosition -> {1.1, -0.7}, 
   LegendSize -> {0.3, 1.5}, 
   LegendLabel -> Style["\[CapitalDelta]t (s)", 14]}];
aa = ShowLegend[
  ContourPlot[
   LinVeikdsTTI[1, 0.05, eta (1 + 2*0.05) + 0.05, 2, angle, 2][[2]] - 
    LinVeikdsTTI[1, 0.05, 0.05, 2, 0, 2][[2]], {angle, 0, 
    90 Pi/180}, {eta, -0.1, 1.0}, Contours -> 20, 
   FrameLabel -> {{\[Eta],}, {\[Theta][rad], "a"}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow", 
   PlotRange -> {{0, Pi/2}, {0, 1}, {-0.3, 0}}], {ColorData[
      "Rainbow"][1 - #1] &, 13, "0.0", "-0.3", 
   LegendPosition -> {1.1, -0.7}, LegendSize -> {0.3, 1.5}, 
   LegendLabel -> Style["\[CapitalDelta]t (s)", 14]}];
cc = ShowLegend[
  ContourPlot[
   LinVeikdsTTI[4, 0.05, eta (1 + 2*0.05) + 0.05, 2, angle, 2][[2]] - 
    LinVeikdsTTI[4, 0.05, 0.05, 2, 0, 2][[2]], {angle, 0, 
    90 Pi/180}, {eta, -0.1, 1.0}, Contours -> 20, 
   FrameLabel -> {{\[Eta],}, {\[Theta][rad], "c"}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow", 
   PlotRange -> {{0, Pi/2}, {0, 1}, {-2.0, 0}}], {ColorData[
      "Rainbow"][1 - #1] &, 13, "0.0", "-2.0", 
   LegendPosition -> {1.1, -0.7}, LegendSize -> {0.3, 1.5}, 
   LegendLabel -> Style["\[CapitalDelta]t (s)", 14]}];
mm = GraphicsRow[{aa, bb, cc}, AspectRatio -> 0.4, Spacings -> 0, 
  ImageSize -> {1000, 300}];
Export["junk_ma.eps", mm, "EPS"];
