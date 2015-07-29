angleValue[v_, vz_, n_, khx_, kmx_, w_] :=
  Block[{hx, mx},
   hx = khx/w;
   mx = kmx/w;
   ccc = 1/(
     2 v^46) (hx - mx)^4 (hx + 
       mx)^4 (-2 + (hx - mx)^2 n v^2)^3 (-2 + (hx + 
          mx)^2 n v^2)^4 (-4 + (hx - mx)^2 (1 + 
          4 n) v^2)^2 vz^12 (-64 + 32 (hx^2 + mx^2) (1 + 4 n) v^2 - 
       4 (hx^2 - mx^2)^2 (1 + 4 n)^2 v^4 + (hx - mx)^2 (hx + 
          mx)^2 (4 - 
          4 (hx^2 + mx^2) n v^2 + (hx^2 - mx^2)^2 n^2 v^4) vz^4)^2;
   bbb = -(1/
      v^46) (hx - mx)^4 (hx + 
       mx)^4 (-2 + (hx - mx)^2 n v^2)^3 (-2 + (hx + 
          mx)^2 n v^2)^4 (-4 + (hx - mx)^2 (1 + 
          4 n) v^2)^2 vz^12 (64 + 
       4 (1 + 4 n) v^2 (-8 (hx^2 + mx^2) + (hx^2 - mx^2)^2 (1 + 
             4 n) v^2) + (hx - mx)^2 (hx + mx)^2 (4 - 
          4 (hx^2 + mx^2) n v^2 + (hx^2 - 
             mx^2)^2 n^2 v^4) vz^4) (-8 + (hx - mx)^2 (-2 vz^2 + 
          v^2 (2 + n (8 + (hx - mx)^2 vz^2)))) (-8 + (hx + 
          mx)^2 (-2 vz^2 + v^2 (2 + n (8 + (hx + mx)^2 vz^2))));
   aaa = 1/(
     2 v^46 ) (hx - mx)^4 (hx + 
       mx)^4 (-2 + (hx - mx)^2 n v^2)^3 (-2 + (hx + 
          mx)^2 n v^2)^4 (-4 + (hx - mx)^2 (1 + 
          4 n) v^2)^2 vz^12 (-8 + (hx - mx)^2 (-2 vz^2 + 
          v^2 (2 + n (8 + (hx - mx)^2 vz^2))))^2 (-8 + (hx + 
          mx)^2 (-2 vz^2 + v^2 (2 + n (8 + (hx + mx)^2 vz^2))))^2;
   (*cc2=If[Abs[hx]<Abs[mx] ,(-bbb-Sqrt[bbb^2-4 aaa ccc])/(
   2 aaa),(-bbb+Sqrt[bbb^2-4 aaa ccc])/(2 aaa)];*)
   
   cc2 = If[Abs[hx] < 0.01  && Abs[mx] == 0.01, 1, 
     If[Abs[hx] > Abs[mx] , (
      2 ccc)/(-bbb - Sqrt[bbb^2 - 4 aaa ccc]), (
      2 ccc)/(-bbb + Sqrt[bbb^2 - 4 aaa ccc])]];
   thta = Sign[hx] ArcCos[Sqrt[cc2]];
   Return[{Sqrt[cc2], thta}]];
<< PlotLegends`;
a = ShowLegend[
  ContourPlot[
   angleValue[2, 2.0, 0.0, phx, pmx, 60][[2]], {phx, -35, 
    35}, {pmx, -35, 35}, Contours -> 20, 
   FrameLabel -> {{Subscript[k, mx]["\!\(\*SuperscriptBox[\"km\", 
RowBox[{\"-\", \"1\"}]]\)"],}, {Subscript[k, hx][
       "\!\(\*SuperscriptBox[\"km\", 
RowBox[{\"-\", \"1\"}]]\)"],}}, 
   LabelStyle -> Directive[Bold, FontFamily -> "Times", 14], 
   ColorFunction -> "Rainbow", 
   PlotRange -> {{-35, 35}, {-35, 35}, {-Pi/2, Pi/2}}], {ColorData[
      "Rainbow"][1 - #1] &, 13, "\[Pi]/2", "-\[Pi]/2", 
   LegendPosition -> {0.9, -0.7}, LegendSize -> 1.3, 
   LegendLabel -> 
    Style["\!\(\*SuperscriptBox[\"\[Theta]\", \"o\"]\)", 14]}];
Display["junk_ma.eps", a, "EPS"];
