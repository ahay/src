GS = {A -> 3.41, B -> 0.54, C -> 2.27, D -> 1.07};
GS = {A -> 14.47, C -> 9.57, B -> 2.28, D -> 4.51}; 
PhaSV[a_, b_] := 
  1/2((A + B) b^2 + (C + B) a^2 - 
        Sqrt[((A - B) b^2 - (C - B) a^2)^2 + 4 (D + B)^2 a^2 b^2]);
PhaP[a_, b_] := 
  1/2((A + B) b^2 + (C + B) a^2 + 
        Sqrt[((A - B) b^2 - (C - B) a^2)^2 + 4 (D + B)^2 a^2 b^2]);
ParametricPlot[{Sqrt[PhaSV[Cos[a], Sin[a]]] Sin[a], 
      Sqrt[PhaSV[Cos[a], Sin[a]]]Cos[a]} /. GS, {a, 0, 2 Pi}, 
  AspectRatio -> Automatic, PlotStyle->AbsoluteThickness[2]];
ParametricPlot[{Sqrt[PhaP[Cos[a], Sin[a]]] Sin[a], 
      Sqrt[PhaP[Cos[a], Sin[a]]]Cos[a]} /. GS, {a, 0, 2 Pi}, 
  AspectRatio -> Automatic, PlotStyle->AbsoluteThickness[2]];
Show[%,%%,Frame->True,
FrameLabel->{"Horizontal component (km/s)","Vertical component(km/s)",
None,None},
     PlotLabel->"Phase Velocity Profiles"];
Display["junk_ma.eps",%,"EPS"];
