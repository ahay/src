fd[w_, k_] := (-27 I + w (-12 + w (-18 I + 5 w)) + 
	(27 I + 12 w + w^3) Cos[k])/
        (-81 I + w (24 + w (-12 I + 5 w)) + 
        (81 I + w (-24 + w (-6 I + w))) Cos[k]);
F[e_] := Sqrt[(1 + Sqrt[1 + e^2])/
	(2 Sqrt[1 + e^2])] Exp[(1 - Sqrt[1 + e^2])/2];
f[e_] := 1/2 (1 - Sqrt[1 + e^2] + Log[(1 + Sqrt[1 + e^2])/2]);
as[w_, k_] := F[2 k /w]/F[0] Exp[I w  (f[2 k /w]  -  f[0])];
ex[w_, k_] := Hypergeometric0F1[1 - (1 +  I w)/2, -k^2/4];
Needs["Graphics`Legend`"];
plarg[w_] := Plot[{-Arg[ex[w, k]], -Arg[as[w, k]], 
                   -Arg[fd[w, k]]}, {k, 0, Pi},
	Frame->True,PlotLabel->SequenceForm["Phase for frequency=",w],
	FrameLabel->{"Wavenumber (radians)",None,None,None},
	PlotLegend->{"Exact","Asymptotic","Finite-difference"},
	LegendPosition->{-0.7,0.2}, LegendSize->{0.9,0.3},
	PlotStyle->{{GrayLevel[0.5],Dashing[{0.01}]},
		    {GrayLevel[0.5],Dashing[{0.03}]},
	            {}}];
Show[GraphicsArray[{plarg[1], plarg[10]}],AspectRatio->1/(2 GoldenRatio)];
Display["junk_ma.eps",%,"EPS",ImageSize->864];
