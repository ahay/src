vtt[a_,vn_,vx_]:=1/Sqrt[Cos[a]^2 + 
	  	(1/vn - 1/vx) Cos[a]^2 Sin[a]^2 + 1/vx Sin[a]^2];
ptt[vn_,vx_]:=ParametricPlot[{vtt[a,vn,vx] Sin[a],-vtt[a,vn,vx] Cos[a]},
		{a,-Pi/2,Pi/2},
		PlotStyle->Thickness[0.01]];
izo[v_]:=ParametricPlot[{Sqrt[v] Sin[a],-Sqrt[v] Cos[a]},
		{a,-Pi/2,Pi/2},
		PlotStyle->Dashing[{0.01,0.01}]];
$DefaultFont={"Courier",9};
Show[ptt[1,1],
	Ticks->None,
	PlotLabel->"a",
	AspectRatio->Automatic];
Show[ptt[1.4,1.4],izo[1],izo[1.4],
	Ticks->None,
	PlotLabel->"b",
	AspectRatio->Automatic];
Show[ptt[0.6,1.4],izo[1],izo[1.4],
	Ticks->None,
	PlotLabel->"d",
	AspectRatio->Automatic];
Show[ptt[1.,1.4],izo[1],izo[1.4],
	Ticks->None,
	PlotLabel->"c",
	AspectRatio->Automatic];
Show[GraphicsArray[{{%5,%6},{%8,%7}}]];
Display["junk_ma.eps", %, "EPS",ImageSize->400];
