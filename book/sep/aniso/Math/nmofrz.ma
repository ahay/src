vbh[a_,vn_,vx_]:=1/Sqrt[Cos[a]^2 + 
		(1/vn - 1/vx) Cos[a]^2 Sin[a]^2 + 1/vx Sin[a]^2];
pbh[vn_,vx_]:=ParametricPlot[{vbh[a,vn,vx] Sin[a],-vbh[a,vn,vx] Cos[a]},
		{a,-Pi/2,Pi/2},
		PlotStyle->Thickness[0.01]];
vvz[a_,vn_,vx_]:=1/(Cos[a] (1 - vn/vx) + 
		Sqrt[Sin[a]^2 1/vx + Cos[a]^2 vn^2/vx^2]);
pvz[vn_,vx_]:=ParametricPlot[{vvz[a,vn,vx] Sin[a],-vvz[a,vn,vx] Cos[a]},
		{a,-Pi/2,Pi/2},
		PlotStyle->Dashing[{0.01,0.01}]];
tbh[x_,vx_,vn_]:=Sqrt[1 + x^2/vn - x^4/(1 + x^2) (1/vn - 1/vx)];
pltbh[vx_,vn_]:=Plot[-tbh[x,vx,vn],
		{x,-3,3},
		PlotStyle->Thickness[0.01]];
tvz[x_,vx_,vn_]:=(1 - vn/vx) + vn/vx Sqrt[1 + vx x^2 /vn^2];
pltvz[vx_,vn_]:=Plot[-tvz[x,vx,vn],
		{x,-3,3},
		PlotStyle->Dashing[{0.01,0.01}]];
$DefaultFont={"Courier",9};
Show[pbh[1.2,1.4],pvz[1.2,1.4],
	Frame->True,
	FrameTicks->{None,None},
	PlotLabel->"a",
	FrameLabel->{None,None},
	AspectRatio->90/162];
Show[pltbh[1.2,1.4],pltvz[1.2,1.4],
		Frame->True,
		PlotLabel->"b",
		PlotRange->{-3.,0.},
		FrameLabel->{None,None},
		FrameTicks->{None,None},
		AspectRatio->90/162];
Show[GraphicsArray[{%%,%}]];
Display["junk_ma.eps", %, "EPS",ImageSize->360];
