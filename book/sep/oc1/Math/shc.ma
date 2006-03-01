f[r_]:=Sqrt[(r-0.5) (r+0.5)];
bf[r_]:=Sqrt[(0.5-r) (r+0.5)];
t[s_,r_]:=Sqrt[(r-s)^2-(f[s]+f[r])^2];
bt[s_,r_]:=Sqrt[(r-s)^2+(bf[s]+bf[r])^2];
pt[s_]:=Plot[-t[s,r],{r,-0.5,0.5}];
Show[pt[-0.4],pt[-0.2],pt[0.],pt[0.2],pt[0.4]];
ptm[s_]:=Plot[-bt[s,r],{r,-1.5,-0.5}];
Show[ptm[0.6],ptm[0.8],ptm[1.0]];
ptp[s_]:=Plot[-bt[s,r],{r,0.5,1.5}];
Show[ptp[-0.6],ptp[-0.8],ptp[-1.0]];
$DefaultFont={"Courier",9};
Show[%6,%8,%10,
	Graphics[Text["-0.4",{-0.625,-0.3}]],
	Graphics[Text["-0.2",{-0.625,-0.55}]],
	Graphics[Text["0",{-0.575,-0.7}]],
	Graphics[Text["0.2",{-0.625,-0.85}]],
	Graphics[Text["0.4",{-0.625,-0.95}]],
	Graphics[Text["0.6",{-0.4,-1.075}]],
	Graphics[Text["0.8",{-0.4,-1.15}]],
	Graphics[Text["1.0",{-0.4,-1.225}]],
	Graphics[Text["-0.6",{0.375,-1.075}]],
	Graphics[Text["-0.8",{0.375,-1.15}]],
	Graphics[Text["-1.0",{0.375,-1.225}]],
	PlotRange->{0,-1.5},
	AxesLabel->{"r",""},
	AspectRatio->1/2,
	Ticks->{Automatic,
	{{-0.25,"0.25"},{-0.5,"0.5"},{-0.75,"0.75"},
	{-1,"1"},{-1.25,"1.25"}}}];
Display["junk_ma.eps", %,"EPS"];
