f[t_,p_,h_]:=(t + Sqrt[t^2 + 4 h^2 p^2])/2;
t2[t_,p_,s_,h_]:=t f[t,p,h] + s h^2;
tp[x_,p_,s_,h_]:=t2[p x,p,s,h]/Sqrt[f[t2[p x,p,s,h],p f[p x,p,h],h]];
xp[x_,p_,s_,h_]:=x + h^2 p/f[p x,p,h] * 
	(1 - f[p x,p,h]^2/f[t2[p x,p,s,h],p f[p x,p,h],h]);
fp[v_]:=ParametricPlot[{xp[x,0.5,1-1/v^2,2], -v tp[x,0.5,1-1/v^2,2]},
	{x,0,2.1},PlotStyle->Thickness[0.01]];
Show[fp[0.8],fp[0.9],fp[1],fp[1.1],fp[1.2]];
rp[x_]:=ParametricPlot[{xp[x,0.5,1-1/v^2,2], -v tp[x,0.5,1-1/v^2,2]},
	{v,0.8,1.2},PlotStyle->Dashing[{0.01,0.01}]];
Show[rp[0],rp[0.5],rp[1],rp[1.5],rp[2]];
Show[%6,%8,PlotRange->{{0,2.5},{0,-2}},
	Ticks->{{0.5,1,1.5,2,2.5},{{-0.5,"0.5"},{-1,"1"},{-1.5,"1.5"}}}];
d[x_]:=Sqrt[1 + x^2];
td[x_,s_,h_]:=t2[d[x],x/d[x],s,h]/Sqrt[f[t2[d[x],x/d[x],s,h],
	x/d[x] f[d[x],x/d[x],h],h]];
xd[x_,s_,h_]:=x + h^2 x/(d[x] f[d[x],x/d[x],h]) * 
	(1 - f[d[x],x/d[x],h]^2/f[t2[d[x],x/d[x],s,h],
	x/d[x] f[d[x],x/d[x],h],h]);
fd[v_]:=ParametricPlot[{xd[x,1-1/v^2,2], -v td[x,1-1/v^2,2]},{x,-2.1,2.1},
	PlotStyle->Thickness[0.01]];
Show[fd[0.9],fd[0.95],fd[1],fd[1.05],fd[1.1]];
rd[x_]:=ParametricPlot[{xd[x,1-1/v^2,2],-v td[x,1-1/v^2,2]},{v,0.9,1.1},
	PlotStyle->Dashing[{0.01,0.01}]];
Show[rd[-1.75],rd[-1.25],rd[-0.75],rd[-0.25],
	rd[0.25],rd[0.75],rd[1.25],rd[1.75]];
Show[%14,%16,PlotRange->{{-2.5,2.5},{0,-3}},
	Ticks->{{-2,-1,0,1,2},{{-1,"1"},{-2,"2"},{-3,"3"}}}];
Show[GraphicsArray[{%17,%9}],AspectRatio->1/2];
Display["junk_ma.eps", %, "EPS", ImageSize->{396,198}];

