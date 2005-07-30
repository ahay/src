sinc [x_]          := Sin[Pi x]/(Pi x);
wind[x_, n_Integer]:=    (Pi x)/(2 n Tan[(Pi x)/(2 n)]);
muir[x_, n_Integer]:= Sin[Pi x]/(2 n Tan[(Pi x)/(2 n)]);
ps[n_Integer]:=Plot [sinc [x],  {x,-n,n},PlotStyle->{Thickness[0.01]},PlotRange->{-0.25,1.05},Frame->True];
pw[n_Integer]:=Plot [wind [x,n],{x,-n,n},PlotStyle->{Thickness[0.01]}];
pm[n_Integer]:=Plot [muir [x,n],{x,-n,n},PlotStyle->{Thickness[0.01]},PlotRange->{-0.25,1.05},Frame->True];
Show[GraphicsArray[{{ps[2],pw[2],pm[2]},
		    {ps[6],pw[6],pm[6]}}]];
Display["junk_ma.eps",%,"EPS",ImageSize->648];



