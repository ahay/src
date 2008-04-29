d={1,2,3,5,5};
d={1,1,1,2,5};
d={1,1,2,3,5};

A [r_,n_,m_] := Sum [Log [1 + (d[[i]] - m)^2/r^2],{i,1,n}];

a1 [r_] := Plot [A[r,1,m],{m,0,2},PlotRange->{0.,Automatic}];
an [r_] := Plot [A[r,5,m],{m,0,6},PlotRange->{0.,Automatic}];

a [r_] := Show[GraphicsArray[{{a1 [r]},{an [r]}}]];

gr = Show[GraphicsArray[{a[2],a[1],a[0.2]}]];

Export["junk_ma.eps", gr, "EPS", ImageSize->1000];




