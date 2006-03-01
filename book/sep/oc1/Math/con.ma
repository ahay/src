(0.25 - h^2/r)*(1 - r) == x^2;
Solve[%,r];
t1[x_,h_]:=Sqrt[2.*(0.25 + 1.*h^2 - 1.*x^2 + 
    1.*(0.0625 - 0.5*h^2 + 1.*h^4 - 0.5*x^2 - 2.*h^2*x^2 + 
        1.*x^4)^(1/2))];
t2[x_,h_]:=Sqrt[2.*(0.25 + 1.*h^2 - 1.*x^2 - 
    1.*(0.0625 - 0.5*h^2 + 1.*h^4 - 0.5*x^2 - 2.*h^2*x^2 + 
        1.*x^4)^(1/2))];
p1[h_]:=Plot[-t1[x,h],{x,-0.999999 Abs[h-0.5],0.999999 Abs[h-0.5]}];
p2[h_]:=Plot[-t2[x,h],{x,-0.999999 Abs[h-0.5],0.999999 Abs[h-0.5]}];
Show[p1[0],p1[0.1],p1[0.2],p1[0.3],p1[0.4]];
Show[p2[0.6],p2[0.7],p2[0.8],p2[0.9],p2[0.9],p2[1]];
Show[%7,%8,PlotRange->{0,-2.25}, 
	AspectRatio->1, Frame->True,
	FrameLabel->{"midpoint",None},
	FrameTicks->{Automatic,
	{{-0.5,"0.5"},{-1,"1"},{-1.5,"1.5"},{-2,"2"}},None,None}];
Show[p2[0.1],p2[0.2],p2[0.3],p2[0.4]];
Show[p1[0.6],p1[0.7],p1[0.8],p1[0.9],p1[0.9],p1[1]];
Show[%10,%11,%9,PlotRange->{0,-2.25},
	AspectRatio->1, Frame->True,
	FrameLabel->{"midpoint",None},
	FrameTicks->{Automatic,
	{{-0.5,"0.5"},{-1,"1"},{-1.5,"1.5"},{-2,"2"}},None,None}];
Show[GraphicsArray[{%12,%9}],AspectRatio->1/2];
Display["junk_ma.eps", %, "EPS", ImageSize -> 432];
