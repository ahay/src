f[a_,y_]:=(2*(-(1 + a^2*Cos[y]^2 + Tan[y]^2 + 
 	(a^4*Cos[y]^4*Sin[y]^2)/(1 + a^2*Tan[y]^2))^(1/2) + 
 	((1 + (-a + Tan[y])^2)^(1/2) + 
 	(1 + (a + Tan[y])^2)^(1/2))/2))/
 	((1 + (-a + Tan[y])^2)^(1/2) + (1 + (a + Tan[y])^2)^(1/2));
$DefaultFont={"Courier",9};
lo[a_]:=Plot[f[a,y Pi/180],{y,0,90},
	AxesLabel->{"angle",None},
	PlotLabel->"Relative Error"];
lo[1];
Display["junk_ma.eps", %, "EPS"];
