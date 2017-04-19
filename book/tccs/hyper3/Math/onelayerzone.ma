Homogeneous Orthorhombic Layer
ClearAll[c11,c12,c13,c12,c23,c22,c33,c44,c55,c66,n1,n2,n3,two,trace,G,H,angle,vhelbig,vpphase,vecdiff,vpgroup,vpg,pn12,pn22,pn32]

Exact velocity
p1=.;p2=.;p3=.;
mat = {{c11,c12,c13,0,0,0},{c12,c22,c23,0,0,0},{c13,c23,c33,0,0,0},{0,0,0,c44,0,0},{0,0,0,0,c55,0},{0,0,0,0,0,c66}};
matn = {{n1,0,0,0,n3,n2},{0,n2,0,n3,0,n1},{0,0,n3,n2,n1,0}} ;
chris = matn.mat.Transpose[matn];
chrisp = chris /.{n1->p1,n2->p2,n3->p3};
chrisp//MatrixForm
det = Det[chrisp-IdentityMatrix[3]];
(c11 p1^2+c66 p2^2+c55 p3^2	c12 p1 p2+c66 p1 p2	c13 p1 p3+c55 p1 p3
c12 p1 p2+c66 p1 p2	c66 p1^2+c22 p2^2+c44 p3^2	c23 p2 p3+c44 p2 p3
c13 p1 p3+c55 p1 p3	c23 p2 p3+c44 p2 p3	c55 p1^2+c44 p2^2+c33 p3^2

)
Hamiltonian to derive coefficients
F = det(*Hamiltonian*)
(c13 p1 p3+c55 p1 p3) ((c12 p1 p2+c66 p1 p2) (c23 p2 p3+c44 p2 p3)-(c13 p1 p3+c55 p1 p3) (-1+c66 p1^2+c22 p2^2+c44 p3^2))-(c23 p2 p3+c44 p2 p3) (-(c12 p1 p2+c66 p1 p2) (c13 p1 p3+c55 p1 p3)+(c23 p2 p3+c44 p2 p3) (-1+c11 p1^2+c66 p2^2+c55 p3^2))+(-1+c55 p1^2+c44 p2^2+c33 p3^2) (-(c12 p1 p2+c66 p1 p2)^2+(-1+c66 p1^2+c22 p2^2+c44 p3^2) (-1+c11 p1^2+c66 p2^2+c55 p3^2))
DON'T FORGET TO CHECK WHICH ROOT IS THE RIGHT ONE FOR qP
p3sqf=p3sq/.Solve[0==F/.p3->Sqrt[p3sq],p3sq][[1]] ;(*Orthorhombic qP*)
(*p3sqf=p3sq/.Solve[0F/.p3Sqrt[p3sq],p3sq][[2]];*)(*VTI qP*)
p3sqf1=p3sqf/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
DON'T FORGET A FACTOR OF  ΔZ  MULTIPLYING THE DERIVATIVES
dxdσ = D[F,p1]//Simplify;
dydσ = D[F,p2]//Simplify;
dzdσ = D[F,p3]//Simplify;
(*dx =Δd dxdσ /dzdσ/.{p3Sqrt[p3sqf]}/.{c119,c335.9375,c551.6,c132.25,c229.84,c442,c123.6,c232.4,c662.182};
dy = Δd dydσ /dzdσ/.{p3Sqrt[p3sqf]}/.{c119,c335.9375,c551.6,c132.25,c229.84,c442,c123.6,c232.4,c662.182};
dt =Δd dtdσ /dzdσ/.{p3Sqrt[p3sqf]}/.{c119,c335.9375,c551.6,c132.25,c229.84,c442,c123.6,c232.4,c662.182};*)
dx =Δd dxdσ /dzdσ /.{p3->Sqrt[p3sqf1]}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
dy = Δd dydσ /dzdσ /.{p3->Sqrt[p3sqf1]}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
dt =Δd ((p1 dxdσ+p2 dydσ)/dzdσ +p3)/.{p3->Sqrt[p3sqf1]}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
depth = 1;
maxx = 0.8/Sqrt[c11]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
maxy = 0.8/Sqrt[c22]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
maxr = R maxx maxy/Sqrt[(maxy^2 Cos[theta]^2 + maxx^2 Sin[theta]^2)]; (*Ellipse in polar*)
xoffset = (2 dx)/.{Δd->depth}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};//Simplify
yoffset = (2dy)/.{Δd->depth}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};//Simplify
exacttime = (2dt)/.{Δd->depth}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};
NMO Ellipse
rortho = Sqrt[xoffset^2+yoffset^2];
c = xoffset/rortho;
s = yoffset/rortho;

δ2=((c13+c55)^2-(c33-c55)^2)/(2*c33*(c33-c55));
δ1=((c23+c44)^2-(c33-c44)^2)/(2*c33*(c33-c44));
Q2 = (1/(c33(1+2δ2)))c^2 +(1/(c33(1+2δ1)))s^2;

orthoell = t0^2 + Q2 rortho^2/.t0->2 depth/Sqrt[c33] /.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
nmorms[radius_,angle_]:= 100Abs[Sqrt[orthoell]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};

Map[SetOptions[#,Prolog->{{EdgeForm[],Texture[{{{0,0,0,0}}}],Polygon[#,VertexTextureCoordinates->#]&[{{0,0},{1,0},{1,1}}]}}]&,{Graphics3D,ContourPlot3D,ListContourPlot3D,ListPlot3D,Plot3D,ListSurfacePlot3D,ListVectorPlot3D,ParametricPlot3D,RegionPlot3D,RevolutionPlot3D,SphericalPlot3D,VectorPlot3D,BarChart3D}];

Sripanich and Fomel (2014)
Q32 = (c11(c33-c55))/((c13+c55)^2+c55(c33-c55));
Q12 = Q32/(0.8373406986+0.1580976683 Q32);
Q31 = (c22(c33-c44))/((c23+c44)^2+c44(c33-c44));
Q21 = Q31/(0.8373406986+0.1580976683 Q31);
Q13 = (c22(c11-c66))/((c12+c66)^2+c66(c11-c66));
Q23 = Q13;

(*Q12=(c33(c11-c55))/((c13+c55)^2+c55(c11-c55));
Q21=(c33(c22-c44))/((c23+c44)^2+c44(c22-c44));
Q23=(c11(c22-c66))/((c12+c66)^2+c66(c22-c66));*)

S32=((c11-c33) (-1+Qh2c) (-1+Qv2c)^2)/(2 c11 (-1+Qh2c^2+2 Qv2c+Qh2c Qv2c (-1+(-2+Qv2c) Qv2c))-2 c33 (-1+Qh2c^2-2 Qh2c Qv2c+Qv2c (2+(-1+Qv2c)^2 Qv2c)))/.{Qv2c->Q32,Qh2c->Q12};
S12=((c11-c33) (-1+Qh2c)^2 (-1+Qv2c))/(2 (c33+c11 (-1+Qh2c (2+(-1+Qh2c)^2 Qh2c)-2 Qh2c Qv2c+Qv2c^2)-c33 (2 Qh2c+Qh2c (-1+(-2+Qh2c) Qh2c) Qv2c+Qv2c^2)))/.{Qv2c->Q32,Qh2c->Q12};
S31=((c22-c33) (-1+Qh1c) (-1+Qv1c)^2)/(2 c22 (-1+Qh1c^2+2 Qv1c+Qh1c Qv1c (-1+(-2+Qv1c) Qv1c))-2 c33 (-1+Qh1c^2-2 Qh1c Qv1c+Qv1c (2+(-1+Qv1c)^2 Qv1c)))/.{Qv1c->Q31,Qh1c->Q21};
S21=((c22-c33) (-1+Qh1c)^2 (-1+Qv1c))/(2 (c33+c22 (-1+Qh1c (2+(-1+Qh1c)^2 Qh1c)-2 Qh1c Qv1c+Qv1c^2)-c33 (2 Qh1c+Qh1c (-1+(-2+Qh1c) Qh1c) Qv1c+Qv1c^2)))/.{Qv1c->Q31,Qh1c->Q21};
S13=((c11-c22) (-1+Qh3c) (-1+Qv3c)^2)/(2 (c22-c22 (Qh3c^2+2 Qv3c+Qh3c Qv3c (-1+(-2+Qv3c) Qv3c))+c11 (-1+Qh3c^2-2 Qh3c Qv3c+Qv3c (2+(-1+Qv3c)^2 Qv3c))))/.{Qv3c->Q13,Qh3c->Q23};
S23=((c11-c22) (-1+Qh3c)^2 (-1+Qv3c))/(2 (c22-c22 (Qh3c (2+(-1+Qh3c)^2 Qh3c)-2 Qh3c Qv3c+Qv3c^2)+c11 (-1+2 Qh3c+Qh3c (-1+(-2+Qh3c) Qh3c) Qv3c+Qv3c^2)))/.{Qv3c->Q13,Qh3c->Q23};
Hortho =t0^2 + x^2/(Q32 vnmox^2) + y^2/(Q31 vnmoy^2);
S1h = (S13 y^2/(Q31 vnmoy^2) + S12 t0^2)/(y^2/(Q31 vnmoy^2) + t0^2);
S2h = (S23 x^2/(Q32 vnmox^2) + S21 t0^2)/(x^2/(Q32 vnmox^2) + t0^2);
S3h = (S32 x^2/(Q32 vnmox^2) + S31 y^2/(Q31 vnmoy^2))/(x^2/(Q32 vnmox^2) + y^2/(Q31 vnmoy^2));
Shat =   (S1h x^2/(Q32 vnmox^2)  + S2h y^2 /(Q31 vnmoy^2) + S3h t0^2)/(x^2/(Q32 vnmox^2)  + y^2 /(Q31 vnmoy^2) + t0^2);

Q1h = (Q21 y^2/(Q31 vnmoy^2) + Q31 t0^2)/(y^2/(Q31 vnmoy^2) + t0^2);
Q2h = (Q12 x^2/(Q32 vnmox^2) + Q32 t0^2)/(x^2/(Q32 vnmox^2) + t0^2);
Q3h = (Q13 x^2/(Q32 vnmox^2) + Q23 y^2/(Q31 vnmoy^2))/(x^2/(Q32 vnmox^2) + y^2/(Q31 vnmoy^2));

zone = Hortho(1-Shat) + Shat Sqrt[Hortho^2 +2/Shat(((Q1h-1)t0^2 y^2 )/(Q31 vnmoy^2) +((Q2h-1)t0^2 x^2)/(Q32 vnmox^2) +((Q3h-1)x^2 y^2)/((Q32 vnmox^2)(Q31 vnmoy^2)))]/.{x->xoffset,y->yoffset}/.{vnmox ->Sqrt[c33(1+2δ2)],vnmoxy->Infinity,vnmoy->Sqrt[c33(1+2δ1)]}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
zonerms[radius_,angle_]:= 100Abs[Sqrt[zone]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};

onelayerzone= Rasterize[ParametricPlot3D[{xoffset,yoffset,100Abs[Sqrt[zone]/exacttime-1]}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]},{R,0,1},{theta,0,2Pi},PlotRange->{All,All,All},Boxed->False,BoxRatios->{1,1,0.3},Axes->{True,True,False},AxesLabel->{"x/H","y/H"},Ticks->{{-3,-2,-1,0,1,2,3},{-3,-2,-1,0,1,2,3}},ColorFunctionScaling->False,ColorFunction->ColorData[{"DarkRainbow",{0,2.5}}],Lighting->{{"Ambient",White,{0,0,5}}},ViewPoint->{0,0,Infinity},LabelStyle->{Black,FontSize->22,FontFamily->"Times",FontWeight->"Bold"},TicksStyle->Directive[Thick,20],Mesh->19,ImageSize->500,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["Error %",Above],LabelStyle->Directive[FontSize->22,FontFamily->"Times",Black,Bold]]]]
(*Export["~/Desktop/onelayerzone.eps",%]*)

Export["junk_ma.eps",onelayerzone,"EPS"]
