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

GMA
extfunc = t0^2 + B1 x^2 + B2 x y + B3 y^2 + Sqrt[t0^4 + 2 t0^2 ( B1 x^2 + B2 x y + B3 y^2) + (D1 x^4 + D2 x^3 y + D3 x^2 y^2 + D4 x y^3 + D5 y^4) ];
exttt = t0^2  + x^2/vnmox^2 + y^2/vnmoy^2 + x y/vnmoxy^2 + (A1 x^4 + A2 x^3 y + A3 x^2 y^2 + A4 x y^3 + A5 y^4)/extfunc;

B1form =(t0^2 (-px T vnmox^2+x))/(vnmox^2 x (-T^2+t0^2+px T x))+(A1 vnmox^2 x^2)/((T-t0) (T+t0) vnmox^2-x^2);
D1form = (px^2 T^2 t0^4)/((T^2-t0^2)^2 x^2)+(2 t0^4 (-px T^3+px T t0^2+px^3 T^3 vnmox^2))/((T^2-t0^2)^3 vnmox^2 x)+(t0^4 (-T^2+t0^2+px^2 T^2 vnmox^2)^2)/((T^2-t0^2)^2 vnmox^4 (-T^2+t0^2+px T x)^2)-(2 px^2 T^2 t0^4 (-T^2+t0^2+px^2 T^2 vnmox^2))/((T^2-t0^2)^3 vnmox^2 (-T^2+t0^2+px T x))+(2 A1 t0^2 vnmox^2)/(-T^2 vnmox^2+t0^2 vnmox^2+x^2);
B3form = (t0^2 (-py T vnmoy^2+y))/(vnmoy^2 y (-T^2+t0^2+py T y))+(A5 vnmoy^2 y^2)/((T-t0) (T+t0) vnmoy^2-y^2);
D5form = (py^2 T^2 t0^4)/((T^2-t0^2)^2 y^2)+(2 t0^4 (-py T^3+py T t0^2+py^3 T^3 vnmoy^2))/((T^2-t0^2)^3 vnmoy^2 y)+(t0^4 (-T^2+t0^2+py^2 T^2 vnmoy^2)^2)/((T^2-t0^2)^2 vnmoy^4 (-T^2+t0^2+py T y)^2)-(2 py^2 T^2 t0^4 (-T^2+t0^2+py^2 T^2 vnmoy^2))/((T^2-t0^2)^3 vnmoy^2 (-T^2+t0^2+py T y))+(2 A5 t0^2 vnmoy^2)/(-T^2 vnmoy^2+t0^2 vnmoy^2+y^2);
(*Need truemax for varying shooting rays  *)
factor = 0.85;
varmaxx = factor/Sqrt[c11]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
varmaxy = factor/Sqrt[c22]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
varxoffset = (2 dx)/.{Δd->depth}/.{p1->varmaxx,p2->0};//Simplify
varyoffset = (2dy)/.{Δd->depth}/.{p2->varmaxy,p1->0};//Simplify
varxtime = (2dt)/.{Δd->depth}/.{p1->varmaxx,p2->0};
varytime = (2dt)/.{Δd->depth}/.{p2->varmaxy,p1->0};
A1gmaone = -0.043640640935251676*2t0^2;
A2gmaone=0;
A3gmaone=-0.030729542044832935*2 t0^2;
A4gmaone = 0;
A5gmaone=-0.014055263078564768`*2t0^2;
B1gmaone = B1form/. {x->varxoffset,T->varxtime,px->varmaxx};
D1gmaone= D1form/. {x->varxoffset,T->varxtime,px->varmaxx};
B3gmaone= B3form/. {y->varyoffset,T->varytime,py->varmaxy};
D5gmaone= D5form/. {y->varyoffset,T->varytime,py->varmaxy};

(*For x=y*)
xoffset1 = (2 dx)/.{Δd->depth}/.{p1->p1r1,p2->p2r1}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}; 
yoffset1 = (2 dy)/.{Δd->depth}/.{p1->p1r1,p2->p2r1}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
Txyf = (2 dt)/.{Δd->depth}/.{p1->p1r1,p2->p2r1}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
roffset1 = xoffset1-yoffset1;

(*For x=-y*)
xoffset2 =(2 dx)/.{Δd->depth}/.{p1->p1r2,p2->p2r2}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}/.{p1r2->p1r1,p2r2->p2r1};
yoffset2 = (2 dy)/.{Δd->depth}/.{p1->p1r2,p2->p2r2}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}/.{p1r2->p1r1,p2r2->p2r1};
Tyxf =(2 dt)/.{Δd->depth}/.{p1->p1r2,p2->p2r2}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}/.{p1r2->p1r1,p2r2->p2r1};
roffset2 = xoffset2-yoffset2; (*Actually xoffset2 + yoffset2 but this is equal to the x=y case in this example*)

p1r1f= 0.2
p2r1f= p2r1/.FindRoot[0==Abs[Re[roffset1/.{p1r1->p1r1f}]],{p2r1,0.2}][[1]]
{xoffset1,yoffset1}/.{p1r1->p1r1f,p2r1->p2r1f}
eq1gmaone = ( 2T py - D[exttt,y]/.y->0)/.{x->varxoffset,T->varxtime,px->varmaxx,py->0}/.{B1->B1gmaone,D1->D1gmaone,B3->B3gmaone,D5->D5gmaone}/.{A1->A1gmaone,A2->A2gmaone,A3->A3gmaone,A4->A4gmaone,A5->A5gmaone}/.{vnmox ->Sqrt[c33(1+2δ2)],vnmoxy->Infinity,vnmoy->Sqrt[c33(1+2δ1)]}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}//N;
eq2gmaone = ( 2T px - D[exttt,x]/.x->0)/. {y->varyoffset,T->varytime,py->varmaxy,px->0}/.{B1->B1gmaone,D1->D1gmaone,B3->B3gmaone,D5->D5gmaone}/.{A1->A1gmaone,A2->A2gmaone,A3->A3gmaone,A4->A4gmaone,A5->A5gmaone}/.{vnmox ->Sqrt[c33(1+2δ2)],vnmoxy->Infinity,vnmoy->Sqrt[c33(1+2δ1)]}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}//N;
eq3gmaone = (Txy^2 - exttt/.{x->y})/.Txy->Txyf/.{y->-xoffset1}/.{B1->B1gmaone,D1->D1gmaone,B3->B3gmaone,D5->D5gmaone}/.{A1->A1gmaone,A2->A2gmaone,A3->A3gmaone,A4->A4gmaone,A5->A5gmaone}/.{vnmox ->Sqrt[c33(1+2δ2)],vnmoxy->Infinity,vnmoy->Sqrt[c33(1+2δ1)]}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}/.{p1r1->p1r1f,p2r1->p2r1f}//N;
eq4gmaone = (Tyx^2- exttt/.{x->-y})/.Tyx->Tyxf/.{y->-xoffset2}/.{B1->B1gmaone,D1->D1gmaone,B3->B3gmaone,D5->D5gmaone}/.{A1->A1gmaone,A2->A2gmaone,A3->A3gmaone,A4->A4gmaone,A5->A5gmaone}/.{vnmox ->Sqrt[c33(1+2δ2)],vnmoxy->Infinity,vnmoy->Sqrt[c33(1+2δ1)]}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}/.{p1r1->p1r1f,p2r1->p2r1f}//N;
{B2gmaone,D2gmaone,D3gmaone,D4gmaone}={B2,D2,D3,D4}/.FindRoot[{Re[eq1gmaone]==0,Re[eq2gmaone]==0,Re[eq3gmaone]==0,Re[eq4gmaone]==0},{{B2,0},{D2,0},{D3,0.001},{D4,0}},MaxIterations->200,AccuracyGoal->5]
0.2
0.168905
{1.95451 -2.1679*10^-16 I,1.95451 -1.63115*10^-16 I}
{0.,0.,-0.0749452,0.}

extttortho = exttt/.{x->xoffset,y->yoffset}/.{B1->B1gmaone,D1->D1gmaone,B3->B3gmaone,D5->D5gmaone}/.{A1->A1gmaone,A2->A2gmaone,A3->A3gmaone,A4->A4gmaone,A5->A5gmaone}/.{B2->B2gmaone,D2->D2gmaone,D3->D3gmaone,D4->D4gmaone}/.{vnmox ->Sqrt[c33(1+2δ2)],vnmoxy->Infinity,vnmoy->Sqrt[c33(1+2δ1)]}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
gmarms[radius_,angle_]:= 100Abs[Sqrt[extttortho]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};

onelayergma= Rasterize[ParametricPlot3D[{xoffset,yoffset,100Abs[Sqrt[extttortho]/exacttime-1]}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]},{R,0,1},{theta,0,2Pi},PlotPoints->100,PlotRange->{All,All,All},Boxed->False,BoxRatios->{1,1,0.3},Axes->{True,True,False},AxesLabel->{"x/H","y/H"},Ticks->{{-3,-2,-1,0,1,2,3},{-3,-2,-1,0,1,2,3}},ColorFunctionScaling->False,ColorFunction->ColorData[{"DarkRainbow",{0,2.5}}],Lighting->{{"Ambient",White,{0,0,5}}},ViewPoint->{0,0,Infinity},LabelStyle->{Black,FontSize->22,FontFamily->"Times",FontWeight->"Bold"},TicksStyle->Directive[Thick,20],Mesh->19,ImageSize->500,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["Error %",Above],LabelStyle->Directive[FontSize->22,FontFamily->"Times",Black,Bold]]]]
(*Export["~/Desktop/onelayergma.eps",%]*)

Export["junk_ma.eps",onelayergma,"EPS"]
