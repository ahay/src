Homogeneous Orthorhombic Layer
ClearAll[c11,c12,c13,c12,c23,c22,c33,c44,c55,c66,n1,n2,n3,two,trace,G,H,angle,vhelbig,vpphase,vecdiff,vpgroup,vpg,pn12,pn22,pn32]

Exact velocity
p1=.;p2=.;p3=.;
azi = 30 Pi/180;
mat = {{c11,c12,c13,0,0,0},{c12,c22,c23,0,0,0},{c13,c23,c33,0,0,0},{0,0,0,c44,0,0},{0,0,0,0,c55,0},{0,0,0,0,0,c66}};
matn = {{n1,0,0,0,n3,n2},{0,n2,0,n3,0,n1},{0,0,n3,n2,n1,0}} ;
bond = {{Cos[azi]^2,Sin[azi]^2,0,0,0,-Sin[2azi]},{Sin[azi]^2,Cos[azi]^2,0,0,0,Sin[2azi]},{0,0,1,0,0,0},{0,0,0,Cos[azi],Sin[azi],0},{0,0,0,-Sin[azi],Cos[azi],0},{Sin[2azi]/2,-Sin[2azi]/2,0,0,0,Cos[2azi]}};
matb= bond.mat.Transpose[bond];
chris = matn.matb.Transpose[matn];
chrisp = chris /.{n1->p1 ,n2->p2 ,n3->p3};
chrisp//MatrixForm
det = Det[chrisp-IdentityMatrix[3]];
(p2 ((1/4 Sqrt[3] ((3 c11)/4+c12/4)-1/4 Sqrt[3] ((3 c12)/4+c22/4)-(Sqrt[3] c66)/4) p1+(1/4 Sqrt[3] ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)-1/4 Sqrt[3] ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)+c66/4) p2)+p1 ((3/4 ((3 c11)/4+c12/4)+1/4 ((3 c12)/4+c22/4)+(3 c66)/4) p1+(3/4 ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)+1/4 ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)-(Sqrt[3] c66)/4) p2)+(c44/4+(3 c55)/4) p3^2	p1 ((1/4 Sqrt[3] ((3 c11)/4+c12/4)-1/4 Sqrt[3] ((3 c12)/4+c22/4)-(Sqrt[3] c66)/4) p1+(1/4 Sqrt[3] ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)-1/4 Sqrt[3] ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)+c66/4) p2)+p2 ((1/4 ((3 c11)/4+c12/4)+3/4 ((3 c12)/4+c22/4)-(3 c66)/4) p1+(1/4 ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)+3/4 ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)+(Sqrt[3] c66)/4) p2)+(-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p3^2	(c44/4+(3 c55)/4) p1 p3+(-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p2 p3+(((3 c13)/4+c23/4) p1+((Sqrt[3] c13)/4-(Sqrt[3] c23)/4) p2) p3
p1 ((3/4 ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)+1/4 ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)-(Sqrt[3] c66)/4) p1+(3/4 (c11/4+(3 c12)/4)+1/4 (c12/4+(3 c22)/4)-(3 c66)/4) p2)+p2 ((1/4 Sqrt[3] ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)-1/4 Sqrt[3] ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)+c66/4) p1+(1/4 Sqrt[3] (c11/4+(3 c12)/4)-1/4 Sqrt[3] (c12/4+(3 c22)/4)+(Sqrt[3] c66)/4) p2)+(-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p3^2	p2 ((1/4 ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)+3/4 ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)+(Sqrt[3] c66)/4) p1+(1/4 (c11/4+(3 c12)/4)+3/4 (c12/4+(3 c22)/4)+(3 c66)/4) p2)+p1 ((1/4 Sqrt[3] ((Sqrt[3] c11)/4-(Sqrt[3] c12)/4)-1/4 Sqrt[3] ((Sqrt[3] c12)/4-(Sqrt[3] c22)/4)+c66/4) p1+(1/4 Sqrt[3] (c11/4+(3 c12)/4)-1/4 Sqrt[3] (c12/4+(3 c22)/4)+(Sqrt[3] c66)/4) p2)+((3 c44)/4+c55/4) p3^2	(-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p1 p3+((3 c44)/4+c55/4) p2 p3+(((Sqrt[3] c13)/4-(Sqrt[3] c23)/4) p1+(c13/4+(3 c23)/4) p2) p3
((3 c13)/4+c23/4) p1 p3+((Sqrt[3] c13)/4-(Sqrt[3] c23)/4) p2 p3+((c44/4+(3 c55)/4) p1+(-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p2) p3	((Sqrt[3] c13)/4-(Sqrt[3] c23)/4) p1 p3+(c13/4+(3 c23)/4) p2 p3+((-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p1+((3 c44)/4+c55/4) p2) p3	p2 ((-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p1+((3 c44)/4+c55/4) p2)+p1 ((c44/4+(3 c55)/4) p1+(-((Sqrt[3] c44)/4)+(Sqrt[3] c55)/4) p2)+c33 p3^2

)
Hamiltonian to derive coefficients
F = det;(*Hamiltonian*)
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
maxx = 0.85/(Sqrt[c11])/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
maxy = 0.85/(Sqrt[c22])/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
maxr = R maxx maxy/Sqrt[(maxy^2 Cos[theta-azi]^2 + maxx^2 Sin[theta-azi]^2)]; (*Ellipse in polar*)
xoffset = (2 dx)/.{Δd->depth}/.{p1->maxr Cos[theta-azi],p2->maxr Sin[theta-azi]};//Simplify
yoffset = (2dy)/.{Δd->depth}/.{p1->maxr Cos[theta-azi],p2->maxr Sin[theta-azi]};//Simplify
exacttime = (2dt)/.{Δd->depth}/.{p1->maxr Cos[theta-azi],p2->maxr Sin[theta-azi]};
NMO Ellipse
rortho = Sqrt[xoffset^2+yoffset^2];
cos = xoffset/rortho;
sin = yoffset/rortho;

c = cos Cos[azi]+ sin Sin[azi];
s = sin Cos[azi]-cos Sin[azi];

δ2=((c13+c55)^2-(c33-c55)^2)/(2*c33*(c33-c55));
δ1=((c23+c44)^2-(c33-c44)^2)/(2*c33*(c33-c44));
Q2 = (1/(c33(1+2δ2)))c^2 +(1/(c33(1+2δ1)))s^2;

orthoell = t0^2 + Q2 rortho^2/.t0->2 depth/Sqrt[c33] /.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
nmorms[radius_,angle_]:= 100Abs[Sqrt[orthoell]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};

Map[SetOptions[#,Prolog->{{EdgeForm[],Texture[{{{0,0,0,0}}}],Polygon[#,VertexTextureCoordinates->#]&[{{0,0},{1,0},{1,1}}]}}]&,{Graphics3D,ContourPlot3D,ListContourPlot3D,ListPlot3D,Plot3D,ListSurfacePlot3D,ListVectorPlot3D,ParametricPlot3D,RegionPlot3D,RevolutionPlot3D,SphericalPlot3D,VectorPlot3D,BarChart3D}];

GMA
(*Need truemax for varying shooting rays  *)
factor = 0.85;
varmaxx = factor/Sqrt[matb[[1,1]]]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
varmaxy = factor/Sqrt[matb[[2,2]]]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};

extfunc = t0^2 + B1 x^2 + B2 x y + B3 y^2 + Sqrt[t0^4 + 2 t0^2 ( B1 x^2 + B2 x y + B3 y^2) + (D1 x^4 + D2 x^3 y + D3 x^2 y^2 + D4 x y^3 + D5 y^4) ];
exttt = t0^2  + x^2/vnmox^2 + y^2/vnmoy^2 + x y/vnmoxy^2 + (A1 x^4 + A2 x^3 y + A3 x^2 y^2 + A4 x y^3 + A5 y^4)/extfunc;

B1form =(t0^2 (-px T vnmox^2+x))/(vnmox^2 x (-T^2+t0^2+px T x))+(A1 vnmox^2 x^2)/((T-t0) (T+t0) vnmox^2-x^2);
D1form = (px^2 T^2 t0^4)/((T^2-t0^2)^2 x^2)+(2 t0^4 (-px T^3+px T t0^2+px^3 T^3 vnmox^2))/((T^2-t0^2)^3 vnmox^2 x)+(t0^4 (-T^2+t0^2+px^2 T^2 vnmox^2)^2)/((T^2-t0^2)^2 vnmox^4 (-T^2+t0^2+px T x)^2)-(2 px^2 T^2 t0^4 (-T^2+t0^2+px^2 T^2 vnmox^2))/((T^2-t0^2)^3 vnmox^2 (-T^2+t0^2+px T x))+(2 A1 t0^2 vnmox^2)/(-T^2 vnmox^2+t0^2 vnmox^2+x^2);
B3form = (t0^2 (-py T vnmoy^2+y))/(vnmoy^2 y (-T^2+t0^2+py T y))+(A5 vnmoy^2 y^2)/((T-t0) (T+t0) vnmoy^2-y^2);
D5form = (py^2 T^2 t0^4)/((T^2-t0^2)^2 y^2)+(2 t0^4 (-py T^3+py T t0^2+py^3 T^3 vnmoy^2))/((T^2-t0^2)^3 vnmoy^2 y)+(t0^4 (-T^2+t0^2+py^2 T^2 vnmoy^2)^2)/((T^2-t0^2)^2 vnmoy^4 (-T^2+t0^2+py T y)^2)-(2 py^2 T^2 t0^4 (-T^2+t0^2+py^2 T^2 vnmoy^2))/((T^2-t0^2)^3 vnmoy^2 (-T^2+t0^2+py T y))+(2 A5 t0^2 vnmoy^2)/(-T^2 vnmoy^2+t0^2 vnmoy^2+y^2);

vnmoxazi=1/Sqrt[0.18562680456768524`];
vnmoxyazi=1/Sqrt[0.04740254613545505`];
vnmoyazi=1/Sqrt[0.15825893179610648`];
A1gmaoneazi = -0.03118810360189571`*2t0^2;
A2gmaoneazi=-0.03729846606084341`*2t0^2;
A3gmaoneazi=-0.06106669925993767`*2 t0^2;
A4gmaoneazi = -0.013944911548059411` *2t0^2;
A5gmaoneazi=-0.0163954146735523`*2t0^2;

Recalculate the ray parameter for shooting along x and y axis
p1x= varmaxx
p2x= p2/.FindRoot[0==Abs[Re[2dy/.{Δd->depth}/.{p1->p1x}]],{p2,0}][[1]]
{2dx/.{Δd->depth},2dy/.{Δd->depth}}/.{p1->p1x,p2->p2x}
p2y= varmaxy
p1y= p1/.FindRoot[0==Abs[Re[2dx/.{Δd->depth}/.{p2->p2y}]],{p1,0}][[1]]
{2dx/.{Δd->depth},2dy/.{Δd->depth}}/.{p1->p1y,p2->p2y}
0.288775
0.00390271
{4.01114 -1.0197*10^-15 I,-2.14322*10^-16-4.68193*10^-17 I}
0.28202
0.0324634
{4.03615*10^-14-2.59838*10^-16 I,4.02191 -2.06882*10^-15 I}
varxoffset = (2 dx)/.{Δd->depth}/.{p1->p1x,p2->p2x};//Simplify
varyoffset = (2dy)/.{Δd->depth}/.{p2->p2y,p1->p1y};//Simplify
varxtime = (2dt)/.{Δd->depth}/.{p1->p1x,p2->p2x};
varytime = (2dt)/.{Δd->depth}/.{p2->p2y,p1->p1y};
B1gmaoneazi = B1form/. {x->varxoffset,T->varxtime,px->p1x};
D1gmaoneazi= D1form/. {x->varxoffset,T->varxtime,px->p1x};
B3gmaoneazi= B3form/. {y->varyoffset,T->varytime,py->p2y};
D5gmaoneazi= D5form/. {y->varyoffset,T->varytime,py->p2y};
xoffset1 = (2 dx)/.{Δd->depth}/.{p1->p1r1,p2->p2r1}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}; (*For x=y*)
yoffset1 = (2 dy)/.{Δd->depth}/.{p1->p1r1,p2->p2r1}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
Txyf = (2 dt)/.{Δd->depth}/.{p1->p1r1,p2->p2r1}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
roffset1 = xoffset1-yoffset1;
xoffset2 =(2 dx)/.{Δd->depth}/.{p1->p1r2,p2->p2r2}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};(*For x=-y*)
yoffset2 = (2 dy)/.{Δd->depth}/.{p1->p1r2,p2->p2r2}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
Txynf =(2 dt)/.{Δd->depth}/.{p1->p1r2,p2->p2r2}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
roffset2 = xoffset2+yoffset2;
xoffset3 = (2 dx)/.{Δd->depth}/.{p1->p1r3,p2->p2r3}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}; (*For -x=y*)
yoffset3 = (2 dy)/.{Δd->depth}/.{p1->p1r3,p2->p2r3}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
Txnyf = (2 dt)/.{Δd->depth}/.{p1->p1r3,p2->p2r3}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
roffset3 = xoffset3+yoffset3;
xoffset4 =(2 dx)/.{Δd->depth}/.{p1->p1r4,p2->p2r4}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};(*For -x=-y*)
yoffset4 = (2 dy)/.{Δd->depth}/.{p1->p1r4,p2->p2r4}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
Txnynf =(2 dt)/.{Δd->depth}/.{p1->p1r4,p2->p2r4}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
roffset4 = xoffset4-yoffset4;
p1r1f= 0.2
p2r1f= p2r1/.FindRoot[0==Abs[Re[roffset1/.{p1r1->p1r1f}]],{p2r1,0.2},AccuracyGoal->5][[1]]
{xoffset1,yoffset1}/.{p1r1->p1r1f,p2r1->p2r1f}

p1r2f= 0.2
p2r2f= p2r2/.FindRoot[0==Abs[Re[roffset2/.{p1r2->p1r2f}]],{p2r2,-0.2},AccuracyGoal->5][[1]]
{xoffset2,yoffset2}/.{p1r2->p1r2f,p2r2->p2r2f}

eq1gmaoneazi = ( 2T py - D[exttt,y]/.y->0)/.{x->varxoffset,T->varxtime,px->p1x,py->p2x}/.{B1->B1gmaoneazi,D1->D1gmaoneazi,B3->B3gmaoneazi,D5->D5gmaoneazi}/.{A1->A1gmaoneazi,A2->A2gmaoneazi,A3->A3gmaoneazi,A4->A4gmaoneazi,A5->A5gmaoneazi}/.{vnmox ->vnmoxazi,vnmoxy->vnmoxyazi,vnmoy->vnmoyazi}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}//N;
eq2gmaoneazi = ( 2T px - D[exttt,x]/.x->0)/. {y->varyoffset,T->varytime,py->p2y,px->p1y}/.{B1->B1gmaoneazi,D1->D1gmaoneazi,B3->B3gmaoneazi,D5->D5gmaoneazi}/.{A1->A1gmaoneazi,A2->A2gmaoneazi,A3->A3gmaoneazi,A4->A4gmaoneazi,A5->A5gmaoneazi}/.{vnmox ->vnmoxazi,vnmoxy->vnmoxyazi,vnmoy->vnmoyazi}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}//N;
eq3gmaoneazi = (Txy^2 - exttt/.{x->y})/.Txy->Txyf/.{y->yoffset1}/.{B1->B1gmaoneazi,D1->D1gmaoneazi,B3->B3gmaoneazi,D5->D5gmaoneazi}/.{A1->A1gmaoneazi,A2->A2gmaoneazi,A3->A3gmaoneazi,A4->A4gmaoneazi,A5->A5gmaoneazi}/.{vnmox ->vnmoxazi,vnmoxy->vnmoxyazi,vnmoy->vnmoyazi}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}/.{p1r1->p1r1f,p2r1->p2r1f}//N;
eq4gmaoneazi = (Txyn^2- exttt/.{x->-y})/.Txyn->Txynf/.{y->-yoffset2}/.{B1->B1gmaoneazi,D1->D1gmaoneazi,B3->B3gmaoneazi,D5->D5gmaoneazi}/.{A1->A1gmaoneazi,A2->A2gmaoneazi,A3->A3gmaoneazi,A4->A4gmaoneazi,A5->A5gmaoneazi}/.{vnmox ->vnmoxazi,vnmoxy->vnmoxyazi,vnmoy->vnmoyazi}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25}/.{p1r2->p1r2f,p2r2->p2r2f}//N;

{B2gmaoneazi,D2gmaoneazi,D3gmaoneazi,D4gmaoneazi}={B2,D2,D3,D4}/.FindRoot[{Re[eq1gmaoneazi]==0,Re[eq2gmaoneazi]==0,Re[eq3gmaoneazi]==0,Re[eq4gmaoneazi]==0},{{B2,0.015},{D2,0.02},{D3,0.01},{D4,-0.01}},MaxIterations->200,AccuracyGoal->5]
0.2
0.205865
{2.94045 -3.55106*10^-16 I,2.94045 -3.49486*10^-16 I}
0.2
-0.163014
{2.29382 -5.7724*10^-16 I,-2.29382+4.50231*10^-16 I}
{0.177078,0.0483517,0.0775516,-0.0191342}

extttortho = exttt/.{x->xoffset,y->yoffset}/.{B1->B1gmaoneazi,D1->D1gmaoneazi,B3->B3gmaoneazi,D5->D5gmaoneazi}/.{A1->A1gmaoneazi,A2->A2gmaoneazi,A3->A3gmaoneazi,A4->A4gmaoneazi,A5->A5gmaoneazi}/.{B2->B2gmaoneazi,D2->D2gmaoneazi,D3->D3gmaoneazi,D4->D4gmaoneazi}/.{vnmox ->vnmoxazi,vnmoxy->vnmoxyazi,vnmoy->vnmoyazi}/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
gmarms[radius_,angle_]:= 100Abs[Sqrt[extttortho]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};
onelayergmaazi= Rasterize[ParametricPlot3D[{xoffset,yoffset,100Abs[Sqrt[extttortho]/exacttime-1]}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]},{R,0,1},{theta,azi,2Pi+azi},PlotPoints->120,PlotRange->{All,All,All},Boxed->False,BoxRatios->{1,1,0.3},Axes->{True,True,False},AxesLabel->{"x/H","y/H"},Ticks->{{-3,-2,-1,0,1,2,3},{-4,-3,-2,-1,0,1,2,3,4}},ColorFunctionScaling->False,ColorFunction->ColorData[{"DarkRainbow",{0,2.5}}],Lighting->{{"Ambient",White,{0,0,5}}},ViewPoint->{0,0,Infinity},LabelStyle->{Black,FontSize->22,FontFamily->"Times",Bold},TicksStyle->Directive[Thick,20],Mesh->19,ImageSize->500,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["Error %",Above],LabelStyle->Directive[FontSize->22,FontFamily->"Times",Black,Bold]]]]
(*Export["~/Desktop/onelayergmaazi.eps",%]*)

Export["junk_ma.eps",onelayergmaazi,"EPS"]
