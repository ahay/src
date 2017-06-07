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
onelayernmo=Rasterize[ParametricPlot3D[{xoffset,yoffset,100Abs[Sqrt[orthoell]/exacttime-1]}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]},{R,0,1},{theta,0,2Pi},PlotRange->{All,All,All},Boxed->False,BoxRatios->{1,1,0.3},Axes->{True,True,False},AxesLabel->{"x/H","y/H"},Ticks->{{-3,-2,-1,0,1,2,3},{-3,-2,-1,0,1,2,3}},ColorFunction->ColorData["DarkRainbow"],Lighting->{{"Ambient",White,{0,0,5}}},ViewPoint->{0,0,Infinity},LabelStyle->{Black,FontSize->22,FontFamily->"Times",FontWeight->"Bold"},TicksStyle->Directive[Thick,20],Mesh->19,ImageSize->500,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["Error %",Above],LabelStyle->Directive[FontSize->22,FontFamily->"Times",Black,Bold]]]]
(*Export["~/Desktop/onelayernmo.eps",%]*)

Export["junk_ma.eps",onelayernmo,"EPS"]
