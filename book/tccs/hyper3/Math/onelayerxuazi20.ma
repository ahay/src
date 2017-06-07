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

Al-Dajaini (1998)
η2 = (c11(c33-c55))/(2c13(c13+2c55)+2c33 c55)-1/2;
η1 = (c22(c33-c44))/(2c23(c23+2c44)+2c33 c44)-1/2;
η3  =(c22(c11-c66))/(2c12(c12+2c66)+2c11 c66)-1/2; 
ηxy = Sqrt[(1+2η1)(1+2η2)/(1+2η3)]-1 ;

Q4 =  ((-2η2) (1/(c33(1+2δ2)))^2c^4+  (-2ηxy) (1/(c33(1+2δ2))) (1/(c33(1+2δ1)))c^2 s^2 +  (-2η1 ) (1/(c33(1+2δ1)))^2s^4)/t0^2;
vhor = (1/2)((c11+c66)c^2+(c22+c66)s^2)+(1/2)Sqrt[((c11-c66)c^2-(c22-c66)s^2)^2+4(c12+c66)^2c^2s^2];(*phase velocity suggested by the authors*) 
Q=Q4/(1/vhor - Q2);
al = t0^2 + Q2 rortho^2 + Q4 rortho^4/(1+Q rortho^2)/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
alrms[radius_,angle_]:= 100Abs[Sqrt[al]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};

Xu et al. (2005)
η2 = (c11(c33-c55))/(2c13(c13+2c55)+2c33 c55)-1/2;
η1 = (c22(c33-c44))/(2c23(c23+2c44)+2c33 c44)-1/2;
η3  =(c22(c11-c66))/(2c12(c12+2c66)+2c11 c66)-1/2; 
ηxy = Sqrt[(1+2η1)(1+2η2)/(1+2η3)]-1 ;
ηav =  (η2)c^2-  (η3)c^2 s^2 +  (η1 )s^2;
xu = t0^2 + Q2 rortho^2 - 2Q2 ηav rortho^4/(t0^2/Q2+(1+2ηav) rortho^2)/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
xurms[radius_,angle_]:= 100Abs[Sqrt[xu]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};


If the azimuth kowledge is off by 20 %
cwr = cos Cos[1.2azi]+ sin Sin[1.2azi];
swr = sin Cos[1.2azi]-cos Sin[1.2azi];
ηavwr =  (η2)cwr^2-  (η3)cwr^2 swr^2 +  (η1 )swr^2;
Q2wr = (1/(c33(1+2δ2)))cwr^2 +(1/(c33(1+2δ1)))swr^2;
xuwr = t0^2 + Q2wr rortho^2 - 2Q2wr ηavwr rortho^4/(t0^2/Q2wr+(1+2ηavwr) rortho^2)/.t0->2 depth/Sqrt[c33]/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
xuwrrms[radius_,angle_]:= 100Abs[Sqrt[xuwr]/exacttime-1]/.{theta->angle}/.{p1->maxr Cos[angle],p2->maxr Sin[angle],R->radius};
onelayerxuazi20= Rasterize[ParametricPlot3D[{xoffset,yoffset,100Abs[Sqrt[xuwr]/exacttime-1]}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]},{R,0,1},{theta,azi,2Pi+azi},PlotRange->{All,All,All},Boxed->False,BoxRatios->{1,1,0.3},Axes->{True,True,False},AxesLabel->{"x/H","y/H"},ColorFunctionScaling->False,ColorFunction->ColorData[{"DarkRainbow",{0,2.5}}],Ticks->{{-3,-2,-1,0,1,2,3},{-4,-3,-2,-1,0,1,2,3,4}},Lighting->{{"Ambient",White,{0,0,5}}},ViewPoint->{0,0,Infinity},LabelStyle->{Black,FontSize->22,FontFamily->"Times",Bold},TicksStyle->Directive[Thick,20],Mesh->19,ImageSize->500,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["Error %",Above],LabelStyle->Directive[FontSize->22,FontFamily->"Times",Black,Bold]]]]
(*Export["~/Desktop/onelayerxuazi20.eps",%]*)


Export["junk_ma.eps",onelayerxuazi20,"EPS"]
