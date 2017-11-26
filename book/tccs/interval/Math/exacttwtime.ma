%Exact velocity
mat = {{c11,c12,c13,0,0,0},{c12,c22,c23,0,0,0},{c13,c23,c33,0,0,0},{0,0,0,c44,0,0},{0,0,0,0,c55,0},{0,0,0,0,0,c66}};
matn = {{n1,0,0,0,n3,n2},{0,n2,0,n3,0,n1},{0,0,n3,n2,n1,0}} ;
chris = matn.mat.Transpose[matn];
chrisp = chris /.{n1->p1,n2->p2,n3->p3};
det = Det[chrisp-IdentityMatrix[3]];
%Hamiltonian to derive coefficients
F = det;(*Hamiltonian*)

p3sqf=p3sq/.Solve[0==F/.p3->Sqrt[p3sq],p3sq][[1]]; (*Orthorhombic qP*)

p3sqf1=p3sqf/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
p3sqf2 = p3sqf/.{c11->11.7,c22->13.5,c33->9,c44->1.728,c55->1.44,c66->2.246,c12->8.824,c23->5.981,c13->5.159};
p3sqf3=p3sqf/.{c11->12.6,c22->13.94,c33->8.9125,c44->2.5,c55->2,c66->2.18182,c12->2.7,c23->3.425,c13->3.15};

%DON'T FORGET A FACTOR OF  ΔZ  MULTIPLYING THE DERIVATIVES
dxdσ = D[F,p1]//Simplify;
dydσ = D[F,p2]//Simplify;
dzdσ = D[F,p3]//Simplify;

depth1 = 0.25; depth2 = 0.45;depth3=0.3; (*depth in km*)
depth = depth1 + depth2+depth3;
dx1 =Δd1 dxdσ /dzdσ /.{p3->Sqrt[p3sqf1]}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
dx2 =Δd2 dxdσ /dzdσ /.{p3->Sqrt[p3sqf2]}/.{c11->11.7,c22->13.5,c33->9,c44->1.728,c55->1.44,c66->2.246,c12->8.824,c23->5.981,c13->5.159};
dx3 =Δd3 dxdσ /dzdσ /.{p3->Sqrt[p3sqf3]}/.{c11->12.6,c22->13.94,c33->8.9125,c44->2.5,c55->2,c66->2.18182,c12->2.7,c23->3.425,c13->3.15};
dy1 = Δd1 dydσ /dzdσ /.{p3->Sqrt[p3sqf1]}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
dy2 = Δd2 dydσ /dzdσ /.{p3->Sqrt[p3sqf2]}/.{c11->11.7,c22->13.5,c33->9,c44->1.728,c55->1.44,c66->2.246,c12->8.824,c23->5.981,c13->5.159};
dy3 = Δd3 dydσ /dzdσ /.{p3->Sqrt[p3sqf3]}/.{c11->12.6,c22->13.94,c33->8.9125,c44->2.5,c55->2,c66->2.18182,c12->2.7,c23->3.425,c13->3.15};
dt1 =Δd1 ((p1 dxdσ+p2 dydσ)/dzdσ +p3)/.{p3->Sqrt[p3sqf1]}/.{c11->9,c22->9.84,c33->5.9375,c44->2,c55->1.6,c66->2.182,c12->3.6,c23->2.4,c13->2.25};
dt2 =Δd2 ((p1 dxdσ+p2 dydσ)/dzdσ +p3)/.{p3->Sqrt[p3sqf2]}/.{c11->11.7,c22->13.5,c33->9,c44->1.728,c55->1.44,c66->2.246,c12->8.824,c23->5.981,c13->5.159};
dt3 =Δd3 ((p1 dxdσ+p2 dydσ)/dzdσ +p3)/.{p3->Sqrt[p3sqf3]}/.{c11->12.6,c22->13.94,c33->8.9125,c44->2.5,c55->2,c66->2.18182,c12->2.7,c23->3.425,c13->3.15};
dx = (dx1 + dx2 + dx3)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3};
dy = (dy1 + dy2 + dy3)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3};
dt = (dt1 + dt2 + dt3)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3};
vert=({2dt1,2dt2,2dt3})/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->0,p2->0};
vertt= (2dt)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->0,p2->0};
vertn= vert/vertt;

vyav = Sqrt[vertn.{9.84,13.5,13.94}];
vxav = Sqrt[vertn.{9,11.7,12.6}];

maxx = 0.56/vxav ;
maxy = 0.56/vyav;
maxr = R maxx maxy/Sqrt[(maxy^2 Cos[theta]^2 + maxx^2 Sin[theta]^2)]; (*Ellipse in polar*)
xoffset =(2 dx)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};//Simplify
yoffset = (2dy)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};//Simplify
exacttime = (2dt)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};

twtime = Rasterize[ParametricPlot3D[{xoffset/depth,yoffset/depth,exacttime^2},{R,0,1},{theta,0,Pi/2},PlotRange->{All,All,All},BoxRatios->{1,1,0.3},
AxesLabel->{"x1(km)","x2(km)","4 t^2 (s^2)              "},ColorFunction->ColorData[{"DarkRainbow",{0.5,0.65}}],ColorFunctionScaling->False,LabelStyle->{Black,FontWeight->"Bold",FontSize->24},TicksStyle->Directive[Thick,20],Mesh->9,ImageSize->700]]

Export["junk_ma.eps",twtime,"EPS"]
junk_ma.eps
