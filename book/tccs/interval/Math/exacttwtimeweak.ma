%Exact velocity
mat = {{c11,c12,c13,0,0,0},{c12,c22,c23,0,0,0},{c13,c23,c33,0,0,0},{0,0,0,c44,0,0},{0,0,0,0,c55,0},{0,0,0,0,0,c66}};
matn = {{n1,0,0,0,n3,n2},{0,n2,0,n3,0,n1},{0,0,n3,n2,n1,0}} ;
chris = matn.mat.Transpose[matn];
chrisp = chris /.{n1->p1,n2->p2,n3->p3};
det = Det[chrisp-IdentityMatrix[3]];
%Hamiltonian to derive coefficients
F = det;(*Hamiltonian*)

p3sqf=p3sq/.Solve[0==F/.p3->Sqrt[p3sq],p3sq][[1]]; (*Orthorhombic qP*)
p3sqf1=p3sqf/.{c11->3.512016`,c22->3.875328`,c33->3.0276`,c44->0.1669390243902439`,c55->0.1521`,c66->0.27378`,c12->3.298434960000476`,c23->2.9819605693294067`,c13->2.8709922298203208`};
p3sqf2 = p3sqf/.{c11->4.51632`,c22->4.817408`,c33->3.7636`,c44->0.64896`,c55->0.6084`,c66->0.97344`,c12->2.995458418291907`,c23->2.3894757760926923`,c13->2.6577570205977543`};
p3sqf3=p3sqf/.{c11->8.7725`,c22->9.68`,c33->7.5625`,c44->1.4896636363636364`,c55->2.3409`,c66->3.27726`,c12->2.3050157470049792`,c23->5.297420162994092`,c13->3.1749113238217275`};

dxdσ = D[F,p1]//Simplify;
dydσ = D[F,p2]//Simplify;
dzdσ = D[F,p3]//Simplify;

depth1 = 0.25; depth2 = 0.45;depth3=0.3; (*depth in km*)
depth = depth1 + depth2+depth3;
dx1 =Δd1 dxdσ /dzdσ /.{p3->Sqrt[p3sqf1]}/.{c11->3.512016`,c22->3.875328`,c33->3.0276`,c44->0.1669390243902439`,c55->0.1521`,c66->0.27378`,c12->3.298434960000476`,c23->2.9819605693294067`,c13->2.8709922298203208`};
dx2 =Δd2 dxdσ /dzdσ /.{p3->Sqrt[p3sqf2]}/.{c11->4.51632`,c22->4.817408`,c33->3.7636`,c44->0.64896`,c55->0.6084`,c66->0.97344`,c12->2.995458418291907`,c23->2.3894757760926923`,c13->2.6577570205977543`};
dx3 =Δd3 dxdσ /dzdσ /.{p3->Sqrt[p3sqf3]}/.{c11->8.7725`,c22->9.68`,c33->7.5625`,c44->1.4896636363636364`,c55->2.3409`,c66->3.27726`,c12->2.3050157470049792`,c23->5.297420162994092`,c13->3.1749113238217275`};
dy1 = Δd1 dydσ /dzdσ /.{p3->Sqrt[p3sqf1]}/.{c11->3.512016`,c22->3.875328`,c33->3.0276`,c44->0.1669390243902439`,c55->0.1521`,c66->0.27378`,c12->3.298434960000476`,c23->2.9819605693294067`,c13->2.8709922298203208`};
dy2 = Δd2 dydσ /dzdσ /.{p3->Sqrt[p3sqf2]}/.{c11->4.51632`,c22->4.817408`,c33->3.7636`,c44->0.64896`,c55->0.6084`,c66->0.97344`,c12->2.995458418291907`,c23->2.3894757760926923`,c13->2.6577570205977543`};
dy3 = Δd3 dydσ /dzdσ /.{p3->Sqrt[p3sqf3]}/.{c11->8.7725`,c22->9.68`,c33->7.5625`,c44->1.4896636363636364`,c55->2.3409`,c66->3.27726`,c12->2.3050157470049792`,c23->5.297420162994092`,c13->3.1749113238217275`};
dt1 =Δd1 ((p1 dxdσ+p2 dydσ)/dzdσ +p3)/.{p3->Sqrt[p3sqf1]}/.{c11->3.512016`,c22->3.875328`,c33->3.0276`,c44->0.1669390243902439`,c55->0.1521`,c66->0.27378`,c12->3.298434960000476`,c23->2.9819605693294067`,c13->2.8709922298203208`};
dt2 =Δd2 ((p1 dxdσ+p2 dydσ)/dzdσ +p3)/.{p3->Sqrt[p3sqf2]}/.{c11->4.51632`,c22->4.817408`,c33->3.7636`,c44->0.64896`,c55->0.6084`,c66->0.97344`,c12->2.995458418291907`,c23->2.3894757760926923`,c13->2.6577570205977543`};
dt3 =Δd3 ((p1 dxdσ+p2 dydσ)/dzdσ +p3)/.{p3->Sqrt[p3sqf3]}/.{c11->8.7725`,c22->9.68`,c33->7.5625`,c44->1.4896636363636364`,c55->2.3409`,c66->3.27726`,c12->2.3050157470049792`,c23->5.297420162994092`,c13->3.1749113238217275`};
dx = (dx1 + dx2 + dx3)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3};
dy = (dy1 + dy2 + dy3)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3};
dt = (dt1 + dt2 + dt3)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3};
depth1 = 0.25; depth2 = 0.45;depth3=0.3; (*depth in km*)
depth = depth1 + depth2+depth3;
vert=({2dt1,2dt2,2dt3})/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->0,p2->0};
vertt= (2dt)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->0,p2->0};
vertn= vert/vertt;

vyav = Sqrt[vertn.{3.876,4.818,9.68}];
vxav = Sqrt[vertn.{3.52,4.52,8.78}];

maxx = 0.64/vxav ;
maxy = 0.64/vyav;
maxr = R maxx maxy/Sqrt[(maxy^2 Cos[theta]^2 + maxx^2 Sin[theta]^2)]; (*Ellipse in polar*)
xoffset =(2 dx)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};//Simplify
yoffset = (2dy)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};//Simplify
exacttime = (2dt)/.{Δd1->depth1,Δd2->depth2,Δd3->depth3}/.{p1->maxr Cos[theta],p2->maxr Sin[theta]};

twtime =Rasterize[ParametricPlot3D[{xoffset/depth,yoffset/depth,exacttime^2},{R,0,1},{theta,0,Pi/2},PlotRange->{All,All,All},BoxRatios->{1,1,0.3},
AxesLabel->{"x1(km)","x2(km)","4 t^2 (s^2)              "},ColorFunction->ColorData[{"DarkRainbow",{0.9,1.65}}],ColorFunctionScaling->False,LabelStyle->{Black,FontWeight->"Bold",FontSize->24},TicksStyle->Directive[Thick,20],Mesh->9,ImageSize->700]]

Export["junk_ma.eps",twtime,"EPS"]
junk_ma.eps
