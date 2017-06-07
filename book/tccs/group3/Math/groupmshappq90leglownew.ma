%Orthorhombic approximations
%1. Exact Expression
ClearAll[c11,c12,c13,c12,c23,c22,c33,c44,c55,c66,n1,n2,n3,two,trace,G,H,angle,vhelbig,vpphase,vecdiff,vpgroup,vpg,pn12,pn22,pn32]
mat = {{c11,c12,c13,0,0,0},{c12,c22,c23,0,0,0},{c13,c23,c33,0,0,0},{0,0,0,c44,0,0},{0,0,0,0,c55,0},{0,0,0,0,0,c66}};
matn = {{n1,0,0,0,n3,n2},{0,n2,0,n3,0,n1},{0,0,n3,n2,n1,0}} ;
chris = matn.mat.Transpose[matn];
Exact qP phase velocity from schoenberg and helbig (1997)
two = (chris[[1,1]]chris[[2,2]] + chris[[1,1]]chris[[3,3]] + chris[[2,2]]chris[[3,3]] - chris[[1,2]]^2 - chris[[1,3]]^2 - chris[[2,3]]^2) //FullSimplify;
trace = (chris[[1,1]] +chris[[2,2]] + chris[[3,3]]) ;
G = (trace^2)/9 - two/3 ;
H = (trace^3)/27 -(two *trace/6) + Det[chris]/2;
angle = ArcCos[H/Sqrt[G^3]];
vhelbig = 2 Sqrt[G]Cos[angle/3] + trace/3;
n1=.;n2=.; n3=.;
vecn = {n1,n2,n3} ;
vpphase = Sqrt[vhelbig];
vecdiff = {D[vpphase,n1],D[vpphase,n2],D[vpphase,n3]};
vpgroup = vpphase vecn  +(IdentityMatrix[3]-Transpose[{vecn}].{vecn}).vecdiff;
vgp = vpgroup[[1]]^2 + vpgroup[[2]]^2 +vpgroup[[3]]^2 ;
pn12 = (vpgroup[[1]]^2)/vgp;
pn22 = (vpgroup[[2]]^2)/vgp;
pn32 = (vpgroup[[3]]^2)/vgp;

%2. Modified Muir-Dellinger Approximation with shifted hyperbola
α=1;
β=1;
mqv1=(c23^2+2 c23 c44+c33 c44)/(c22 c33-c22 c44);
mqh1=(c23^2+c22 c44+2 c23 c44)/(c22 c33-c33 c44);
mqv2=(c13^2+2 c13 c55+c33 c55)/(c11 c33-c11 c55);
mqh2=(c13^2+c11 c55+2 c13 c55)/(c11 c33-c33 c55);
mqv3=(c12^2+c11 c66+2 c12 c66)/(c11 c22-c22 c66);
mqh3=(c12^2+2 c12 c66+c22 c66)/(c11 c22-c11 c66);
mQv1=1/mqv1;
mQh1=1/mqh1;
mQv2=1/mqv2;
mQh2=1/mqh2;
mQv3=1/mqv3;
mQh3=1/mqh3;
mSv2Q=((c11-c33) (-1+Qh2c) (-1+Qv2c)^2)/(2 c11 (-1+Qh2c^2+2 Qv2c+Qh2c Qv2c (-1+(-2+Qv2c) Qv2c))-2 c33 (-1+Qh2c^2-2 Qh2c Qv2c+Qv2c (2+(-1+Qv2c)^2 Qv2c)));
mSh2Q=((c11-c33) (-1+Qh2c)^2 (-1+Qv2c))/(2 (c33+c11 (-1+Qh2c (2+(-1+Qh2c)^2 Qh2c)-2 Qh2c Qv2c+Qv2c^2)-c33 (2 Qh2c+Qh2c (-1+(-2+Qh2c) Qh2c) Qv2c+Qv2c^2)));
mSv1Q=((c22-c33) (-1+Qh1c) (-1+Qv1c)^2)/(2 c22 (-1+Qh1c^2+2 Qv1c+Qh1c Qv1c (-1+(-2+Qv1c) Qv1c))-2 c33 (-1+Qh1c^2-2 Qh1c Qv1c+Qv1c (2+(-1+Qv1c)^2 Qv1c)));
mSh1Q=((c22-c33) (-1+Qh1c)^2 (-1+Qv1c))/(2 (c33+c22 (-1+Qh1c (2+(-1+Qh1c)^2 Qh1c)-2 Qh1c Qv1c+Qv1c^2)-c33 (2 Qh1c+Qh1c (-1+(-2+Qh1c) Qh1c) Qv1c+Qv1c^2)));
mSv3Q=((c11-c22) (-1+Qh3c) (-1+Qv3c)^2)/(2 (c22-c22 (Qh3c^2+2 Qv3c+Qh3c Qv3c (-1+(-2+Qv3c) Qv3c))+c11 (-1+Qh3c^2-2 Qh3c Qv3c+Qv3c (2+(-1+Qv3c)^2 Qv3c))));
mSh3Q=((c11-c22) (-1+Qh3c)^2 (-1+Qv3c))/(2 (c22-c22 (Qh3c (2+(-1+Qh3c)^2 Qh3c)-2 Qh3c Qv3c+Qv3c^2)+c11 (-1+2 Qh3c+Qh3c (-1+(-2+Qh3c) Qh3c) Qv3c+Qv3c^2)));
%Group velocity
El=1/c11*pn12+1/c22*pn22+1/c33*pn32;
MS =( (pn12/c11)^α*(mSh2f*(pn32/c33)^α+mSv3f*(pn22/c22)^α)/((pn22/c22)^α +(pn32/c33)^α) + (pn22/c22)^α*(mSh1f*(pn32/c33)^α+mSh3f*(pn12/c11)^α)/((pn12/c11)^α + (pn32/c33)^α)+(pn32/c33)^α*(mSv1f*(pn22/c22)^α+mSv2f*(pn12/c11)^α)/((pn22/c22)^α+ (pn12/c11)^α))/((pn12/c11)^α+(pn22/c22)^α+(pn32/c33)^α);

MQ1 = (mQv1f*(pn32/c33)^β+mQh1f*(pn22/c22)^β)/((pn22/c22)^β +(pn32/c33)^β);
MQ2 = (mQv2f*(pn32/c33)^β+mQh2f*(pn12/c11)^β)/((pn12/c11)^β + (pn32/c33)^β);
MQ3 = (mQv3f*(pn12/c11)^β+mQh3f*(pn22/c22)^β)/((pn22/c22)^β + (pn12/c11)^β); 

MSH =(1-MS)*El+MS*Sqrt[El^2+2*((MQ1-1)*(1/(c22*c33))*pn22*pn32+(MQ2-1)*(1/(c11*c33))*pn12*pn32+(MQ3-1)*(1/(c11*c22))*pn12*pn22)/MS];
MSHsub = MSH/.{n1->Sin[theta]*Cos[phi],n2->Sin[theta]*Sin[phi],n3->Cos[theta]};
mS =( Sh*(a n^2)^α+(b(1-n^2))^α*Sv)/((a n^2)^α+(b(1-n^2))^α);
mQ =(Qh*(a n^2)^β+(b(1-n^2))^β*Qv)/((a n^2)^β+(b(1-n^2))^β);
Fmsh = (1-mS)*(a*n^2-(n^2-1)*b)+(mS)*Sqrt[2*(1-n^2)*(mQ-1)*a*b*n^2/(mS)+(a*n^2-(n^2-1)*b)^2];
groupx = Sign[Cos[phi*Pi/180]]Sqrt[pn12]*ArcSin[Sqrt[pn12+pn22]]*180/Pi/.{n1->Sin[theta*Pi/180]Cos[phi*Pi/180],n2->Sin[theta*Pi/180]Sin[phi*Pi/180],n3->Cos[theta*Pi/180],c11->9,c33->5.9375,c55->1.6,c13->2.25,c22->9.84,c44->2,c12->3.6,c23->2.4,c66->2.182};
groupy =Sign[Sin[phi*Pi/180]] Sqrt[pn22]*ArcSin[Sqrt[pn12+pn22]]*180/Pi/.{n1->Sin[theta*Pi/180]Cos[phi*Pi/180],n2->Sin[theta*Pi/180]Sin[phi*Pi/180],n3->Cos[theta*Pi/180],c11->9,c33->5.9375,c55->1.6,c13->2.25,c22->9.84,c44->2,c12->3.6,c23->2.4,c66->2.182};

%Approximated q
groupmshappq =100*Abs[Sqrt[1/(vgp*MSH)]-1]/.{mSv1f->mSv1Q,mSh1f->mSh1Q,mSv2f->mSv2Q,mSh2f->mSh2Q,mSv3f->mSv3Q,mSh3f->mSh3Q,mQv1f->mQv1,mQh1f->mQv1/(0.8373406986+0.1580976683*mQv1),mQv2f->mQv2,mQh2f->mQv2/(0.8373406986+0.1580976683*mQv2),mQv3f->mQv3,mQh3f->mQv3}/.{Qv1c->mQv1,Qh1c->mQv1/(0.8373406986+0.1580976683*mQv1),Qv2c->mQv2,Qh2c->mQv2/(0.8373406986+0.1580976683*mQv2),Qv3c->mQv3,Qh3c->mQv3}/.{n1->Sin[theta*Pi/180]Cos[phi*Pi/180],n2->Sin[theta*Pi/180]Sin[phi*Pi/180],n3->Cos[theta*Pi/180],c11->9,c33->5.9375,c55->1.6,c13->2.25,c22->9.84,c44->2,c12->3.6,c23->2.4,c66->2.182};
ga =Rasterize[ParametricPlot3D[{groupx,groupy,groupmshappq},{theta,0.1,89.9},{phi,0,90},PlotRange->{Full,Full,All},AxesLabel->{"Theta  ","  Theta","Relative\nError (%)   "},LabelStyle->{FontWeight->"Bold",FontSize->22},ColorFunction->ColorData[{"DarkRainbow",{0,0.5}}],TicksStyle->Directive[Thick,20],ViewAngle->Automatic,PlotLegends->Automatic,Mesh->10,ColorFunctionScaling->False,BoxRatios->{1,1,1/4},ImageSize->800],RasterSize->4000,ImageSize->600]

Export["junk_ma.eps",ga,"EPS"]
