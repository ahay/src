%Orthorhombic approximations
%1. Exact Expression
ClearAll[c11,c12,c13,c12,c23,c22,c33,c44,c55,c66,n1,n2,n3,two,trace,G,H,angle,vhelbig,vpphase,vecdiff,vpgroup,vpg,pn12,pn22,pn32]
mat = {{c11,c12,c13,0,0,0},{c12,c22,c23,0,0,0},{c13,c23,c33,0,0,0},{0,0,0,c44,0,0},{0,0,0,0,c55,0},{0,0,0,0,0,c66}};
matn = {{n1,0,0,0,n3,n2},{0,n2,0,n3,0,n1},{0,0,n3,n2,n1,0}} ;
chris = matn.mat.Transpose[matn];
%Exact qP phase velocity from schoenberg and helbig (1997)
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
Phase velocity
α=1;
β=1;
ms = ((c11 n1^2)^α*( msh2f*(c33 n3^2)^α+ msv3f*(c22 n2^2)^α)/((c22 n2^2 )^α+(c33 n3^2)^α) + (c22 n2^2)^α*( msh1f*(c33 n3^2)^α+msh3f*( c11 n1^2)^α)/((c11 n1^2)^α +(c33 n3^2)^α)+(c33 n3^2)^α*( msv1f*(c22 n2^2)^α+ msv2f*(c11 n1^2)^α)/((c11 n1^2)^α+ (c22 n2^2)^α))/((c11 n1^2)^α + (c22 n2^2)^α +(c33 n3^2)^α);

mq1 = ( mqv1f*(c33 n3^2)^β+mqh1f*(c22 n2^2)^β)/((c22 n2^2)^β +(c33 n3^2)^β);
mq2 = (mqv2f*(c33 n3^2)^β+mqh2f*(c11 n1^2)^β)/((c11 n1^2)^β+(c33 n3^2)^β);
mq3 = (mqv3f*(c11 n1^2)^β+mqh3f*(c22 n2^2)^β)/((c11 n1^2)^β + (c22 n2^2)^β); 

el=c11*n1^2+c22*n2^2+c33*n3^2;
msh =(1-ms)*el+ms*Sqrt[el^2+2*((mq1-1)*c22*c33*n2^2*n3^2+(mq2-1)*c11*c33*n1^2*n3^2+(mq3-1)*c11*c22*n1^2*n2^2)/ms];
mshsub = msh /.{n1->Sin[theta]*Cos[phi],n2->Sin[theta]*Sin[phi],n3->Cos[theta]};

mqv1=(c23^2+2 c23 c44+c33 c44)/(c22 c33-c22 c44);
mqh1=(c23^2+c22 c44+2 c23 c44)/(c22 c33-c33 c44);
mqv2=(c13^2+2 c13 c55+c33 c55)/(c11 c33-c11 c55);
mqh2=(c13^2+c11 c55+2 c13 c55)/(c11 c33-c33 c55);
mqv3=(c12^2+c11 c66+2 c12 c66)/(c11 c22-c22 c66);
mqh3=(c12^2+2 c12 c66+c22 c66)/(c11 c22-c11 c66);

msv2q=((c11-c33) (-1+qh2c) (-1+qv2c)^2)/(2 (-c33 (-1+2 qv2c+qh2c (1+qh2c+(-4+qv2c) qv2c))+c11 (-1+qh2c^2-2 qh2c qv2c+qv2c (3+(-2+qv2c) qv2c))));
msh2q=((c11-c33) (-1+qh2c)^2 (-1+qv2c))/(2 (c33-c33 (qh2c (3+(-2+qh2c) qh2c)-2 qh2c qv2c+qv2c^2)+c11 (-1+2 qh2c+qv2c+(-4+qh2c) qh2c qv2c+qv2c^2)));
msv1q=((c22-c33) (-1+qh1c) (-1+qv1c)^2)/(2 (-c33 (-1+2 qv1c+qh1c (1+qh1c+(-4+qv1c) qv1c))+c22 (-1+qh1c^2-2 qh1c qv1c+qv1c (3+(-2+qv1c) qv1c))));
msh1q=((c22-c33) (-1+qh1c)^2 (-1+qv1c))/(2 (c33-c33 (qh1c (3+(-2+qh1c) qh1c)-2 qh1c qv1c+qv1c^2)+c22 (-1+2 qh1c+qv1c+(-4+qh1c) qh1c qv1c+qv1c^2)));
msv3q=((c11-c22) (-1+qh3c) (-1+qv3c)^2)/(2 (c22+c11 (-1+2 qv3c+qh3c (1+qh3c+(-4+qv3c) qv3c))-c22 (qh3c^2-2 qh3c qv3c+qv3c (3+(-2+qv3c) qv3c))));
msh3q=((c11-c22) (-1+qh3c)^2 (-1+qv3c))/(2 (c22+c11 (-1+qh3c (3+(-2+qh3c) qh3c)-2 qh3c qv3c+qv3c^2)-c22 (2 qh3c+qv3c+(-4+qh3c) qh3c qv3c+qv3c^2)));
%True q
mshsubtrue=msh/.{theta->theta*Pi/180,phi->phi*Pi/180};
phasex =theta Cos[phi*Pi/180];
phasey =theta Sin[phi*Pi/180];
%Approximated q
shift =0;
orthoplot=Graphics[{Text[Style[phi,FontSize->16],{shift+110Cos[Pi/4],shift+110Sin[Pi/4]}],Text[Style[theta,FontSize->16],{shift,shift+100}],Text[Style[theta,FontSize->16],{shift+100,shift}],{Thick,Arrow[BSplineCurve[Table[{shift+100Cos[i],shift+ 100 Sin[i]},{i,Pi/9,7Pi/18,Pi/90}]]]},{Dashed,Table[Line[{{shift,shift},{shift + 90Cos[i],shift+90 Sin[i]}}],{i,0,2 Pi,Pi/6}]},PointSize[Large],Table[Point[{shift+ 90Cos[i],shift+ 90 Sin[i]}],{i,0,2 Pi,Pi/6}],{Thick,Circle[{shift,shift},90]}},PlotRange->{{-90+shift,90+shift},{-90+shift,90+shift}},Axes->True,AxesOrigin->{shift,shift},AxesStyle->Directive[Bold,Black,Thick],LabelStyle->{FontWeight->"Bold",FontSize->22},Ticks->{{{-60+shift,"60"},{-30+shift,"30"},{0+shift,"0"},{30+shift,"30"},{60+shift,"60"}},{{-60+shift,"60"},{-30+shift,"30"},{0+shift,"0"},{30+shift,"30"},{60+shift,"60"}}},ImagePadding->50];
mshsubtrue=msh/.{theta->theta*Pi/180,phi->phi*Pi/180};
phasex =theta Cos[phi*Pi/180];
phasey =theta Sin[phi*Pi/180];
phasemshappq =100*Abs[Sqrt[mshsubtrue/vhelbig]-1]/.{msv1f->msv1q,msh1f->msh1q,msv2f->msv2q,msh2f->msh2q,msv3f->msv3q,msh3f->msh3q,mqv1f->mqv1,mqh1f->0.8373406986*mqv1+0.1580976683,mqv2f->mqv2,mqh2f->0.8373406986*mqv2+0.1580976683,mqv3f->mqv3,mqh3f->mqv3}/.{qv1c->mqv1,qh1c->0.8373406986*mqv1+0.1580976683,qv2c->mqv2,qh2c->0.8373406986*mqv2+0.1580976683,qv3c->mqv3,qh3c->((c11(c22-c66))mqv3+(c13/c23-1)c66(c22+c12))/(c22(c11-c66))}/.{n1->Sin[theta*Pi/180]Cos[phi*Pi/180],n2->Sin[theta*Pi/180]Sin[phi*Pi/180],n3->Cos[theta*Pi/180],c11->9,c33->5.9375,c55->1.6,c13->2.25,c22->9.84,c44->2,c12->3.6,c23->2.4,c66->2.182};
ga =Show[Rasterize[ParametricPlot3D[{phasex,phasey,phasemshappq},{theta,0,90},{phi,0,360},PlotRange->{Full,Full,All},Axes->{False,False,False},Boxed->False,LabelStyle->{FontWeight->"Bold",FontSize->22},ColorFunction->ColorData[{"DarkRainbow",{0,0.5}}],ViewPoint->{0,0,Infinity},Mesh->None,PlotPoints->200,ColorFunctionScaling->False,ImagePadding->40],Background->None,ImageSize->600],Rasterize[orthoplot,Background->None,ImageSize->600]]
Export["junk_ma.eps",ga,"EPS"]
junk_ma.eps
c11->6.3,c33->5.411,c55->0.8,c13->2.25,c22->6.871,c44->1,c12->2.7,c23->2.393,c66->1.5
