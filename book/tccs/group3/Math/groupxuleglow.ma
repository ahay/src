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

vp0 = Sqrt[c33];
eps2 = (c11 - c33)/(2*c33);
delta2 = ((c13 + c55)^2 - (c33 - c55)^2)/(2*c33*(c33 - c55));
eps1 = (c22 - c33)/(2*c33);
delta1 = ((c23 + c44)^2 - (c33 - c44)^2)/(2*c33*(c33 - c44));
delta3 = ((c12 + c66)^2 - (c11 - c66)^2)/(2*c11*(c11 - c66));

%Group velocity
t02=.;eta=.;eta2=.;
vnmosq = 1/(sin2/vnmo12 + cos2/vnmo22);
eta = eta2* cos2 - eta3 *cos2 * sin2+eta1 *sin2;
tsq= t02 +xhat2/vnmosq -2*eta*xhat2^2/(vnmosq*(t02*vnmosq + xhat2*(1+2* eta)));
tsqsub =tsq /.{sin2->pn22/(pn12+pn22),cos2->pn12/(pn12+pn22)};

t02f =pn32/vp0^2 /.{theta->theta*Pi/180,phi->phi*Pi/180}; 
xhat2f =pn12+pn22/.{theta->theta*Pi/180,phi->phi*Pi/180};
vnmo1sq = vp0^2 *(1+2 *delta1);
vnmo2sq =  vp0^2 *(1+2 *delta2);
eta1f = (eps1-delta1)/(1+2 *delta1);
eta2f = (eps2-delta2)/(1+2 *delta2);
eta3f = (eps1 - eps2 - delta3*(1+2* eps2))/((1+2*delta3)(1+2 *eps2));
tsqsubnew = tsqsub /.{t02->t02f,xhat2->xhat2f};
vgroupsq = 1/tsqsubnew;
vgroupsqxu = vgroupsq /.{vnmo12->vnmo1sq,vnmo22->vnmo2sq,eta1->eta1f,eta2->eta2f,eta3->eta3f};
groupxu =100*Abs[Sqrt[vgroupsqxu/vgp]-1]/.{n1->Sin[theta*Pi/180]Cos[phi*Pi/180],n2->Sin[theta*Pi/180]Sin[phi*Pi/180],n3->Cos[theta*Pi/180],c11->9,c33->5.9375,c55->1.6,c13->2.25,c22->9.84,c44->2,c12->3.6,c23->2.4,c66->2.182};

shift =0;
orthogroupplot=Graphics[{Text[Style[Phi,FontSize->16],{shift+110Cos[Pi/4],shift+110Sin[Pi/4]}],Text[Style[Theta,FontSize->16],{shift,shift+100}],Text[Style[Theta,FontSize->16],{shift+100,shift}],{Thick,Arrow[BSplineCurve[Table[{shift+100Cos[i],shift+ 100 Sin[i]},{i,Pi/9,7Pi/18,Pi/90}]]]},{Dashed,Table[Line[{{shift,shift},{shift + 90Cos[i],shift+90 Sin[i]}}],{i,0,2 Pi,Pi/6}]},PointSize[Large],Table[Point[{shift+ 90Cos[i],shift+ 90 Sin[i]}],{i,0,2 Pi,Pi/6}],{Thick,Circle[{shift,shift},90]}},PlotRange->{{-90+shift,90+shift},{-90+shift,90+shift}},Axes->True,AxesOrigin->{shift,shift},AxesStyle->Directive[Bold,Black,Thick],LabelStyle->{FontWeight->"Bold",FontSize->22},Ticks->{{{-60+shift,"60"},{-30+shift,"30"},{0+shift,"0"},{30+shift,"30"},{60+shift,"60"}},{{-60+shift,"60"},{-30+shift,"30"},{0+shift,"0"},{30+shift,"30"},{60+shift,"60"}}},ImagePadding->50];
groupx = Sign[Cos[phi*Pi/180]]Sqrt[pn12]*ArcSin[Sqrt[pn12+pn22]]*180/Pi/.{n1->Sin[theta*Pi/180]Cos[phi*Pi/180],n2->Sin[theta*Pi/180]Sin[phi*Pi/180],n3->Cos[theta*Pi/180],c11->9,c33->5.9375,c55->1.6,c13->2.25,c22->9.84,c44->2,c12->3.6,c23->2.4,c66->2.182};
groupy =Sign[Sin[phi*Pi/180]] Sqrt[pn22]*ArcSin[Sqrt[pn12+pn22]]*180/Pi/.{n1->Sin[theta*Pi/180]Cos[phi*Pi/180],n2->Sin[theta*Pi/180]Sin[phi*Pi/180],n3->Cos[theta*Pi/180],c11->9,c33->5.9375,c55->1.6,c13->2.25,c22->9.84,c44->2,c12->3.6,c23->2.4,c66->2.182};

ga =Show[Rasterize[ParametricPlot3D[{groupx,groupy,groupxu},{theta,0.1,89.9},{phi,0,360},PlotRange->{Full,Full,All},Axes->{False,False,False},Boxed->False,LabelStyle->{FontWeight->"Bold",FontSize->22},ColorFunction->ColorData[{"DarkRainbow",{0,2.5}}],ViewPoint->{0,0,Infinity},Mesh->None,PlotPoints->200,ColorFunctionScaling->False,ImagePadding->40],Background->None,ImageSize->600],Rasterize[orthogroupplot,Background->None,ImageSize->600]]
Export["junk_ma.eps",ga,"EPS"]
junk_ma.eps
