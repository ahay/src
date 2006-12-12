eq1=(x-x2)^2+(z-z2)^2==R2^2;
eq2=(x2-x)(x3-x1)+(z2-z)(z3-z1)==R2 (R3-R1);
c[x_,z_,r_, a_]:= Graphics[{Circle[{x,z},r], Disk[{x,z},0.2], 
			   Text[a,{x+0.5,z+0.5}]}];
p[x_,z_,a_]:= Graphics[{Disk[{x,z},0.2], Text[a,{x+0.5,z+0.5}]}];  
cr[x_,z_,r_]:=Graphics[{AbsoluteThickness[3], Circle[{x,z},r]}];
$TextStyle = {FontFamily->"Helvetica", FontWeight->Bold, FontSize->14};
ln[x_, x1_,x2_,z1_,z2_]:=x (z1-z2)/(x1-x2) + (z2 x1 - z1 x2)/(x1-x2);
abline = Plot[ln[x,1, 4  , 1, -2],{x,1,4  }];
bcline = Plot[ln[x,4, 4.5,-2,  2],{x,4,4.5}];
<< Graphics`Arrow`
obarrow = Show[Graphics[Arrow[{8,4},{4,-2}]]];
Show[p[8,4,"O"], c[1,1, 1, "A"], c[4, -2, 2, "B"],  c[4.5,2,3, "C"],
abline, bcline, obarrow,  AspectRatio->Automatic];
Display["junk_ma.eps",%,"EPS"];
