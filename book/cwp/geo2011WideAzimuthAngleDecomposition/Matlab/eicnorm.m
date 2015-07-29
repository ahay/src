function bb=eicnorm(aa)

ax=aa(1);
ay=aa(2);
az=aa(3);

am=sqrt(aa*aa');

%bx=ax/am;
%by=ay/am;
%bz=az/am;
%bb=[bx by bz];

bb=aa./am;
