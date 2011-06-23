function eicdemo(ldis,lang,tsf,prefl,pqvec,psouw,precw)

clf
eicframe();
eicgeom();
hold on
view(119,-12)

% reflector normal
nn = eicnorm(nn);

% azimuth reference;
vv = eicnorm(vv);
aa=cross(cross(nn,vv),nn);
aa = eicnorm(aa);

% source plane
ns = eicnorm(ns);

% receiver plane
nr = eicsnell(ns,nn);
nr = eicnorm(nr);

% in-plane "reflection" vector
qq = cross(cross(nn,ns),nn);
qq = eicnorm(qq);

% space-lag vector (in the reflection plane)
R=rotationmat3D(pi*lang/180,nn);
ll=ldis*(R*qq')';

% reflection angle
sth=norm(cross(nn,ns));
th=180*asin(sth)/pi;

% time shift
if tsf=='y'
  tt=-sth*(qq*ll')/v;
else
  tt=0;
end

zsft=[0 0 0.05];

if(prefl=='y')
  eicvect(2.0*nn,black,'n');      % reflector normal
  eicvect(4.0*vv,brown,'v');      % azimuth vector
  eicvect(2.0*aa+zsft,brown,'a');   % azimuth projection

  eicplane(nn,0,[0 0 0],5.0,1.5,'n'); % reflector
end

if(pqvec=='y')
  eicvect(2.0*qq+zsft,brown,'q');
end

if(psouw=='y')
  % source plane
  eicvect(2.0*ns,red,'n_s');
  if(ldis>0)
    eicvect(+ll,yellow,'\lambda');
  end
  eicplane(ns,+tt,+ll,1.5,0.5,'n');
end

if(precw=='y')
  % receiver plane
  eicvect(2.0*nr,blue,'n_r');
  if(ldis>0)
    eicvect(-ll,yellow,'\lambda');
  end
  eicplane(nr,-tt,-ll,1.5,2.5,'n');
end

