function eicplane(aa,tt,ll,rmax,colornum,small)
eicdims;

% disk normal
aa=eicnorm(aa);
ax=aa(1);
ay=aa(2);
az=aa(3);

% space lag
lx=ll(1);
ly=ll(2);
lz=ll(3);

if(small=='y')
  kx=(nx-1)/4; 
  ky=(ny-1)/4; 
  kz=(nz-1)/4;
else
  kx=0;
  ky=0;
  kz=0;
end

if az==max([ax ay az])
        %'loop over x and y'
	for ix=1+kx:nx-kx;
	  for iy=1+ky:ny-ky;
    		x(iy,ix)=ox+(ix-1)*dx;
    		y(iy,ix)=oy+(iy-1)*dy;
  	  end
	end
	z = lz+(v*tt -ax*(x-lx) -ay*(y-ly) )/az;
end

if ax==max([ax ay az])
        %'loop over y and z'
        for iy=1+ky:ny-ky;
          for iz=1+kz:nz-kz;
                y(iz,iy)=oy+(iy-1)*dy;
                z(iz,iy)=oz+(iz-1)*dz;
          end
        end
        x = lx+(v*tt -ay*(y-ly) -az*(z-lz) )/ax;
end

if ay==max([ax ay az])
        %'loop over z and x'
        for iz=1+kz:nz-kz;
          for ix=1+kx:nx-kx;
                z(ix,iz)=oz+(iz-1)*dz;
                x(ix,iz)=ox+(ix-1)*dx;
          end
        end
        y = ly+(v*tt -az*(z-lz) -ax*(x-lx) )/ay;
end

c = 0*z+colornum;
r = sqrt( (y-ly).^2+...
          (z-lz).^2+...
          (x-lx).^2 );
%c(r>rmax)=NaN;

x(r>rmax)=NaN;
y(r>rmax)=NaN;
z(r>rmax)=NaN;


%s=surface(x,y,z,c);
s=surf(x,y,z,c);
set(s,'edgecolor','none');
shading flat
