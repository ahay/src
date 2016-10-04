function bipolar3dtest_hti(cc,threshold,mode)
clf

% polar coordinates
R=3.14;
ot=0; dt=5; nt=(180-2*ot)/dt+1;  jt=1; kt=13;%latitude
oh=0;   dh=10; nh=360/dh+1; jh=1;%longitude
for     it=1:nt; t=pi*(ot+(it-1)*dt)/180;
    for ih=1:nh; h=pi*(oh+(ih-1)*dh)/180.;
        kx(it,ih) = R*sin(t) * cos(h);
        ky(it,ih) = R*sin(t) * sin(h);
        kz(it,ih) = R*cos(t);
    end
end
min(min(kz));
max(max(kz));
kk=sqrt(kx.^2 + ky.^2 + kz.^2);


% sort eigenvectors according to corresponding eignevalues, and P vector is kinda aligned
% with wave vector 'k'
k = [kx(1,1) ky(1,1) kz(1,1)];
[p q r] = cseed(cc,k);


abnit=[];abnih=[];nabn=1;
lat=[];lon=[];
for it=1:nt;
   for ih=1:nh;
        k = [kx(it,ih) ky(it,ih) kz(it,ih)];

%         p = [ax(it,ih-1) ay(it,ih-1) az(it,ih-1)];
%         q = [bx(it,ih-1) by(it,ih-1) bz(it,ih-1)];
%         r = [cx(it,ih-1) cy(it,ih-1) cz(it,ih-1)];

        [a b c dd] = christofel3dtest(cc,k);%a, b,c are sorted by their corresponding eigenvalues
        
%       (find(abs(b)==1) | find (abs(c)==1))  &
%        if(abs(dd(1)-dd(2))<threshold | abs(dd(1)-dd(3))<threshold | abs(dd(2)-dd(3))<threshold)
        if(   (abs(dd(1)-dd(2))<threshold | abs(dd(1)-dd(3))<threshold | abs(dd(2)-dd(3))<threshold)  )
            lat=[lat pi*(ot+(it-1)*dt)/180];
            lon=[lon pi*(oh+(ih-1)*dh)/180];
            abnit=[abnit it];
            abnih=[abnih ih];
%            b=[0 ;0 ;0];
%            c=[0 ;0 ;0];
            nabn=nabn+1;
        end
 
        [a b c] = vsign(a',b',c',p,q,r);
        %[a b c] = cneighbor(cc,k,p,q,r);

        a=a.*kk(it,ih);
        ax(it,ih)=a(1);
        ay(it,ih)=a(2);
        az(it,ih)=a(3);

        b=b.*kk(it,ih);
        bx(it,ih)=b(1);
        by(it,ih)=b(2);
        bz(it,ih)=b(3);

        c=c.*kk(it,ih);
        cx(it,ih)=c(1);
        cy(it,ih)=c(2);
        cz(it,ih)=c(3);

        p=a;
        q=b;
        r=c;
    end

    % re-seed
    %p=[ax(it,1) ay(it,1) az(it,1)];
    %q=[bx(it,1) by(it,1) bz(it,1)];
    %r=[cx(it,1) cy(it,1) cz(it,1)];

    %p=[ax(1,ih) ay(1,ih) az(1,ih)];
    %q=[bx(1,ih) by(1,ih) bz(1,ih)];
    %r=[cx(1,ih) cy(1,ih) cz(1,ih)];

end

% ax(kt,:)=0; ay(kt,:)=0; az(kt,:)=0;
% bx(kt,:)=0; by(kt,:)=0; bz(kt,:)=0;
% cx(kt,:)=0; cy(kt,:)=0; cz(kt,:)=0;


% --------------------------------------------------
% PLOTS
% --------------------------------------------------
ccll=ones(64,3);
colormap(ccll)

[sa,sb,sc]=sphere();
surf(sa*R,sb*R,sc*R);
hold on
quiver3(1*R,0,0,1,0,0,.25*R,'y','LineWidth',3);text(1.25*R,0,0,'x','BackgroundColor',[.7 .9 .7])
quiver3(0,1*R,0,0,1,0,.25*R,'y','LineWidth',3);text(0,1.25*R,0,'y','BackgroundColor',[.7 .9 .7])
quiver3(0,0,1*R,0,0,1,.25*R,'y','LineWidth',3);text(0,0,1.25*R,'z','BackgroundColor',[.7 .9 .7])

%nabn=length(abnit);
%for ii=1:nabn
%      plot3(kx(abnit(ii),abnih(ii)),ky(abnit(ii),abnih(ii)),kz(abnit(ii),abnih(ii)),'oy');        
%end


if(mode=='k')
    h=quiver3(kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh),kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh));
    set(h, 'Color', 'black','LineWidth',1);
end

if( mode=='p' )
    h=quiver3(kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh),ax(1:jt:nt,1:jh:nh),ay(1:jt:nt,1:jh:nh),az(1:jt:nt,1:jh:nh));
    set(h, 'Color', 'blue'  ,'LineWidth',1);
end

if( mode=='s1' )
    h=quiver3(kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh),...
        bx(1:jt:nt,1:jh:nh),by(1:jt:nt,1:jh:nh),bz(1:jt:nt,1:jh:nh));
    set(h, 'Color', 'blue','LineWidth',1);
end

if( mode=='s2' )
     h=quiver3(kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh),cx(1:jt:nt,1:jh:nh),cy(1:jt:nt,1:jh:nh),cz(1:jt:nt,1:jh:nh));
    set(h, 'Color', 'blue','LineWidth',1);
end

if( mode=='sv' )
    h=quiver3(kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh),...
           ay(1:jt:nt,1:jh:nh).*ay(1:jt:nt,1:jh:nh)+az(1:jt:nt,1:jh:nh).*az(1:jt:nt,1:jh:nh),-ax(1:jt:nt,1:jh:nh).*ay(1:jt:nt,1:jh:nh),-ax(1:jt:nt,1:jh:nh).*az(1:jt:nt,1:jh:nh))
    set(h, 'Color', 'blue','LineWidth',1);
end

if( mode=='sh' )
    h=quiver3(kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh),zeros(nt/jt,nh/jh),-az(1:jt:nt,1:jh:nh),ay(1:jt:nt,1:jh:nh))
    %h=quiver3(kx(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh),kz(1:jt:nt,1:jh:nh),zeros(nt/jt,nh/jh),-kz(1:jt:nt,1:jh:nh),ky(1:jt:nt,1:jh:nh))
    set(h, 'Color', 'blue','LineWidth',1);
end


axis([-R R -R R -R R])
daspect([1 1 1])
view(120,40)

singlat=lat;
singlon=lon;

xlabel('k_x(radians)'); set(get(gca,'XLabel'),'Rotation',50,'fontsize',20,'position',[0,4,-4])
ylabel('k_y(radians)'); set(get(gca,'YLabel'),'Rotation',-18,'fontsize',20,'position',[4,0,-4])
zlabel('k_z(radians)'); set(get(gca,'ZLabel'),'fontsize',20)

set(gca,'xtick',[-3  -1 1  3])
set(gca,'ytick',[-3  -1 1  3])
set(gca,'ztick',[-3  -1 1  3])
