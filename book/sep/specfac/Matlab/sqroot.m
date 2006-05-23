function sqroot(s, s0, niter,b)
clf; hold on;

c=['k ';'r ';'b ';'g ';'c ';'m '];
p=['ko-';'rx-';'b*-';'gd-';'cd-';'m>-'];

ii=1:1:niter;
tit=['muir  ';'secant';'newton';'ideal '];

N=4; 
xx=[];
xx(1)=s0;
for k=1:N
    for i=2:niter
	if	k==1 a=1;                           % muir
	elseif  k==2 xx(2)=xx(1); a=xx(max(1,i-2)); % secant
	elseif  k==3 a=xx(i-1);                     % newton
	else         a=1.25*sqrt(s);
	end
	xx(i)=(xx(i-1)*a+s)/(xx(i-1)+a);
    end
    plot(ii,xx,p(k,:)); 
end

grid on; axis([min(ii) max(ii) 5 15]);
legend('muir','secant','newton','ideal');











