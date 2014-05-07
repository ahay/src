function test_dlct
clc,clear, close all

Fs = 100;
dt = 1 / Fs ;
N= 500;n=1:N;
t=n'*dt;
x1=exp(j*2*pi* t.^2);%sin(2*pi* t.^2)
x2=exp(j*2*pi*5*t.^2);%sin(2*pi*5*t.^2)
x =x1 +x2;




% -------------decomposition-----------------------
L=100;
C=0.01;
XX=dlct(x, C, L,N);
u=fftshift(abs(XX),1);
u=abs(XX);
figure(1),clf
subplot(121)
mesh(u),shading interp
%pcolor(u), shading interp
xlabel('l:[-L/2,L/2-1]')
ylabel('k:[-N/2:N/2-1]')
title('(a) mesh plot of DLCT')

mask=ones(N,L);
for ix=1:L
    for iy=1:N
        if (iy> -500*(ix-100)/60 || iy<-10*(ix-55))
            mask(iy,ix)=0;
        end
    end
end
XX=ifftshift(mask.*fftshift(XX,1),1);
u1=fftshift(abs(XX),1);
subplot(122)
mesh(u1),shading interp
title('(b) masking the DLCT')
set(gcf,'PaperPosition',[0 0 12 6])
% print -depsc mesh_dlct.eps



figure(2),clf
subplot(311)
plot(t,real(x))
xlabel('Time:s'),ylabel('Amplitude')
title('(a) original chirp signal')

subplot(312)
plot(t,real(x1))
xlabel('Time:s'),ylabel('Amplitude')
title('(b) signal component 1')
set(gca,'CLim',[-1,1])

subplot(313)
plot(t,real(x2))
xlabel('Time:s'),ylabel('Amplitude')
title('(c) signal component 2')
set(gca,'CLim',[-1,1])


set(gcf,'PaperPosition',[0 0 12 6])
% print -depsc chirps.eps



%------------------- reconstruct ---------------------
x_rec=idlct(XX,C,L,N);

figure(4),clf
subplot(411)
plot(t,real(x1))
xlabel('Time:s'),ylabel('Amplitude')
title('(a) signal component 1')
set(gca,'CLim',[-1,1])

subplot(412)
plot(t,real(x_rec))
xlabel('Time:s'),ylabel('Amplitude')
title('(b) estimated component 1')
set(gca,'yLim',[-1,1])

subplot(413)
plot(t,real(x2))
xlabel('Time:s'),ylabel('Amplitude')
title('(c) signal component 2')
set(gca,'CLim',[-1,1])

subplot(414)
plot(t,real(x-x_rec))
xlabel('Time:s'),ylabel('Amplitude')
title('(d) estimated component 2')
set(gca,'yLim',[-1,1])

set(gcf,'PaperPosition',[0 0 12 6])
% print -depsc estimated.esp



function y=dlct(x,C, L,N)
%  forward discrete linear chirp transform
n=1:N;
h=zeros(N,1);
y=zeros(N,L);
for l=-L/2:L/2-1
    h=x.*exp(-j*2*pi*C*l*n'.^2/N);
    y(:,l+L/2+1)=fft(h);
end

function x=idlct(y,C,L,N)
%  forward discrete linear chirp transform
n=1:N;
g=zeros(N,L);
for l=-L/2:L/2-1
    g(:,l+L/2+1)=ifft(y(:,l+L/2+1));
end
x=sum(g.*exp(j*2*pi*C*n'.^2*[-L/2:L/2-1]/N),2)/L;
x=real(x);





