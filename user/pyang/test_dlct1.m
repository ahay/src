function test_dlct

% Reference: Alkishriwo, Osama A., and Luis F. Chaparro. "A discrete 
%  linear chirp transform (DLCT) for data compression." Information 
%  Science, Signal Processing and their Applications (ISSPA), 2012 
%  11th International Conference on. IEEE, 2012.

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


mask=ones(N,L);
for ix=1:L
    for iy=1:N
        if (iy> -500*(ix-100)/60 || iy<-10*(ix-55))
            mask(iy,ix)=0;
        end
    end
end

XX1=dlct(x, C, L,N);
[Y,ind]=sort(abs(XX1),'ascend');
thr=Y(floor(0.1*N*L));
XX=thresholding(XX1,thr);
XX=ifftshift(mask.*fftshift(XX1,1),1);


u=fftshift(abs(XX1),1);
figure(2),clf
subplot(121)
mesh(u),shading interp
xlabel('l:[-L/2,L/2-1]')
ylabel('k:[-N/2:N/2-1]')
title('mesh plot of DLCT')
u1=fftshift(abs(XX),1);
subplot(122)
mesh(u1),shading interp




%------------------- reconstruct ---------------------
x_rec=idlct(XX,C,L,N);
% figure(3),clf
% %plot(real(x_rec))
% plot(real(x-x_rec))
% title('reconstruction error')

figure(4),clf
subplot(511)
plot(real(x1))
xlabel('Time:s'),ylabel('Amplitude')
title('signal component 1')
set(gca,'CLim',[-1,1])

subplot(512)
plot(real(x_rec))
xlabel('Time:s'),ylabel('Amplitude')
title('estimated component 1')
set(gca,'yLim',[-1,1])

subplot(513)
plot(real(x2))
xlabel('Time:s'),ylabel('Amplitude')
title('signal component 2')
set(gca,'CLim',[-1,1])

subplot(514)
plot(real(x-x_rec))
xlabel('Time:s'),ylabel('Amplitude')
title('estimated component 2')
set(gca,'yLim',[-1,1])

subplot(515)
plot(real(x))
xlabel('Time:s'),ylabel('Amplitude')
title('original chirp signal')


function y=thresholding(x,thr)
y=x.*(abs(x)>thr);


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
