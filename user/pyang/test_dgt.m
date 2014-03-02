function test_dgt
%dgt.m: This programme is used to compute dual frame and reconstruct 
%           the original signal.
%Copyright(c) Xi'an Jiaotong University, Pengliang YANG,2010.5.1
%
% Reference: Qian, Shie, and Dapang Chen. "Discrete gabor transform."
%   Signal Processing, IEEE Transactions on 41.7 (1993): 2429-2438.
clc, clear,close all

Ls = 1024;          % length of the original signal s[n]
Ts = 1/1000;        % Sampling time: Ts = 1/Fs
t =[0:Ls-1]*Ts;
s = [sin(400*pi*t(1:256))   sin(800*pi*t(257:512))...
    sin(200*pi*t(513:768))    sin(600*pi*t(769:1024))];  % s[n]
figure(1),clf
subplot(221)
plot(s),axis tight
title('Original Signal')

L = 128;            % length of the frame basic function
[u,dM] = mydivider(Ls+L);% dM:step length in time field
M = (Ls+L)/dM-1;
Lo = Ls+L-dM;         
[ ,dN] = mydivider(M);
N = Lo/dN;         
r = Lo/(dM*dN);     % oversampling rate
fprintf('oversampling rate:r = %f\n',r)
sigma2 = dM*N/(2*pi);   
k = [0:L-1];
h = (pi*sigma2)^(-0.25)*exp(-(k-0.5*(L-1)).^2/(2*sigma2));

s = [zeros(1,Lo-Ls) s ]; % auxiliary periodic sequence for s(t)->s[n]
h = [h zeros(1,Lo-L)];  % auxiliary periodic sequence for h(t)->h[n]
h = h/norm(h);

H = zeros(dM*dN,Lo);
for q = 0:dN-1
    for p = 0:dM-1
        for k = 0:Lo-1
            indice = mymod(k+q*N+1,Lo);
            H(p+q*dM+1,k+1) = h(indice)*exp(-j*2*pi*p*k/dM);
        end
    end
end
u = [dM/N,zeros(1,dM*dN-1)]';
g = conj(pinv(H)*u);    % g(t) is the dual function of h(t).g(t)->g[n] 
gg = g/norm(g);         % normalization
subplot(222)
plot(h(1:L),'b')
hold on 
plot(real(gg(1:L)),'r'),axis tight
xlabel('k'),ylabel('Amplitude')
title('Basic Frame Fuction and Dual Frame Function')
legend('Basic Fuction','Dual Function')

    % ------------------Gabor Transform--------------------------
for m = 1:M
    for k = 1:Lo
        temp(k) = s(k)*conj(g(mymod(k-(m-1)*dM,Lo)));
    end
    temp = reshape(temp,N,dN);  % ----change dimension:Lo->N*dN
    C(m,:) = sum(fft(temp),2);  % ----line:N-FFT��sum in row.
end    
    % -----Inverse Transform(Reconstruct the original signal----- 
ss = zeros(1,Lo); 
for k = 1:Lo
    for m = 1:M
        InvC(m,:) = N*ifft(C(m,:));
        ss(k) = ss(k)+h(mymod(k-(m-1)*dM,Lo))*InvC(m,mymod(k,N));
    end
end
subplot(223)
plot(real(ss(Lo-Ls+1:Lo))),axis tight
title('Reconstructed Signal')
err = abs(s(Lo-Ls+1:Lo)-ss(Lo-Ls+1:Lo));
subplot(224)
plot(err),axis tight
title('Reconstruction Error')

figure(2),clf
pcolor(abs(C(floor((Lo-Ls+1)/dM):M,1:N/2)).')
shading interp
title('Spectrogram')


function y = mymod(x,N)
%mymod: like mod(x,N) or rem(x,N),if y == 0,substitute 0 with N
%   Example:
%           >>y=mymod([1:6],4)
%           y =
%                1     2     3     4     1     2
%Copyright(c) Xi'an Jiaotong University,Pengliang YANG,2010.5.1
y = mod(x,N);
y(find(y==0)) = N;

function [M,N] = mydivider(L)
%mydivider: Find dividers of an integer.  
%  	[M,N] = divider(L) find two integers M and N, such that 
%    M*N = L, M >= N, M and N as close as possible from sqrt(L).
%  	Example :
%       L = 258; [M,N] = mydivider(L)
%Copyright(c) Xi'an Jiaotong University, Pengliang YANG,2010.5.1
N = floor(sqrt(L));
flag = 1;
while flag == 1
    Nold = N;
    M = ceil(L/N);
    N = floor(L/M);
    flag = (N~=Nold);
end