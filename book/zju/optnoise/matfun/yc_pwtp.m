function [s2] = yc_pwtp(s1,sigma,p1,p2,type)
% yc_pwtp: plane-wave trace prediction (PWTP) in 2D
% predicting trace s2 from trace s1
% p1,p2: positions of s1,s2 (if p2>p1, forward prediction; else backward prediction)
%
% BY Yangkang Chen, Oct, 2017
%
% INPUT
% s1:       input trace
% sigma:    slope [n1xn2]
% p1,p2:    positions of s1,s2
% type:     discretization type (default: 1) compare 2: Plane-wave orthogonal polynomial transform for amplitude-preserving noise attenuation
%
% OUTPUT
% s2:       output trace
%
% DEMO
% test/test_yc_pwtp.m
%
% REFERENCE
% Chen et al., 2018, Plane-wave orthogonal polynomial transform for amplitude-preserving noise attenuation, GJI, 214, 2207-2223.
% Zhou et al., 2017, Spike-Like Blending Noise Attenuation Using Structural Low-Rank Decomposition, IEEE GRSL, 14, 1633-1637;
% Xie et al., 2018, Estimation of the anisotropy parameter $\delta$ and elastic parameter $C_{13}$ of organic-rich shale from laser ultrasonic technique (LUT) measurement, Geophysics, 83, C137?C152
% 
[n1,n2]=size(sigma);

if p1<0 || p2<0 || p1>n2 || p2>n2
    error('i1 and i2 are not correct !');
end

if length(s1)~=n1
    error('Wrong size for s1 !');
end

if nargin==4
    type=1;
end

if p2~=p1
    U=zeros(n1,n2);
    U(:,p1)=s1(:);
end
if p2>p1
    if type==1 %this discretization method is used in Zhou et al., Spike-Like Blending Noise Attenuation Using Structural Low-Rank Decomposition, IEEE GRSL. 
        for i2=p1:p2-1
            for i1=1:n1-1
                U(i1+1,i2+1)=(1-sigma(i1,i2))*U(i1+1,i2)+sigma(i1,i2)*U(i1,i2);
                % U(i1+1,i2+1)=(1+sigma(i1,i2))*U(i1+1,i2)-sigma(i1,i2)*U(i1,i2);
                %       U(i1+1,i2+1)=1/(1+sigma(i1,i2))*((1+sigma(i1,i2))*U(i1,i2)+(sigma(i1,i2)-1)*U(i1,i2+1)+(1-sigma(i1,i2))*U(i1+1,i2));
            end
        end
    end
else
    if p2<p1
        if type==1
            for i2=p1-1:-1:p2
                for i1=1:n1-1
                    U(i1+1,i2)=(1+sigma(i1,i2))*U(i1+1,i2+1)-sigma(i1,i2)*U(i1,i2+1);
                    % U(i1+1,i2+1)=(1+sigma(i1,i2))*U(i1+1,i2)-sigma(i1,i2)*U(i1,i2);
                    %       U(i1+1,i2+1)=1/(1+sigma(i1,i2))*((1+sigma(i1,i2))*U(i1,i2)+(sigma(i1,i2)-1)*U(i1,i2+1)+(1-sigma(i1,i2))*U(i1+1,i2));
                end
            end
        end
    else
        %p2==p1
        s2=s1;
    end
end
if p2~=p1
    s2=U(:,p2);
end

return

