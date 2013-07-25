clear;clc;close all

%% Description
% In this exercise, we solve Ax = b using 2 different methods:
% Steepest decent
% Conjugate direction

%% Initialization
A = [3. 1; 1 3];
b = [1 2]';
at = A\b;
a0 = [40 -90]';

%% Compute the objective function values
x1 = -100:1:100;
x2 = -100:1:100;
[X1 X2] = meshgrid(x1,x2);

obj = zeros(length(x1), length(x2));
for i = 1:length(x1)
    for j = 1:length(x2)
        r = A * [x1(i) x2(j)]' - b; 
        obj(i,j) = sum(r.^2);
    end
end

%% Steepest decent iteration
maxit = 20;
a1 = zeros(2,maxit); r1 = a1; g1 = a1;
alpha1 = zeros(1,maxit);
a1(:,1) = a0;
r1(:,1) = A  * a1(:,1) - b;

for iter=1:maxit
    g1(:,iter) = A' * r1(:,iter); 
    dd = A * g1(:,iter);
    alpha1(iter) = -r1(:,iter)'*dd / (dd'*dd);
    a1(:,iter+1) = a1(:,iter) + alpha1(iter)*g1 (:,iter);
    r1(:,iter+1) = r1(:,iter) + alpha1(iter)*dd;
end
obj1 = sum(r1.^2,1);

%% Conjugate direction iteration
maxit = 21;
a3 = zeros(2,maxit); r3 = a3; g3 = a3;
alpha3 = zeros(1,maxit); beta3 = alpha3;
a3(:,1) = a0;
for iter=1:maxit-1
    r3(:,iter) = -A  * a3(:,iter) + b;
    g3(:,iter) = A' * r3(:,iter); 
    dd = A * g3(:,iter);
    if (iter==1)
        alpha3(iter)= r3(:,iter)'*dd / (dd'*dd);
        beta3(iter) = 0.;
        s3  = 0.;
        ss3 = 0.;
    else
        gdg = dd'*dd;
        sds = ss3'*ss3;
        gds = dd'* ss3;
        determ = gdg * sds * max(1.d0 - (gds/gdg)*(gds/sds), 1.d-12);
        gdr = r3(:,iter)'* dd;
        sdr = r3(:,iter)'* ss3;
        alpha3(iter) = ( sds * gdr - gds * sdr ) / determ;
        beta3 (iter) = (-gds * gdr + gdg * sdr ) / determ;
    end
    s3  = alpha3(iter) * g3(:,iter) + beta3(iter) * s3;            % update solution step
    ss3 = alpha3(iter) * dd + beta3(iter) * ss3;            % update residual step
    a3(:,iter+1)  = a3(:,iter) +  s3;                          % update solution
    r3(:,iter+1)  = r3(:,iter) + ss3;                          % update residual
end
obj3 = sum(r3.^2,1);

%% plot path
h = figure(1);
v = interp(obj1,4);
contour(X1,X2,obj,v,'--k','LineWidth',0.5)
hold on

plot(a1(2,:),a1(1,:),'-.ks','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',2)
plot(a3(2,:),a3(1,:),'-gd','LineWidth',1,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',2)     
plot(at(2),at(1),'p','MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',10)
title('Search paths')
xlabel('m1');ylabel('m2')
print(h,'-deps','junk_ml.eps')