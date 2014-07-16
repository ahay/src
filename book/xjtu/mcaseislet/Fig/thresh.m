function thresholding_op
close,clc,clear all

x=[-5:0.02:5];
thr=1.5;
normp=0.5;
y1=HardThresh(x,thr);
y2=SoftThresh(x,thr);
y3=pThresh(x,thr,normp);
y4=SteinThresh(x,thr);
y5=pexpThresh(x,thr,normp);
plot(x,y1,'k','linewidth',1)
hold on
plot(x,y2,'g','linewidth',1)
hold on 
plot(x,y3,'b','linewidth',1)
hold on
plot(x,y4,'r','linewidth',1)
hold on
plot(x,y5,'b--','linewidth',1)
grid on
title('Thresholding operators')
legend('Hard',...
    'Soft',...
    'pThresh,p=0.5',...
    'Stein',...
    'exponential,p=0.5',4)

function y=HardThresh(x,thr)
y=   x .* (abs(x) > thr);

function y=SoftThresh(x, thr)
y = x.* max(1.0-thr./(abs(x)+(x==0)),0.0);

function y=SteinThresh(x,thr)
y = x.* max(1.0-(thr./(abs(x)+(x==0))).^2,0.0);

function y=pThresh(x, thr, normp)
y = x.* max(1.0-(thr./(abs(x)+(x==0))).^(2-normp),0.0);

function y=pexpThresh(x, thr, normp)
y=x.*exp(-(thr./(abs(x)+(x==0))).^(2-normp));
