function operator()
clf; hold on;

k = linspace(-pi,pi,1800);

% ------------------------------------------------------------
% exact
s0 = k;
t0= s0./k;

% ------------------------------------------------------------
% Joe
c = 0.5.*(1+cos(k));
sj= s0.*c;
tj= sj./k;

% ------------------------------------------------------------
% second order
s2 = sin(k);
t2 = s2./k;
c2 = [0 0 0 -1/2 0 1/2 0 0 0];
% ------------------------------------------------------------
% fourth order
s4 = 4/3*sin(k) - 1/6*sin(2*k);
t4 = s4./k;
c4 = [0 0 1/12 -2/3 0 2/3 -1/12 0 0];
% ------------------------------------------------------------
% sixth order
s6 = 3/2*sin(k) - 3/10*sin(2*k) + 1/30*sin(3*k);
t6 = s6./k;
c6 = [0 -1/60 3/20 -3/4 0 3/4 -3/20 1/60 0];
% ------------------------------------------------------------
% eigth order
s8 = 8/5*sin(k) - 2/5*sin(2*k) + 8/105*sin(3*k) - 1/140*sin(4*k);
t8 = s8./k;
c8 = [1/140 -4/105 1/5 -8/5 0 8/5 -1/5 4/105 -1/140];

x=[-4:1:4];
subplot(3,1,1)
plot(x,c2,'g','LineWidth',2);hold on;
plot(x,c4,'b','LineWidth',2);
plot(x,c6,'c','LineWidth',2);
plot(x,c8,'m','LineWidth',2);
xlabel('x');
ylabel('Coefficient ');
%legend('2nd','4th','6th','8th','Location','SouthEastOutside')

subplot(3,1,2)
plot(k/pi,s0,'k','LineWidth',2); hold on;
%plot(k/pi,sj,'r','LineWidth',2);
plot(k/pi,s2,'g','LineWidth',2);
plot(k/pi,s4,'b','LineWidth',2);
plot(k/pi,s6,'c','LineWidth',2);
plot(k/pi,s8,'m','LineWidth',2);

axis([-1 +1 -pi +pi]);
xlabel('k [\pi]');
ylabel('Response ');
%legend('exact','2nd','4th','6th','8th','Location','SouthEastOutside')

subplot(3,1,3)
plot(k/pi,t0,'k','LineWidth',2); hold on;
%plot(k/pi,tj,'r','LineWidth',2); 
plot(k/pi,t2,'g','LineWidth',2);
plot(k/pi,t4,'b','LineWidth',2);
plot(k/pi,t6,'c','LineWidth',2);
plot(k/pi,t8,'m','LineWidth',2);

axis([-1 +1 0 1.1]);
xlabel('k [\pi]');
ylabel('Weight');
legend('exact','2nd','4th','6th','8th','Location','SouthOutside','orientation','horizontal')

