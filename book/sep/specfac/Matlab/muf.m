function muf(s)
x=-3:.001:3;

f1=abs(x+sqrt(s)).^((sqrt(s)-1)/(2*sqrt(s)));
f2=abs(x-sqrt(s)).^((sqrt(s)+1)/(2*sqrt(s)));

f=f1.*f2;

clf
hold on
plot(x,f1,'b--');
plot(x,f2,'k--');
plot(x,f ,'k');
xlabel('x');
ylabel('f(x)');

grid on;
