function [va vb vc]=vsign(a,b,c,p,q,r)

va=a;vb=b;vc=c;
if a*p'<0  
    va=-a;
end

if b*q'<=0
    vb=-b;
end

if c*r'<=0
    vc=-c;
end

