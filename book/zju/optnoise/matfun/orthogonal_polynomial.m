function [p,alpha]=orthogonal_polynomial(x,n)
%x is offset coordination;  unit is km;
%p is orthogonal polynomial matrics
%x=0:0.05:3;%offset coordination;  unit is km
if (nargin==1),
  n=length(x);
end
p=zeros(n,length(x));%orthogonal polynomial matrics
alpha=zeros(n,n);%orthogonal polynomial coeffient
temp=0;
for j=1:n  %j present the order of polynomial
    for k=1:j %k present the kth item of the jth polynomial
        if j~=k 
        alpha(j,k)=sum((x.^(j-1)).*p(k,:));
        temp=temp+alpha(j,k)*p(k,:);
        end
     
    end
    alpha(j,j)=sqrt(sum(x.^(2*(j-1)))-sum(alpha(j,:).^2));
    p(j,:)=(x.^(j-1)-temp)/alpha(j,j);
    temp=0;
end

% for jj=1:n  %j present the order of polynomial
%     for k=1:jj %k present the kth item of the jth polynomial
%         if jj~=k 
%             temp=temp+sum((x.^(jj-1)).*p(k,:))/sum(p(k,:).^p(k,:))*p(k,:);
%         end
%      
%     end
%   %  alpha(j,j)=sqrt(sum(x.^(2*(j-1)))-sum(alpha(j,:).^2));
%     p(jj,:)=x.^(jj-1)-temp;
%     temp=0;
% end

 %p=[p0(x0)  p0(x1)   p0(x2) -----p0(xn);
 %  p1(x0)  p1(x1)   p1(x2) -----p1(xn);
 %  p2(x0)  p2(x1)   p2(x2) ------p2(xn);
 %   --------
 %  pn(x0)   pn(x1)  pn(x2)-----pn(xn)]
