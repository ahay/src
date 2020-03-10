function [ D1 ] =process_win(oper, param, D, twin, xwin)
% Processing seismic profile with t and x windows.
% D:            input data
% oper:         operator
% param:        parameters of operator
% twin:         t window length
% xwin:         x window length
%
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Aug, 2013
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%

tic

[sample,trace]=size(D);

k1=1;
k2=k1+twin-1;
l1=1;
l2=l1+xwin-1;
D1=zeros(size(D));

while k2<=sample
    if k1==1
        D_0=zeros(twin,trace);
        while l2<=trace
            if l1==1
                % key operations
                D0 = feval(oper,D(k1:k2,l1:l2),param);
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[ones(1,xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            else
                % key operations
                D0 = feval(oper,D(k1:k2,l1:l2),param);
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[(0:1:xwin/2-1)/(xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            end
            l1=l1+xwin/2;
            l2=l1+xwin-1;
            
        end
        % key operations
        D0 = feval(oper,D(k1:k2,l1:trace),param);               
        D_0(:,l1:trace)=D_0(:,l1:trace)+D0.*(ones(twin,1)*[(0:1:xwin/2-1)/(xwin/2),ones(1,trace-l1+1-xwin/2)]);
        D1(k1:k2,:)=D1(k1:k2,:)+D_0.*([ones(1,twin/2),(twin/2:-1:1)/(twin/2)]'*ones(1,trace));
    else
        D_0=zeros(twin,trace);
        l1=1;
        l2=l1+xwin-1;
        while l2<=trace
            if l1==1
                % key operations
                D0 = feval(oper,D(k1:k2,l1:l2),param);  
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[ones(1,xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            else
                % key operations               
                D0 = feval(oper,D(k1:k2,l1:l2),param);                 
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[(0:1:xwin/2-1)/(xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            end
            l1=l1+xwin/2;
            l2=l1+xwin-1;
        end

        % key operations               
        D0 = feval(oper,D(k1:k2,l1:trace),param);          
        
        D_0(:,l1:trace)=D_0(:,l1:trace)+D0.*(ones(twin,1)*[(0:1:xwin/2-1)/(xwin/2),ones(1,trace-l1+1-xwin/2)]);
        D1(k1:k2,:)=D1(k1:k2,:)+D_0.*([(0:1:twin/2-1)/(twin/2),(twin/2:-1:1)/(twin/2)]'*ones(1,trace));
    end
    k1=k1+twin/2;
    k2=k1+twin-1;
end

D_0=zeros(sample-k1+1,trace);
l1=1;
l2=l1+xwin-1;
while l2<=trace
    if l1==1
        % key operations               
        D0 = feval(oper,D(k1:sample,l1:l2),param); 
        D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(sample-k1+1,1)*[ones(1,xwin/2),(xwin/2:-1:1)/(xwin/2)]);
    else
        % key operations               
        D0 = feval(oper,D(k1:sample,l1:l2),param);         
        D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(sample-k1+1,1)*[(0:1:xwin/2-1)/(xwin/2),(xwin/2:-1:1)/(xwin/2)]);
    end
    l1=l1+xwin/2;
    l2=l1+xwin-1;
end

% key operations               
D0 = feval(oper,D(k1:sample,l1:trace),param); 
D_0(:,l1:trace)=D_0(:,l1:trace)+D0.*(ones(sample-k1+1,1)*[(0:1:xwin/2-1)/(xwin/2),ones(1,trace-l1+1-xwin/2)]);

D1(k1:sample,:)=D1(k1:sample,:)+D_0.*([(0:1:twin/2-1)/(twin/2),ones(1,sample-k1+1-twin/2)]'*ones(1,trace));

t=toc;

return












