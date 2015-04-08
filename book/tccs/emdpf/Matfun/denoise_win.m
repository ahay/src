function [ D1 ] =denoise_win(D,flow,fhigh,mode, dt, N, lf, mu, twin, xwin)
%  denoise with t-win and x-win
%  mode=1,2,3 or other -> fxdecon, fxemd,fxemdpf
%
%  Author      : Yangkang Chen
%                Texas Consortium of Computational Seismology
%                Jackson School of Geosciences
%                The University of Texas at Austin
%         
%  Date        : Sep, 2013
%
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
%  
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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
                if mode==1
                    D0=fx_decon(D(k1:k2,l1:l2),dt,lf,mu,flow,fhigh);
                else 	if mode==2
                    D0=fx_emd(D(k1:k2,l1:l2),flow,fhigh,dt,N);
			else
                    D0=fx_emdpf(D(k1:k2,l1:l2),flow,fhigh,dt,N, lf, mu);
			end
                end
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[ones(1,xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            else
                if mode==1
                    D0=fx_decon(D(k1:k2,l1:l2),dt,lf,mu,flow,fhigh);
                else 	if mode==2
                    D0=fx_emd(D(k1:k2,l1:l2),flow,fhigh,dt,N);	
			else		
                    D0=fx_emdpf(D(k1:k2,l1:l2),flow,fhigh,dt,N, lf, mu);
			end
                end
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[(1:1:xwin/2)/(xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            end
            l1=l1+xwin/2;
            l2=l1+xwin-1;
            
        end
        if mode==1
            D0=fx_decon(D(k1:k2,l1:trace),dt,lf,mu,flow,fhigh);
        else 	if mode==2
            D0=fx_emd(D(k1:k2,l1:trace),flow,fhigh,dt,N);
		else 
            D0=fx_emdpf(D(k1:k2,l1:trace),flow,fhigh,dt,N, lf, mu);
		end
        end
        D_0(:,l1:trace)=D_0(:,l1:trace)+D0.*(ones(twin,1)*[(1:1:xwin/2)/(xwin/2),ones(1,trace-l1+1-xwin/2)]);
        D1(k1:k2,:)=D1(k1:k2,:)+D_0.*([ones(1,twin/2),(twin/2:-1:1)/(twin/2)]'*ones(1,trace));
    else
        D_0=zeros(twin,trace);
        l1=1;
        l2=l1+xwin-1;
        while l2<=trace
            if l1==1
                if mode==1
                    D0=fx_decon(D(k1:k2,l1:l2),dt,lf,mu,flow,fhigh);
                else 	if mode==2
                    	D0=fx_emd(D(k1:k2,l1:l2),flow,fhigh,dt,N);
			else
                    	D0=fx_emdpf(D(k1:k2,l1:l2),flow,fhigh,dt,N, lf, mu);
			end
                end
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[ones(1,xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            else
                if mode==1
                    D0=fx_decon(D(k1:k2,l1:l2),dt,lf,mu,flow,fhigh);
                else 	if mode==2
                    	D0=fx_emd(D(k1:k2,l1:l2),flow,fhigh,dt,N);
			else
                    	D0=fx_emdpf(D(k1:k2,l1:l2),flow,fhigh,dt,N, lf, mu);
			end
                end
                D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(twin,1)*[(1:1:xwin/2)/(xwin/2),(xwin/2:-1:1)/(xwin/2)]);
            end
            l1=l1+xwin/2;
            l2=l1+xwin-1;
        end
        if mode==1
            D0=fx_decon(D(k1:k2,l1:trace),dt,lf,mu,flow,fhigh);
        else	if mode==2
            	D0=fx_emd(D(k1:k2,l1:trace),flow,fhigh,dt,N);
		else 
            	D0=fx_emdpf(D(k1:k2,l1:trace),flow,fhigh,dt,N, lf, mu);
		end
        end
        D_0(:,l1:trace)=D_0(:,l1:trace)+D0.*(ones(twin,1)*[(1:1:xwin/2)/(xwin/2),ones(1,trace-l1+1-xwin/2)]);
        D1(k1:k2,:)=D1(k1:k2,:)+D_0.*([(1:1:twin/2)/(twin/2),(twin/2:-1:1)/(twin/2)]'*ones(1,trace));
    end
    k1=k1+twin/2;
    k2=k1+twin-1;
end

D_0=zeros(sample-k1+1,trace);
l1=1;
l2=l1+xwin-1;
while l2<=trace
    if l1==1
        if mode==1
            D0=fx_decon(D(k1:sample,l1:l2),dt,lf,mu,flow,fhigh);
        else 	if mode==2
            	D0=fx_emd(D(k1:sample,l1:l2),flow,fhigh,dt,N);
		else
            	D0=fx_emdpf(D(k1:sample,l1:l2),flow,fhigh,dt,N, lf, mu);
		end		
        end
        D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(sample-k1+1,1)*[ones(1,xwin/2),(xwin/2:-1:1)/(xwin/2)]);
    else
        if mode==1
            D0=fx_decon(D(k1:sample,l1:l2),dt,lf,mu,flow,fhigh);
        else 	if mode==2
            	D0=fx_emd(D(k1:sample,l1:l2),flow,fhigh,dt,N);
		else
            	D0=fx_emdpf(D(k1:sample,l1:l2),flow,fhigh,dt,N, lf, mu);
		end
        end
        D_0(:,l1:l2)=D_0(:,l1:l2)+D0.*(ones(sample-k1+1,1)*[(1:1:xwin/2)/(xwin/2),(xwin/2:-1:1)/(xwin/2)]);
    end
    l1=l1+xwin/2;
    l2=l1+xwin-1;
end
if mode==1
    D0=fx_decon(D(k1:sample,l1:trace),dt,lf,mu,flow,fhigh);
else 	if mode==2
    	D0=fx_emd(D(k1:sample,l1:trace),flow,fhigh,dt,N);
	else 
    	D0=fx_emdpf(D(k1:sample,l1:trace),flow,fhigh,dt,N, lf, mu);
	end
end

D_0(:,l1:trace)=D_0(:,l1:trace)+D0.*(ones(sample-k1+1,1)*[(1:1:xwin/2)/(xwin/2),ones(1,trace-l1+1-xwin/2)]);
D1(k1:sample,:)=D1(k1:sample,:)+D_0.*([(1:1:twin/2)/(twin/2),ones(1,sample-k1+1-twin/2)]'*ones(1,trace));

t=toc

return













