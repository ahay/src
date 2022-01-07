%function [ndata,a]=splinenormalize(data)
 
% The function SPLINENORMALIZE normalizes with splines,
% the data(n,m) where n specifies the number of time points, 
% and m is the number of IMF components.
% The normalization is carried out with both maxima and minima; hence emax
% is taken as the absolute sign.
%
% Calling sequence-
% [ndata,a]=splinenormalize(data)
%
% Input-
%	data	- 2-D matrix data(n,m) that specifies the IMF components 
% Output-
%	ndata	- normalized data (FM)
%	a	    - splined envelope(AM)
%
% Used by-
%	FA

%written by 
% K. Arnold (for NASA GSFC)	    Jan. 28 2004 Modified
%				(fixed a bug where a "random" point 
%				could get thrown somewhere in the middle
%				of the dataset and the endpoints 
%				were not controlled).
%  Modified Karin Blank March 22, 2005
%footnote: S.C.Su (2009/09/02)
%
% There are two loops in this code ,it's dealing with IMF by IMF .
% 0. set the default value about the  end-process methodtype
%  1.read the data,check input matrix
%  2.Do the normalization  for an IMF ---loop A start   
%        3.find maximun for the asbolute value 
%        4.add check for trend--loop B start  
%           5.Do end-process for bouble sides 
%           6.spline interpolation for the envelope-AM
%           7.normalize the IMF,AM been devided out
%        4. add check for trend--loop B end    
%  2.add stop algorithm-for the normalize procedure-loop A end  
%
function [ndata,a]=splinenormalize(data)

% 0. set the default value about the  end-process methodtype
%----- Get the dimension
[npt,ncol] = size(data);

%----- Initialize and flip data if needed 
flipped=0;
if (ncol > npt)
    data=data';
    [npt,ncol] = size(data);
    flipped=1;
end

te=(1:npt)';  %a series of total number

%2.Do the normalization  for an IMF ---loop A start 
%----- Process each IMF component
for c=1:ncol   % loop-a start here

%3.find maximun for the asbolute value
    %----- Extract the set of max points and their coordinates
    [mx, tx]=emax(abs(data(:,c)));
%4. add check for trend--loop B start    
    %----- Fix the end to prevent wide swaying in spline 
    %----- by assigning the te(1) and te(n) 
    %----- the same values as the first and last tx and mx.
     nExtrema=length(mx);
    if nExtrema > 1
%5.Do end-process for bouble sides    
        tx=[1;tx;npt];
        %----- Fix the ends at the same as the next point
        mx=[mx(1);mx;mx(nExtrema)];  
          
%6.spline interpolation for the envelope-AM 
        %find spline envelope as AM 
        a(:,c)=spline(tx,mx,te);
        %***linear interpolation backup in case spline error. Added by KBB.
        %Experimental
        if(any(a(:, c) < data(:,c)))
            a2 = interp1(tx, mx, te);
            
            mask_sp = ones(length(a), 1);
            
            for(i=1:length(tx)-1)
                if(any(a(tx(i):tx(i+1),c) < data(tx(i):tx(i+1),c)))
                    mask_sp(tx(i):tx(i+1)) = 0;
                end
            end
            
    %        mask_sp = (a(:,c) >= data(:, c)); 
        
            a(:, c) = (a(:,c) .* mask_sp) + (a2 .* ~mask_sp);
        end      
        %end KBB
        
%7.normalize the IMF,AM been devided out          
        %----- %  Normalize the data by splined envelope
        ndata(:,c)=data(:,c)./a(:,c);
        

    else
        %----- Leave data unchanged
        a(:,c)=ones(npt,1);
        ndata(:,c)=0*ones(npt,1);
        disp('Watch! there is a trend in the input Matrix')
    end
%4. add check for trend--loop B end     
end
%2.Do the normalization  for an IMF ---loop A end

%----- Flip data back if needed
if (flipped)
    ndata=ndata';
    a=a';
end
