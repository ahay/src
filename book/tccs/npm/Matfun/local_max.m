%function [local_max_x, local_max_y] = local_max(data_x, data_y)
%
% The function LOCAL_MAX returns local max point, their coordinates and the number of values.
% Identify local maxima points of physical signal
% Note: in the code initially the local variable flag = 0.
% It changes its value: flag = 0 when a slope is rising,
% and flag = 1 when a slope is descending.
%
% Calling sequence-
% [local_max_x,local_max_y] = local_max(data_x, data_y)
%
% Input-
%	data_x		- input vector of coordinates
%	data_y		- input vector of corresponding values
% Output-
%	local_max_x	- vector that specifies the coordinates 
%			          of max values in the order found
%	local_max_y	- corresponding max values
% 
% Jelena Marshak (NASA GSFC)    May 8, 2004 Modified
%          (added variable flag initialization to 0)
% Footnote:S.C.Su

function [local_max_x, local_max_y] = local_max(data_x, data_y)


flag = 0; %flag = 0 when slope rising, 1 when slope descending

local_max_y = [];
local_max_x = [];
num_extrema = 0;
Y = [];
X = [];

%----- Find the extrema values
%loop start ---check every value in the series 
for i=2:length(data_y)-1;
        
      % when slope descending (flag=1) ,the value of i position =left or right hand side value 
        if(((data_y(i) == data_y(i+1)) | (data_y(i) == data_y(i-1))) & (flag == 1))
           %mark this position
            Y = [Y, data_y(i)];
            X = [X, data_x(i)];
        end
      
      %first time start here, because initial setting flag=0
      % when slope rising (flag=0) ,the value of i position >left hand side value  
        if((data_y(i) > data_y(i-1)) & (flag == 0)  & ~(data_y(i) == data_y(i-1)))
           %set slope descending (flag=1),mark the position
            flag = 1;
            Y = data_y(i);
            X = data_x(i);
        end
      
      % when slope descending (flag=1) ,the value of i position <right hand side value  
        if((data_y(i) < data_y(i+1)) & (flag == 1))
            %set slope rising 
            flag = 0;      
            Y = [];
            X = [];
        end
        
      % when slope descending (flag=1) ,the value of i position >right hand side value  
        if((data_y(i) > data_y(i+1) & (flag == 1)))
            %set slope rising,mark the position
            flag = 0;
            local_max_y = [local_max_y, Y];
            local_max_x = [local_max_x, X];
            Y = [];
            X = [];
            num_extrema = num_extrema+1;
        end
        
%loop end ---check every value in the series    
end


