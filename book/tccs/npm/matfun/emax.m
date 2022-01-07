%function [ma, ta]=emax(x)
%
%
% Input-
%	        x	- vector of input data x(n)
% Output-
%	       ma	- vector that specifies max points
%	       ta	- vector that specifies the coordinates of max points
%
%NOTE:
%      give a 1D data-x   
%      The function EMAX returns the  max points and their coordinates.
%
%  Reference:
% Association: NFAM, NFAMSM, NFAMSMMULTI  
%     
%
%  code writer:Norden Huang (NASA GSFC)	November 11, 1999
% footnote:S.C.Su 2009/05/14

function [ma, ta]=emax(x)
 
% Calling sequence-
% [ma, ta]=emax(x)
%

%----- Get dimensions
n=length(x);

%----- Initialize
n_x=1;
ma=0;
ta=0;

%----- Extract the set of max points and their coordinates
for i=2:n-1
   if (x(i-1)<x(i))&(x(i)>=x(i+1))
      ma=[ma x(i)];
      ta=[ta i];
      n_x=n_x+1;

   end
end

%----- Define the output
ma=ma(2:n_x);
ta=ta(2:n_x);

ma=ma';
ta=ta';
