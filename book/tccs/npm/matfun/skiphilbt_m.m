%function [bx, ex] = skiphilbt_m(imf)
%
%  INPUT:
%         IMF: an IMF 1D matrix
%  OUTPUT:
%        bx: when the value=1,need no modification for Gibbs phenomenon in the begginning of this imf
%        ex: when the value=1,need no modification for Gibbs phenomenon in the endding of this imf
% 
%    NOTE:
%     1.check if intermittency is at endpoints and skip appending waves if so 
%     2.This code performs starting and endding segment checking
%       when first or last value of imf are 0,0,0 ,this means no need for modification.
%       [0,0,0]means value=0,slope=0,curvature=0 ,the Gibbs effect is already fixed.
%     3.before using hilbtm.m to perform a extending segment to the original imf
%       in the beginning and endding part, one should use  skiphilbt_m to check the necessarity.
%%
% References:
%
% code writer: Karin Blank (karin.blank@nasa.gov)
% footnote:S.C.Su 2009/06/20
%
% Association: hilbtm.m
% this function is a check procedure before hilbtm.m
%hilbtm,m is an improved edition used for calculating the Hilbert-Transform. 
% matlab original function disadvantage is been solved.
%
% Concerned function: hilbtm.m
%                     above mentioned m file must be put together
%
%

function [bx, ex] = skiphilbt_m(imf)

%initiate with 0
bx = 0;
ex = 0;

%Are 3 values in a row zero?
%the head
if(imf(1) == 0 & imf(2) == 0 & imf(3) == 0)
    bx = 1;
end

%the end
if(imf(end) == 0 & imf(end-1) == 0 & imf(end-2) == 0)
    ex = 1;
end

return;