%function [nt,tscale,fscale]=nspplote(f,a,t0,t1,fres,tres,fw0,fw1,tw0,tw1,lscale)
%
% Calling sequence-
% [nt,tscale,fscale]= nspplote(f,a[,t0][,t1][,fres][,tres][,fw0][,fw1][,tw0][,tw1][,lscale])
% Example, [nt,tscale,fscale]=nspplote(f,a,1,1000,500,300,0,1,1,1000).
% All arguments are required, WITHOUT default values.
%
% Input-
%	  f	          - 2-D matrix that specifies the frequency values
%	  a	          - 2-D matrix that specifies the amplitude values
%	  t0          - true start time
%	  t1          - true end time
%	  fres        - frequency resolution
%	  tres        - time resolution
%	  fw0         - minimum frequency
%	  fw1         - maximum frequency
%	  tw0         - minimum time for zooming
%	  tw1         - maximum time for zooming, if no zooming, tw0=t0, tw1=t1
%   lscale      - when value=0 means linear axis scale,value= other integer means log axis scale  
%        
% Output-
%	  nt          - 2-D matrix that specifies the spectrum
%	  tscale      - vector that specifies the time axis values
%	  fscale      - vector that specifies the frequency axis values
%
% NOTE
% NSPPLOTE.m compute the Hilbert-Huang Spectrum of energy.
% The function NSPPLOTE generates and plots the HHT Energy spectrum of data (or log data) 
% in time-frequency (or time-log frequency) space based on given frequency 
% f(n,k) and amplitude a(n,k), where n specifies the length of time series, and 
% k is the number of IMF components.
%
%
% References:   
%  N. E Huang (2008),NCU hht class lecture 
%  6 Instantaneous Frequency.ppt
%
% code writer:  Kenneth Arnold (NASA GSFC)		Summer, 2003 Initial
% code writer:  Jelena Marshak (NASA GSFC)		July 30, 2004 Modified
%                 added the LOG scale option
% code writer:  Xianyao Chen     September, 2008
%                 nspplote is update version of nspplot. 
%                 It is called by nnspe to compute the Hilbert-Huang spectrum of energy.
%                 The default input of each arguments are removed for better control of the program
% code writer: S.C.su,August,2009
%                 Add log scale calculation back (Norden E. Huang think those values are helpful for users) 
%                 Add check procedure for the matrix dimension,it's helpful when doing transfer function.
%                 When the dimension of the matrixs are stable,transfer function can calculate one by one divison.
% footnote:S.C.Su 2009/08/10         
%
%   This code is used when FM and AM is already been separated.
%   Code nspplote.m plots "Time-Frequency relationship" of those AM and FM signals.
%   The key concept is "Fill in those AM values in the corresponding FM positions"
%
%  The structure of this code:
% 1. read data and check input
% 2.start to construct the 2D grid
% 3.mapping the FM values into the grid 
% 4.when the FM values are correct,fill in the corresponding AM 
% 5. form the x,y axis coordinate values
% 6.plot the 2d figure (optional)
%

function [nt,tscale,fscale]=nspplote(f,a,t0,t1,fres,tres,fw0,fw1,tw0,tw1,lscale)

% 1. read data and check input 
%----- Check the input arguments
if nargin<11
    lscale=[];
end
if isempty(lscale)
    lscale=0;
end
if nargin < 10
    help nspplote
    error('The values of all parameters should be given.')
end
if size(f) ~= size(a)
    error('nspplot: frequency and amplitude matrix sizes differ');
end
%----- check input
if (t1 < t0 | fw1 < fw0 | tw1 < tw0 | tw0 < t0 | tw1 > t1)
    error('check the region: t0 t1, tw0 tw1, and fw0 fw1.')
end

%----- Check the frequency range
if abs(fw0-fw1)/fres<1e-10
    warning('nspplot: frequency is nearly constant; giving an artificial range of +/- 1');
    fw0=fw0-.5;
    fw1=fw1+.5;
end

%----- Get dimensions
[npt,nimf]=size(f);

%----- Set the log scale values if requested
if lscale ~= 0;
if((fw0<=0)|(fw1<=0))
    disp('==========================================================');
    disp('WARNING: min (or max) frequency range is less or equal to 0');
    disp('SUGGESTION: Change the frequency range and run again');
    disp('==========================================================');
end
    fw2=min(min(f));
    if (fw2 <= 0)
       warning('NSPPLOT:zero/negative frequency encountered, check the range');
       fw2
    end
    for i=1:nimf;
        % Apply LOG to a frequency
        for j=1:npt;
        if (f(j,i) > 0);
            f(j,i)=log(f(j,i));
        end
    end
    end
    fw0=log(fw0);
    fw1=log(fw1);
end

%----- Flip frequency and amplitude if necessary
if npt < nimf
    f=f';
    a=a';
    [npt,nimf]=size(f);
end

%2.start to construct the 2D grid
%----- Get local frequency and time
fw=fw1-fw0;
tw=tw1-tw0;

%----- Get time interval
dt=(t1-t0)/(npt-1);

%----- Construct the ploting matrix in time axis(time axis grid)
sidx=floor((tw0-t0)/dt)+1; 
eidx=ceil((tw1-t0)/dt)+1;

% The following  algorithm is needed when doing transfer function for time-freq spectrum   
% confirm the total grid number in time axis !
if ( t0 == tw0 & t1 == tw1)
   sidx=1;
   eidx=npt;   
end
  
nidx=eidx-sidx+1;

%-----check for the time axis grid number
%condition :  data point not enough ,user should  reduce time resolution! 
    
    if tres>nidx
        disp('data point in time axis not enough ! Please stop and reduce time resolution!')
        disp('if you continue doing this,the time resolution is changed ,the value becomes Npt ')  
        tres=nidx;
    end

%----- initial the matrix 
nt=zeros(tres,fres);

%3.mapping the FM values into the grid 
%----- Construct the ploting matrix in frequency axis(freq axis grid)
%P is the mapping position of frequency values into the freq axis grid
p=round((fres-1)*(f-fw0)/fw)+1;

%----- Map vector from point space into block space
%t is the mapping position of time values into the time axis grid
t=ceil((1:nidx)*tres/nidx);

%4.when the FM values are correct,fill in the corresponding AM 
%the most important part-"Fill in those AM values in the corresponding FM positions"
%because this calculate energy ,so AM square! 
% the following loop checking FM about its position
% when the position is correct ,put its corresponding energy in that position.
for x=sidx:eidx                                       %checking every FM values-loop A start  
    for imf=1:nimf                                     %checking the FM values for every IMF-loop B start 
        freqidx=p(x,imf);	                              % use P as the FM position index ,called freqidx
        if (freqidx >= 1 & freqidx <= fres)             %checking the position is 'inside' or 'outside' the grid 
            tx = t(x-sidx+1);                           %tx is the final position of frequency on time axis 
            nt(tx,freqidx)=nt(tx,freqidx)+a(x,imf)^2;   %put energy(AM*AM) in its position
        end
    end                                                %checking the FM values for every IMF-loop B end
end                                                   %checking every FM values-loop A end

%----- Define the output and mask out the neg ampl in case of log scale
nt=abs(nt);
if (lscale ~=0)
    for i=1:tres
    for j=1:fres
        if(nt(i,j) == 0.)
            nt(i,j)=-999.;
        else
            nt(i,j)=log(nt(i,j));
        end
    end
    end
end

% 5. form the x,y axis coordinate values
%form freq-axis grid value
fscale=linspace(fw0,fw1,fres)';
%form time-axis grid value
tscale=linspace(tw0,tw1,tres)';
nt=flipud(rot90((nt)));

%6.plot the 2d figure (optional)
%----- Plot if no output arguments are passed
if nargout == 0
    img(tscale,fscale,nt);
end
