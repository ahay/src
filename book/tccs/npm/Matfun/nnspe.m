%function [nt,tscale,fscale]=nnspe(data,t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale)
%
%
% Calling sequence-
%   [nt,tscale,fscale] = nnspe(data[, t0][, t1][, fres][, tres]
%                                  [, fw0][, fw1][,tw0][, tw1]
%                                  [, ifmethod][, normmethod][,nfilter][,lscale]) 
%
% Input-
%	  data        - 2-D input data(n,k) that specifies IMFs
%	  t0          - true start time
%	  t1          - true end time
%	  fres        - frequency resolution
%	  tres        - time resolution
%	  fw0         - minimum frequency
%	  fw1         - maximum frequency
%	  tw0         - minimum time for zooming
%	  tw1         - maximum time for zooming, if no zooming, tw0=t0, tw1=t1
%	  ifmethod    - nstantaneous frequency method, options: 'hilbert','hilbt','acos','zc','quad'
%	  normmethod  - normalization method, options: 'none','hilbert','spline','linear','pchip'
%	  nfilter     - number of points to filter by median
%   lscale      - when value=0 means linear axis scale,value= other integer means log axis scale  
%
% Output-
%	  nt          - 2-D matrix that specifies the spectrum
%	  tscale      - vector that specifies the time axis values
%	  fscale      - vector that specifies the frequency axis values
%
%
% The function NNSPE generates an improved HHT or zero-crossing spectrum of data (or LOG of data) data(n,k)
% in time-frequency (or time-log frequency) space, where n specifies the length of time series, and 
% k is the number of IMF components.
%
% NNSPE compute the Hilbert-Huang Spectrum of energy.
% Types of transform and normalization method are passed as arguments.
% All arguments except data, which is required, are optional.
% Optional arguments have default values assigned as follows:
%	t0=0; t1=max(size(data)-1); fres=400; tres=400; 
%	fw0= min(min(f)); fw1=max(max(f)); tw0=t0; tw1=t1;lscale=0;
% By default the regular scale is chosen. If the value of lscale is 
% non zero than the logarithmic scale is chosen.
% If no zooming in time is expected then tw0=t0 and tw1=t1.
%
% Non MATLAB Library routines used in the function are: FA.m, NSPPLOTE.m
%
% References:   
%  N. E Huang (2008),NCU hht class lecture 
%  6 Instantaneous Frequency.ppt
%
%
% code writer: Kenneth Arnold (NASA GSFC)	    	Summer, 2003 Initial
%	              (Based on work by Z. Shen (JHU) July, 2 1995.)
% code writer: Jelena Marshak (NASA GSFC)		March 17, 2004 Modified
%					      (Added the LOG scale option.) 
% code writer: Xianyao Chen, September, 2008
%       nnspe is update version of nnsp, some options for ifmethod and normmethod are added, 
%        default values for the input parameters are removed.
%        nnspe calculates the Hilbert-Huang spectrum of energy.
%        The countpart of nnspe is nnpsa, which calculate the spectrum of amplitude.
% code writer: S.C.su,August,2009
%       Add default input parameter values back(Norden E. Huang think those values are helpful for users) 
%       Add log scale calculation back (Norden E. Huang think those values are helpful for users) 
% footnote:S.C.Su 2009/08/8
%
% Many people separate Hilbert-Huang Tansform into 2 main parts
% 1. EMD(Empirical Mode Decomposition)  
% 2. HT (Hilbert Transform)
% This code is the integrated code for Hilbert Transform .
% After putting EMD or EEMD result into nnspe.m, 
%  the Hilbert-Huang Tansform(time-frequency spectrum) is finished.
%  
%  The structure of this code:
%  1. read input data from input arguments
%     1.1 read the IMF data,
%     1.2 give Normalze method and Instantaneous Frequency method  
%     1.3.the plotting parameters about the time-frequency spectrum
%  2.Give default values (about 1.2 and 1.3) 
%  3.call fa.m to get [F,A] 
%  4.call nspplote.m to plot the figure of time-frequency spectrum 
%  5.if no output arguments are passed,plot the figure on screen
%     
% Association: 
%  this function is the user interface of Hilbert Transform
% Concerned function: fa,nspplote
% 
%
% Ifmethod options:
%    'hilbert' :use Hilbert transform (function FAH )
%                  normalization of input data recommended but not required
%    'acos'    :use arccosine method (function FAACOS )
%		               normalization of input data required
%    'zc'      :use Generalized Zero-Crossing method (function FAZ )
%                  normalization of input data not recommended
%    'quad'    :use quadrature method (function calcqf),
%                  normalization of input data required
%    'hilbtm'  :use Hilbert transform, but hilbtm.m is used instead of the standard hilbert.m
%                  Gibbs phenomenom is been fixed ,when hilbtm.m is used for Hilbert-transform
%
% Normmethod options:
%    'none'	   :no normalization, recommended for 'zc' option
%    'spline'	 :spline normalization, recommended for 'hilbert' and 'acos' options
%               not recommended for Ensemble EMD method due to the possible overshot
%    'hilbert' (default) :Hilbert amplitude normalization, recommended when using Ensemble EMD method
%                this method is chosen as default because since 2008, EEMD will be mainly used
%    'linear'  :linear normalization, recommended for normalization when
%               using Ensemble EMD
%    'pchip'   :cubic hermite spline normalization, recommended for normalization when
%               using Ensemble EMD
%    'block'   :define the interval between zero-crossing points as a block
%               finding extreme value in the block,and normalize the whole value with the extreme    
%
% To display the output, use
%	   img(tscale,fscale,nt) or contour(tscale,fscale,nt)
%	or imagesc(tscale,fscale,nt);axis xy

function [nt,tscale,fscale]=nnspe(data,t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale)

%1. read input data from input arguments
%     1.1 read the IMF data,
%     1.2 give Normalze method and Instantaneous Frequency method  
%     1.3.the plotting parameters about the time-frequency spectrum


%2.Give default values (about 1.2 and 1.3)  
%----- Initialize the parameters---------------
%if user don't give, set it as empty 
if nargin<1
    error('nsp: data argument must be passed');
end
if max(size(data)) < 3
    error('nsp: data must have 3 or more points');
end
if max(size(data)) < 400
    disp('nsp: data must less than 400 points--should be careful ');
    disp('because default value of frequency resolution is 400 ');
end
if nargin<3
    t0 = [];
    t1 = [];
end
if nargin<5
    fres=[];
    tres=[];
end
if nargin<7
    fw0 = [];
    fw1 = [];
end
if nargin<9
    tw0=[];
    tw1=[];
end
if nargin<10
    ifmethod = [];
end
if nargin<11
    normmethod = [];
end
if nargin<12
    nfilter = [];
end
if nargin<13
    lscale = [];
end

%assign default values to those empty input parameters 
if isempty(t0)
    t0 = 0;
end
if isempty(t1)
    t1 = max(size(data))-1;
end
if isempty(fres)
    fres=400;
end
if isempty(tres)
    tres=400;
end
if isempty(tw0)
    tw0=t0;
end
if isempty(tw1)
    tw1=t1;
end
if isempty(lscale)
    lscale=0;
end

%
%----- check input
if (t1 < t0 | fw1 < fw0 | tw1 < tw0 | tw0 < t0 | tw1 > t1)
    error('check the region: t0 t1, tw0 tw1, and fw0 fw1.')
end

%----- Define the time interval
dt=(t1-t0)/(max(size(data))-1);

%  3.call fa.m to get [F,A] 
%----- Get frequency and amplitude
[freq,amp] = fa(data,dt,ifmethod,normmethod,nfilter);

%----- Get frequency bounds
if isempty(fw0)
    fw0=min(min(freq));
end
if isempty(fw1)
    fw1=max(max(freq));
end

%[fw0 fw1]

if fw0<0.
%    warning('nsp: negative frequency encountered, setting frequency window to 0');
    fw0=0;
end

%  4.call nspplote.m to plot the figure of time-frequency spectrum 
%----- Get the values to plot the spectrum
[nt,tscale,fscale] = nspplote(freq,amp,t0,t1,fres,tres,fw0,fw1,tw0,tw1,lscale);

%for clear appearance,clear the low-freqency if you want 
setclear=1;
 if setclear==1
  nt(1,:)=0;
 end

%  5.if no output arguments are passed,plot the figure on screen
%----- Plot the spectrum if no output arguments are passed
if nargout == 0
    imagesc(tscale,fscale,nt);axis xy;
end
