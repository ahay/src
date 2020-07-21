%function [f,a] = fa(data,dt,ifmethod,normmethod,nfilter)

% The function FA computes a frequency and amplitude of data(n,k), where 
% n specifies the length of time series, and k is the number of IMFs.
% The user has a choice to choose the instantaneous frequency and
% normalization methods. Nature of the arcosine method suggests not 
% to use the latter to process the residue component.
% First 2 arguments are required. If not passed the rest is set to have 
% default values: the frequency and amplitude is calculated using Hilbert
% method with spline normalization and no smoothing of data.
%
% Calling sequence-
% [f,a] = fa(data,dt[,ifmethod][,normmethod][,nfilter])
%
% Input-
%	  data		  - 2-D matrix of IMF components
%	  dt		    - time increment per point
%	  ifmethod	- method of determining an instantaneous frequency (and amplitude)
%	  normmethod	- normalization method
%	  nfilter		- number of points to use for filter
% Output-
%	  f		    - 2-D matrix f(n,k) that specifies frequency
%	  a		    - 2-D matrix a(n,k) that specifies amplitude
%
% Ifmethod options:
%    'hilbert' :use Hilbert transform (function FAhilbert )
%                  normalization of input data recommended but not required
%    'hilbtm'  :use Hilbert transform, but hilbtm.m is used instead of the standard hilbert.m(function FAimpHilbert )
%                  normalization of input data recommended but not required
%    'acos'    :use arccosine method (function FAacos )
%		               normalization of input data required
%    'zc'      :use Generalized Zero-Crossing method (function FAzc )
%                  normalization of input data not recommended
%    'quad'    :use quadrature method (function FAquadrature),
%                  normalization of input data required
%    'cosfor'  :use cosine-formula method (function FAcosfor),
%                  normalization of input data required 
%
%
% Normmethod options:
%    'none'	   :no normalization, recommended for 'zc' option 
%    'spline'	 :spline normalization, recommended for 'hilbert' and 'acos' options (function splinenormalize )
%               not recommended for Ensemble EMD method due to the possible overshot
%    'splineEP':spline normalization with several kind of end-process,(function splinenormalizeep )
%                  recommended for 'hilbert' and 'acos' options
%               not recommended for Ensemble EMD method due to the possible overshot
%    'hilbert' (default) :Hilbert amplitude normalization, recommended when using Ensemble EMD method (function hilbertnormalize )
%                this method is chosen as default because since 2008, EEMD will be mainly used
%    'linear'  :linear normalization, recommended for normalization when (function linearnormalize )
%               using Ensemble EMD
%    'pchip'   :cubic hermite spline normalization , recommended for normalization when (function pchipnormalize )
%               using Ensemble EMD
%    'block'   :block normalization, a experimental method ,not recommended to use (function blocknormalize )
%               user must use this method carefully 
%
%
% Non MATLAB Library routines used in the function are:
% 	FAACOS, FAH, FAZ,FAHT 
%	BLOCKNORMALIZE, SPLINENORMALIZE, HILBERTNORMALIZE, MEDIANFILTER.
%
%written by  
% Kenneth Arnold (NASA GSFC)	Summer 2003, Initial
% Karin Blank (NASA GSFC)       5/17/05, edited to add Quad
% Xianyao Chen (RCADA, FIO)   Sep. 20, 2008
%   1. add 'linear' and 'hermite' options for the normalization
%      to avoid the overshot of 'spline' method when using Ensemble EMD,
%   2. the default normalization method is given as 'hilbert' for the case
%      of 'ifmethod' is 'hilbert'
%   3. the extra necessary spline normalization when using 'quad' ifmethod
%      is changed from splinenormalize to hermitenormalize
%Sheng-Chung.Su 2009/09/10
%   1.rename and integrate all the mfile functions 
%   2.reset default Normalize---spline(splinenormalize.m
%     reset default IF       ---Improved Hilbert(FAimpHilbert.m)
%     reset default nfilter  ---5(5 point median filter)
%
%footnote: S.C.Su (2009/09/01)
%0.Initial the parameters and default settings
%1.Normalize the data by specific method
%2.Calculate the [f,a] by specific method
%3.use median-filter to process the I.F values 
%

function [f,a] = fa(data,dt,ifmethod,normmethod,nfilter)

%0.Initial the parameters and default settings
  %----- Define default parameters
  if nargin<3
      ifmethod = [];
  end
  if nargin<4
      normmethod = [];
  end
  if nargin<5
      nfilter=[];
  end
  
  if isempty(ifmethod)
      ifmethod = 'hilbtm';
  end
  if isempty(normmethod)
      if (isequal(ifmethod, 'hilbtm'))
          normmethod = 'spline';
      else
          normmethod = 'none';
      end
  end
  if isempty(nfilter)
      nfilter=5;
  end
  
  %----- Get the dimensions
  [npt,nIMF]=size(data);
  flipped=0;
  
  if npt<nIMF
      %----- Flip the data
      data=data';
      flipped=1;
      [npt, nIMF]=size(data);
  end
  
  odata = data;   %keep copy of original data KBB

%1.Normalize the data by specific method
%2009/10/14 Norden E Huang decided to operate normalize procedure 3 times in fa.m
  %----- Normalize data if requested    
  if isequal(normmethod, 'block')
      [data1,na1]=blocknormalize(data);
      [data2,na2]=blocknormalize(data1);
      [data3,na3]=blocknormalize(data2);
      data=data3;
      na=na1.*na2.*na3;
  elseif isequal(normmethod, 'spline')
      [data1,na1]=splinenormalize(data);
      [data2,na2]=splinenormalize(data1);
      [data3,na3]=splinenormalize(data2);
      data=data3;
      na=na1.*na2.*na3;
  elseif isequal(normmethod, 'linear')
      [data1,na1]=linearnormalize(data);
      [data2,na2]=linearnormalize(data1);
      [data3,na3]=linearnormalize(data2);
      data=data3;
      na=na1.*na2.*na3;
  elseif isequal(normmethod, 'pchip')
      [data1,na1]=pchipnormalize(data);  
      [data2,na2]=pchipnormalize(data1);  
      [data3,na3]=pchipnormalize(data2);    
      data=data3;
      na=na1.*na2.*na3;
  elseif isequal(normmethod, 'hilbert')
      [data1,na1]=hilbertnormalize(data);
      [data2,na2]=hilbertnormalize(data1);
      [data3,na3]=hilbertnormalize(data2);
      data=data3;
      na=na1.*na2.*na3;
  elseif isequal(normmethod, 'splineEP')
      [data1,na1]=splinenormalizeep(data);
      [data2,na2]=splinenormalizeep(data1);
      [data3,na3]=splinenormalizeep(data2);
      data=data3;
      na=na1.*na2.*na3;
      %or [data,na]=splinenormalizeep(data,'type');type=A,B,C
  elseif isequal(normmethod, 'none')
      na = [];
  else
      [data1,na1]=splinenormalize(data);
      [data2,na2]=splinenormalize(data1);
      [data3,na3]=splinenormalize(data2);
      data=data3;
      na=na1.*na2.*na3;
      disp ('fa: unknown normalization method,use default ');
  end % ... use spline normalize as default
  
%2.Calculate the [f,a] by specific method
  %----- Calculate the frequency and amplitude
  if (isequal(ifmethod, 'hilbert'))
      [f, a] = FAhilbert(data,dt);
  elseif (isequal(ifmethod, 'hilbtm'))
      [f, a] = FAimphilbert(data,dt);
  elseif (isequal(ifmethod,'acos'))
      [f, a] = FAacos(data,dt);
  elseif (isequal(ifmethod, 'zc'))
      [f,a]  = FAzc(data,dt);
  elseif (isequal(ifmethod, 'quad')) %added quad, KBB
      [f,a]  = FAquadrature(data,dt);
  elseif (isequal(ifmethod, 'cosfor')) %added cosine-formula, SCSu
      [f,a]  = FAcosfor(data,dt);    
  else
      [f, a] = FAimpHilbert(data,dt); 
      disp ('fa: unknown instantaneous frequency method,use default');
  end
  
  if ~isempty(na)
     %----- Throw away Hilbert etc. amplitude, use normalized amplitude
     a=na;
  end
  
%3.use median-filter to process the I.F values 
  %----- Filter the frequency if requested
  if nfilter>0
      for i=1:nIMF
          f(:,i)=medianfilter(f(:,i),nfilter);
      end
  end

  %----- Flip again if data was flipped at the beginning
  if flipped
      f=f';
      a=a';
  end
