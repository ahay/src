%Version 1.00 stable
%**************************************************************************
%*************************** Wavelet Transform ****************************
%**************************************************************************
%-------------------------------Copyright----------------------------------
%
% Author: Dmytro Iatsenko
% Information about these codes (e.g. links to the Video Instructions),
% as well as other MatLab programs and many more can be found at
% http://www.physics.lancs.ac.uk/research/nbmphysics/diats/tfr
% 
% Related articles:
% [1] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part I: Overview, standards of use, related issues and algorithms."
% {preprint:arXiv:1310.7215}
% [2] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part II: Resolution, reconstruction and concentration."
% {preprint:arXiv:1310.7274}
%
%------------------------------Documentation-------------------------------
%
% [WT,freq,Optional:wopt]=wt(sig,fs,Optional:'PropertyName',PropertyValue)
%
% INPUT:
% sig - signal for which to calculate WT
% fs  - sampling frequency of the signal
% 
% Properties: ({{...}} denotes default)
% 'fmin':value (default is the minimal frequency for which at least one WT
%               coefficient is determined up to a specified relative
%               accuracy [RelTol] with respect to boundary errors)
%            minimal frequency for which to calculate WT
% 'fmax':value (default = fs/2, i.e. the Nyquist frequency)
%            maximal frequency for which to calculate WT
% 'nv':{{'auto'}}|'auto-NB'|value
%            "number of voices", which determines frequency discretization,
%            so that the next frequency equals previous one multiplied on
%            [2^(1/nv)]; when set to 'auto-NB' (e.g. 'auto-20') determines
%            [nv] automatically as described in [1], so that 1/nv equals
%            1/NB of the logarithmic frequency region containing 50% of the
%            wavelet function; the default 'auto' is equivalent to 'auto-10'
% 'Wavelet':{{'Lognorm'}}|'Morlet'|'Bump'|'Morse-a'|{@(xi)fwt(xi),[xi1,xi2],@(t)twf(t),[t1,t2]}
%            wavelet used in WT calculation, for a list of all supported
%            names and their properties see Appendix E in [1]. However, you
%            can use any wavelet by specifying its frequency domain form,
%            i.e. wavelet FT (defined by function [fwt], the argument of
%            which is cyclic frequency) and/or time-domain form [twf] (the
%            argument is time) together with corresponding full supports;
%            only [fwt] or [twf] is enough, so if one of them is not known
%            just put empty field [] for it and its support (but it is
%            better to specify both [fwt] and [twf] when available). Thus,
%            the Lognormal wavelet with [f0=1] can be alternatively
%            defined as {@(xi)exp(-(6^2/2)*(log(xi)).^2),[0,Inf],[],[]}.
% 'f0':value (default = 1)
%            wavelet resolution parameter, which determines the tradeoff
%            between the time and frequency resolutions: the higher it is,
%            the closer in frequency components can be resolved in WT,
%            but the slower time-variations, e.g. amplitude/frequency
%            modulation, can be reliably represented; for the way it is
%            introduced for each wavelet see Appendix E in [1], while if the
%            wavelet is user-defined in terms of its function in frequency
%            and/or time (see 'Wavelet' property), then 'f0' obviously does
%            not influence anything.
% 'Preprocess':{{'on'}}|'off'
%            perform or not an initial signal preprocessing, which consists
%            of subtracting 3rd order polynomial fit and then bandpassing
%            the signal in the band of interest [fmin,fmax], for which WT
%            is calculated
% 'Padding':{{'predictive'}}|0|'symmetric'|'none'|'periodic'|value|{padleft,padright}
%            what padding to use when computing transform (for a list of
%            all paddings and their effects see [1]); when set to some
%            numeric value, pads with that values; most useful are the
%            zero-padding, for which boundary errors are well-determined,
%            and 'predictive' padding (default), for which they are
%            most reduced, while other choices are not very useful;
%            alternatively, one can specify it as {padleft,padright}, where
%            [padleft] and [padright] can be any of those (except 'none'),
%            e.g. {0,'symmetric'} will pad signal to the left with zeros,
%            and to the right by symmetric reflection; additionally,
%            [padleft] and [padright] can be vectors to pad with (if there
%            needed more values than padleft/padright lengths,
%            pads additionally with zeros), e.g. if you want to consider
%            only some part of the signal, from indexes n1 to n2, but do
%            not introduce additional boundary effects to it, then use:
%            wt(signal(n1:n2),fs,'Padding',{signal(1:n1-1),signal(n2+1:end)});
%            the same applies if you can predict signal out of the time-limits
% 'RelTol':value (default = 0.01) (in [1] commonly referred as \epsilon)
%            relative tolerance (e.g. 0.01 means 1%), which specifies cone
%            of influence for WT (i.e. range of WT coefficients which are
%            determined up to this accuracy in respect of boundary errors);
%            determines also the minimal number of values to pad signal
%            with, so that relative contribution of effects of implicit
%            periodic signal continuation due to convolution in frequency
%            domain is smaller [RelTol], see [1] for details.
% 'Plot':{{'off'}} - do not plot anything
%          'amp'   - plots WT amplitude (i.e. its absolute value) together
%                    with line denoting the cone of influence
%          'amp+'  - additionally shows time-averaged WT amplitude
%          'amp++' - additionally shows 95% range of WT amplitude
%          'pow','pow+','pow++' - the same but for the WT power (i.e. its 
%                                 squared modulus)
%           IMPORTANT: to avoid plotting huge data (in which case it might
%           be very slow to render and modify the figure, and MatLab can
%           even crash), the plotted WT is resampled to have no more than
%           few data points displayed per pixel (for the current screen
%           resolution). Unless 'Display' is 'off', it will always notify if
%           the WT size exceeds the number of pixels on the plot and the
%           resamling is therefore performed. Note, that in this case considerable
%           zooming in of the parts of displayed plots might not show the full
%           structure of the original WT. If one wants to investigate the
%           resultant plot in fine details (and not only see how it looks),
%           then the original, full WT might be displayed without resampling
%           by adding 'wr' to the end of this option, e.g. 'ampwr' or 'pow++wr'.
% 'Display':{{'on'}} - displays all relevant information about progress etc.
%           'notify' - displays only information if something went wrong
%            'off'   - does not display anything
% 'CutEdges':{{'off'}}|'on'
%           determine should WT coefficients be set to NaNs out of the cone
%           of influence (determined for the current 'RelTol'); if you wish
%           to analyze only WT within the cone of influence, which is
%           recommended if you want to estimate e.g. only the time-averaged
%           quantities, set it to 'on'.
%
%
% OUTPUT:
% WT   - wavelet transform of the signal [sig],
%        with rows corresponding to frequencies and columns - to time;
%        represents FNxL matrix, where FN is the number of frequencies
%        and L is the length of the signal in samples
% freq - frequencies corresponding to rows of WT
% wopt - structure with all parameters of the wavelet and simulation
%
%-------------------------------Examples-----------------------------------
%
% [WT,freq]=wt(sig,100,'fmin',0.1,'fmax',10,'f0',2)
% given a signal [sig] sampled at 100 Hz, returns its WT [WT] based on a
% Lognormal wavelet (default) with [f0=2], calculated for frequencies from
% 0.1 to 10 Hz, returned in [freq].
%
%-----------------------Additional possibilities---------------------------
%
% One can also pass the structure with the properties as a third input
% argument instead of specifying them by pairs, e.g.
% opts=struct; opts.f0=2; opts.fmax=12; opts.Padding='predictive';
% [WT,freq]=wt(sig,fs,opts);
% You can also add further parameters in the usual way, e.g.
% [WT,freq]=wt(sig,fs,opts,'Display','off');
% If you add a parameter contained in [opts], it will be ovewritten, e.g.
% [WT,freq]=wt(sig,fs,opts,'fmax',10);
% will change 'fmax' property from 12 specified in [opts.fmax] to 10.
%
% When one needs to calculate WT using the same wavelet and simulation
% parameters, but for different signals, it is a good idea to do it like
% [WT1,freq1,wopt]=wt(sig1,fs1,...(parameters));
% [WT2,freq2]=wt(sig2,fs2,wopt); [WT3,freq3]=wt(sig3,fs3,wopt); ...
% This will also avoid recalculating wavelet parameters each time, thus
% slightly speeding up the computations.
%
% The same can be done if you calculated WT of a signal but want to obtain
% it using some different parameter (e.g. Padding), for example
% [WT1,freq1,wopt]=wt(sig,fs,'Padding',0,...(other parameters));
% [WT2,freq2]=wt(sig,fs,wopt,'Padding','predictive');
% Unless you overwrite 'f0', 'Wavelet' or 'RelTol' properties, the wavelet
% parameters will not be recalculated for the second time.
%
%--------------------------------------------------------------------------


function [WT,freq,varargout] = wt(signal,fs,varargin)

L=length(signal); signal=signal(:);
p=1; %WT normalization

%Default parameters
Wavelet='Lognorm'; f0=1;
fmin=[]; fmax=fs/2;
nv='auto';
PadMode='predictive';
RelTol=0.01;
Preprocess='on';
DispMode='on';
PlotMode='off';
CutEdges='off';
%Update if user defined
vst=1; recflag=1;
if nargin>2 && isstruct(varargin{1})
    opts=varargin{1}; vst=2;
    if isfield(opts,'Wavelet'), Wavelet=opts.Wavelet; end
    if isfield(opts,'f0'), f0=opts.f0; end
    if isfield(opts,'fmin'), fmin=opts.fmin; end
    if isfield(opts,'fmax'), fmax=opts.fmax; end
    if isfield(opts,'nv'), nv=opts.nv; end
    if isfield(opts,'Padding'), PadMode=opts.Padding; end
    if isfield(opts,'RelTol'), RelTol=opts.RelTol; end
    if isfield(opts,'Preprocess'), Preprocess=opts.Preprocess; end
    if isfield(opts,'Plot'), PlotMode=opts.Plot; end
    if isfield(opts,'Display'), DispMode=opts.Display; end
    if isfield(opts,'CutEdges'), CutEdges=opts.CutEdges; end
    if isfield(opts,'wp')
        recflag=0;
        for vn=vst:nargin-2
            if strcmpi(varargin{vn},'Wavelet'), recflag=1; end
            if strcmpi(varargin{vn},'f0'), recflag=1; end
            if strcmpi(varargin{vn},'RelTol'), recflag=1; end
        end
        wp=opts.wp;
    end
    if isfield(opts,'nvsim') && recflag==1, nv=opts.nvsim; end
end
for vn=vst:2:nargin-2
    if strcmpi(varargin{vn},'Wavelet'), Wavelet=varargin{vn+1};
    elseif strcmpi(varargin{vn},'f0'), f0=varargin{vn+1};
    elseif strcmpi(varargin{vn},'fmin'), fmin=varargin{vn+1};
    elseif strcmpi(varargin{vn},'fmax'), fmax=varargin{vn+1};
    elseif strcmpi(varargin{vn},'nv'), nv=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Padding'), PadMode=varargin{vn+1};
    elseif strcmpi(varargin{vn},'RelTol'), RelTol=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Preprocess'), Preprocess=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Plot'), PlotMode=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Display'), DispMode=varargin{vn+1};
    elseif strcmpi(varargin{vn},'CutEdges'), CutEdges=varargin{vn+1};
    else
        error(['There is no Property "',varargin{vn},'" (which is ',num2str(1+(vn-vst)/2,'%d'),...
            '''th out of ',num2str(ceil((nargin-1-vst)/2),'%d'),' specified)']);
    end
end
    


%========================= Wavelet function ===============================
%[fwt] and [twf] - wavelet function in frequency and time
%wp - structure with wavelet parameters containing fields:
%     xi1,xi2 - wavelet full support in the frequency domain;
%     ompeak,tpeak - wavelet peak frequency and peak time
%     t1,t2 - wavelet full support in the time domain;
%     fwtmax,twfmax - maximum value of abs(fwt) and abs(twf);
%     C,D - coefficients needed for reconstruction
%--------------------------------------------------------------------------
fwt=[]; twf=[];
if recflag==1
    wp=struct;
    wp.fwtmax=[]; wp.twfmax=[]; wp.C=[]; wp.D=[];
    wp.xi1=-Inf; wp.xi2=Inf; wp.ompeak=[];
    wp.t1=-Inf; wp.t2=Inf; wp.tpeak=[];
end
if iscell(Wavelet)
    fwt=Wavelet{1};
    if length(Wavelet{2})==2, wp.xi1=Wavelet{2}(1); wp.xi2=Wavelet{2}(2); end
    twf=Wavelet{3};
    if length(Wavelet{4})==2, wp.t1=Wavelet{4}(1); wp.t2=Wavelet{4}(2); end
elseif strcmpi(Wavelet,'Lognorm')
    q=2*pi*f0;
    fwt=@(xi)exp(-(q^2/2)*(log(xi).^2));
    wp.xi1=0; wp.ompeak=1;
    wp.C=sqrt(pi/2)/q; wp.D=wp.C*exp(1/(2*q^2));
elseif strcmpi(Wavelet,'Morlet')
    om0=2*pi*f0; %just for convenience denote circular central frequency
    fwt=@(xi)(exp(-(1/2)*(om0-xi).^2)-exp(-(1/2)*(om0^2+xi.^2)));
    if f0>=1, twf=@(t)(1/sqrt(2*pi))*(exp(1i*om0*t)-exp(-om0^2/2)).*exp(-t.^2/2); end
    wp.D=Inf;
elseif strcmpi(Wavelet,'Bump')
    q=2.5*f0;
    if q<1, error('For Bump wavelet f0 cannot be lower than 0.4'); end
    fwt=@(xi)exp(1-abs(1./(1-(q^2)*(1-xi).^2)));
    wp.xi1=max([0,1-1/q]); wp.xi2=1+1/q;
    wp.ompeak=1;
elseif ~isempty(strfind(lower(Wavelet),'morse'))
    a=3; if length(Wavelet)>5, a=str2double(Wavelet(7:length(Wavelet))); end
    q=30*f0^2/a;
    B=(exp(1)*a/q).^(q/a);
    fwt=@(xi)exp(-xi.^a+q*(log(xi)+(1/a)*log(exp(1)*a/q))); wp.xi1=0;
    wp.ompeak=(q/a)^(1/a);
    wp.fwtmax=fwt(wp.ompeak); wp.C=(1/2)*(B/a)*gamma(q/a);
    wp.D=Inf; if q>1, wp.D=(wp.ompeak/2)*(B/a)*gamma((q-1)/a); end
else
    error('Invalid wavelet name');
end
%==========================================================================

%Determine parameters of the wavelet function (and nv if needed)
if recflag==1
    if strcmpi(DispMode,'on'), fprintf('Estimating wavelet parameters...\n'); end
    parcalc(RelTol); %calculate/update wavelet parameters
end
if isempty(fmin), fmin=(wp.ompeak/2/pi)*(wp.t2e-wp.t1e)*fs/L; end %if not specified, determine minimum possible WT frequency for accuracy RelTol
if fmin>fmax, error('Minimal frequency %1.2e, either specified by you or determined automatically from the cone-of-influence, exceeds maximal frequency %1.2e!',fmin,fmax); end


%Produce warning if there are no WT in the cone-of-influence
if (wp.t2e-wp.t1e)*wp.ompeak/(2*pi*fmax)>L/fs
    if ~strcmpi(DispMode,'off')
        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
        wstr='';
        wstr=[wstr,sprintf('For used wavelet function parameters and signal time-length there are no WT coefficients determined with specified accuracy %1.2e.\n',RelTol)];
        wstr=[wstr,sprintf('Transform is inaccurate under the specified precision. Either allow lower accuracy by increasing ''RelTol'', or consider changing wavelet resolution parameter(s).\n')];
        if strcmpi(CutEdges,'on')
            wstr=[wstr,sprintf('Continuing without cone-of-influence.\n')];
        end
        fprintf(wstr);
        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
    end
    if strcmpi(CutEdges,'on'), CutEdges='off'; end
elseif (wp.t2e-wp.t1e)*wp.ompeak/(2*pi*fmin)/(L/fs)>1+2*eps
    if ~strcmpi(DispMode,'off')
        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
        wstr='';
        wstr=[wstr,sprintf('At lowest frequency %1.2e there are no WT coefficients determined with specified accuracy %1.2e.\n',fmin,RelTol)];
        wstr=[wstr,sprintf('There is nothing crucial in this, but WT computation might be slower due to using more padding.\n')];
        wstr=[wstr,sprintf('If you want to use only WT from the cone of influence (i.e. determined with specified precision),\n')];
        wstr=[wstr,sprintf('then better increase the minimal frequency to at least %1.2e or change wavelet resolution parameter(s).\n',(wp.t2e-wp.t1e)*(wp.ompeak/2/pi)/(L/fs))];
        fprintf(wstr);
        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
    end
end


%Define frequencies
nvsim=nv; wp.nv=nv;
if ~isempty(strfind(nv,'auto')) %determine number-of-voices [nv] if needed
    Nb=10; if length(nv)>4, Nb=str2double(nv(6:length(nv))); end
    wp.nv=Nb*log(2)/log(wp.xi2h/wp.xi1h);
    nv=ceil(wp.nv);
    if strcmpi(DispMode,'on')
        fprintf('Optimal nv ("number-of-voices" needed for frequency binning) was determined to be %0.2f (rounded to %d)\n',wp.nv,ceil(wp.nv));
    end
end
freq=2.^((ceil(nv*log2(fmin)):floor(nv*log2(fmax)))'/nv); % frequencies (exact values depend only on [nv])
SN=length(freq); %number of frequencies
coib1=ceil(abs(wp.t1e*fs*(wp.ompeak./(2*pi*freq)))); coib2=ceil(abs(wp.t2e*fs*(wp.ompeak./(2*pi*freq)))); %cone of influence edges


%======== Signal preprocessing: detrending, filtering and padding =========
%[dflag] determines to do detrending and filtering before or after padding
dflag=0;
if ~iscell(PadMode)
    if ~ischar(PadMode) && ~isempty(PadMode(PadMode~=0)), dflag=1; end
    if strcmpi(PadMode,'predictive') && fmin<5*fs/L, dflag=1; end
else
    if ~ischar(PadMode{1}) && ~isempty(PadMode{1}(PadMode{1}~=0)), dflag=1; end
    if ~ischar(PadMode{2}) && ~isempty(PadMode{2}(PadMode{2}~=0)), dflag=1; end
    if strcmpi(PadMode{1},'predictive') && fmin<5*fs/L, dflag=1; end
    if strcmpi(PadMode{2},'predictive') && fmin<5*fs/L, dflag=1; end
end
%Detrend (subtract third-order polynomial fit) and filter first for usual padding
if strcmpi(Preprocess,'on') && dflag==0
    %Detrending
    X=(1:L)'/fs; XM=[ones(L,1),X,X.^2,X.^3];
    w=warning('off','all'); p3fit=XM*(XM\signal); warning(w);
    signal=signal-p3fit;
    %Filtering
    fx=fft(signal,L); % Fourier transform of a signal
    Nq=ceil((L+1)/2); ff=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L; ff=ff(:); % frequencies in Fourier transform
    fx(abs(ff)<=max([fmin,fs/L]) | abs(ff)>=fmax)=0; % filter signal in a chosen frequency domain
    signal=ifft(fx);
end
%Padding
NL=2^nextpow2(L+coib1(1)+coib2(1));
if coib1(1)==0 && coib2(1)==0, n1=floor((NL-L)/2); n2=ceil((NL-L)/2);
else
    n1=floor((NL-L)*coib1(1)/(coib1(1)+coib2(1)));
    n2=ceil((NL-L)*coib2(1)/(coib1(1)+coib2(1)));
end
if strcmpi(DispMode,'on')
    if strcmpi(Preprocess,'on')
        fprintf('Signal preprocessing (detrending, then filtering) and padding (%d values to the left and %d to the right)...\n',n1,n2);
    else
        fprintf('Padding (%d values to the left and %d to the right)...\n',n1,n2);
    end
end
if iscell(PadMode)
    p1=PadMode{1}; p2=PadMode{2};
    %Padding to the left
    if strcmpi(p1,'predictive') || strcmpi(p2,'predictive')
        if strcmpi(DispMode,'on'), fprintf('Applying predictive padding (may take some time)...\n'); end
        dflag=1;
    end
    if ~ischar(p1) && length(p1)>1
        p1=p1(:); PL1=length(p1);
        padleft=[zeros(n1-PL1,1);p1(PL1-min([n1,PL1])+1:PL1)];
    elseif strcmpi(p1,'predictive')
        w=2.^(-(2*pi*fmin/wp.ompeak)*(L/fs-(1:L)/fs)/(wp.t2h-wp.t1h));
        padleft=fcast(flipud(signal),fs,n1,[max([fmin,fs/L]),fmax],min([ceil(SN/2)+5,round(L/3)]),w); padleft=flipud(padleft);
    elseif isnumeric(p1), padleft=p1*ones(n1,1);
    elseif strcmpi(p1,'symmetric'), padleft=[zeros(n1-L,1);flipud(signal(1:min([n1,L])))];
    elseif strcmpi(p1,'periodic'), padleft=[zeros(n1-L,1);signal(L-min([n1,L])+1:L)];
    else error('Bad ''Padding'' property');
    end
    %Padding to the right
    if ~ischar(p2) && length(p2)>1
        p2=p2(:); PL2=length(p2);
        padright=[p2(1:min([n2,PL2]));zeros(n2-PL2,1)];
    elseif strcmpi(p2,'predictive')
        w=2.^(-(2*pi*fmin/wp.ompeak)*(L/fs-(1:L)/fs)/(wp.t2h-wp.t1h));
        padright=fcast(signal,fs,n2,[max([fmin,fs/L]),fmax],min([ceil(SN/2)+5,round(L/3)]),w);
    elseif isnumeric(p2), padright=p2*ones(n2,1);
    elseif strcmpi(p2,'symmetric'), padright=[flipud(signal(L-min([n2,L])+1:L));zeros(n2-L,1)];
    elseif strcmpi(p2,'periodic'), padright=[signal(1:min([n2,L]));zeros(n2-L,1)];
    else error('Bad ''Padding'' property');
    end
    signal=[padleft;signal;padright];
elseif strcmpi(PadMode,'predictive')
    if strcmpi(DispMode,'on'), fprintf('Applying predictive padding (may take some time)...\n'); end
    w=2.^(-(2*pi*fmin/wp.ompeak)*(L/fs-(1:L)/fs)/(wp.t2h-wp.t1h));
    padleft=fcast(flipud(signal),fs,n1,[max([fmin,fs/L]),fmax],min([ceil(SN/2)+5,round(L/3)]),w); padleft=flipud(padleft);
    padright=fcast(signal,fs,n2,[max([fmin,fs/L]),fmax],min([ceil(SN/2)+5,round(L/3)]),w);
    signal=[padleft;signal;padright];
    dflag=1; %to detrend one more time
elseif isnumeric(PadMode)
    signal=[PadMode*ones(n1,1);signal;PadMode*ones(n2,1)]; %padding with predefined values (by default with zeros)
elseif strcmpi(PadMode,'symmetric'),
    signal=[zeros(n1-L,1);flipud(signal(1:min([n1,L])));signal;flipud(signal(L-min([n2,L])+1:L));zeros(n2-L,1)]; %symmetric padding
elseif strcmpi(PadMode,'periodic'),
    signal=[zeros(n1-L,1);signal(L-min([n1,L])+1:L);signal;signal(1:min([n2,L]));zeros(n2-L,1)]; %periodic padding
elseif strcmpi(PadMode,'none')
    NL=L; n1=0; n2=0; %no padding (quite equivalent to periodic padding)
else
    error('Bad ''Padding'' property');
end
%Detrend (subtract third-order polynomial fit) after padding for special cases
if strcmpi(Preprocess,'on') && dflag==1
    X=(0:length(signal)-1)'/fs; XM=[ones(length(signal),1),X,X.^2,X.^3];
    w=warning('off','all'); p3fit=XM*(XM\signal); warning(w);
    signal=signal-p3fit;
end
%Filtering of the padded signal
Nq=ceil((NL+1)/2); ff=[(0:Nq-1),-fliplr(1:NL-Nq)]*fs/NL; ff=ff(:); % frequencies in Fourier transform
fx=fft(signal,NL); fx(ff<=0)=0; % Fourier transform of a signal (set to zero at negative frequencies)
if strcmpi(Preprocess,'on')
    fx(ff<=max([fmin,fs/L]) | ff>=fmax)=0; % filter signal in a chosen frequency domain
end
%--------------------------------------------------------------------------


%Wavelet transform by itself
WT=zeros(SN,L)*NaN; ouflag=0; if (wp.t2e-wp.t1e)*wp.ompeak/(2*pi*fmax)>L/fs, coib1(:)=0; coib2(:)=0; end
if strcmpi(DispMode,'on'), pos=0; fprintf('Calculating Wavelet Transform (%d frequencies from %0.3f to %0.3f): ',SN,freq(1),freq(end)); end
for sn=1:SN
    freqwf=ff*wp.ompeak/(2*pi*freq(sn)); %frequencies for the wavelet function
    ii=find(freqwf>wp.xi1/2/pi & freqwf<wp.xi2/2/pi); %take into account only frequencies within the wavelet support
    if ~isempty(fwt)
        fw=conj(fwt(2*pi*freqwf(ii))); nid=find(isnan(fw) | ~isfinite(fw));
        if ~isempty(nid) %to avoid NaNs due to numerics, e.g. sin(0)/0
            fw(nid)=conj(fwt(2*pi*freqwf(ii(nid))+10^(-14)));
            nid=find(isnan(fw) | ~isfinite(fw)); fw(nid)=0;
            if ~isempty(nid), ouflag=1; ouval=2*pi*freqwf(nid(1)); end
        end
    else
        timewf=(2*pi*freq(sn)/wp.ompeak)*(1/fs)*[-(1:ceil((NL-1)/2))+1,NL+1-(ceil((NL-1)/2)+1:NL)]';
        jj=find(timewf>wp.t1 & timewf<wp.t2); tw=zeros(NL,1); %take into account only times within the wavelet support
        tw(jj)=conj(twf(timewf(jj))); nid=find(isnan(tw) | ~isfinite(tw));
        if ~isempty(nid) %to avoid NaNs due to numerics, e.g. sin(0)/0
            tw(nid)=conj(twf(timewf(nid)+10^(-14)));
            nid=find(isnan(tw) | ~isfinite(tw)); tw(nid)=0;
            if ~isempty(nid), ouflag=1; ouval=timewf(nid(1)); end
        end
        fw=(1/fs)*fft(tw); fw=fw(ii);
    end
    cc=zeros(NL,1); cc(ii)=fx(ii).*fw(:); %convolution in the frequency domain
    out=((wp.ompeak/(2*pi*freq(sn)))^(1-p))*ifft(cc,NL); % calculate WT at each time
    WT(sn,1:L)=out(1+n1:NL-n2);
    
    if strcmpi(DispMode,'on') && floor(100*sn/SN)>floor(100*(sn-1)/SN)
        cstr=num2str(floor(100*sn/SN)); fprintf([repmat('\b',1,pos),cstr,'%%']); pos=length(cstr)+1;
    end
end
if strcmpi(DispMode,'on'), fprintf('\n'); end
if ouflag==1
    if ~isempty(fwt)
        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
        fprintf('Possibly overflow/underflow (e.g. Inf/Inf): specified frequency-domain form of wavelet function\n');
        fprintf('returns NaN or Inf, e.g. when its argument is %e. In all such cases it is set to zero.\n',ouval);
        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
    else
        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
        fprintf('Possibly overflow/underflow (e.g. Inf/Inf): specified time-domain form of wavelet function\n');
        fprintf('returns NaN or Inf, e.g. when its argument is %e. In all such cases it is set to zero.\n',ouval);
        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
    end
end

%Set to NaN all WT coefficients outside the cone of influence if specified
if strcmpi(CutEdges,'on')
    icoib=find(L-coib1-coib2<=0); WT(icoib,:)=NaN;
    ovL=ceil(sum(coib1+coib2)-L*length(icoib)); frn=zeros(ovL,1)*NaN; ttn=zeros(ovL,1)*NaN; qn=0;
    for fn=1:SN, cL=coib1(fn)+coib2(fn);
        if cL>0 && cL<L, frn(1+qn:qn+cL)=fn; ttn(1+qn:qn+cL)=[1:coib1(fn),(L-coib2(fn)+1):L]'; qn=qn+cL; end
    end
    frn=frn(1:qn); ttn=ttn(1:qn);
    lid=sub2ind([SN,L],frn,ttn); WT(lid)=NaN;
end

%Plotting WT if specified
if ~strcmpi(PlotMode,'off')
    if strcmpi(DispMode,'on'), fprintf('Plotting...\n'); end
    scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,scrsz(3)/2,6*scrsz(4)/8]);
    axes('Position',[0.15,0.1,0.8,0.5333],'Layer','top','XScale','log','Box','on','FontSize',16);
    hold all;
    
    YY=freq; XX=(0:(L-1))/fs; ZZ=abs(WT); ZZname='WT amplitude';
    if ~isempty(strfind(lower(PlotMode),'pow')), ZZ=ZZ.^2; ZZname='WT power'; end
    MYL=round(scrsz(3)); MXL=round(scrsz(4)); %maximum number of points seen in plots
    if isempty(strfind(lower(PlotMode),'wr')) && (size(ZZ,1)>MYL || size(ZZ,2)>MXL)
        if ~strcmpi(DispMode,'off')
            fprintf('Warning: WT contains more data points (%d x %d) than pixels in the plot, so for a\n',size(ZZ,1),size(ZZ,2));
            fprintf('         better performance its resampled version (%d x %d) will be displayed instead.\n',min([MYL,size(ZZ,1)]),min([MXL,size(ZZ,2)]));
        end
        if size(ZZ,1)>MYL, YY=exp(linspace(log(freq(1)),log(freq(end)),MYL)); end
        if size(ZZ,2)>MXL, XX=linspace(0,(L-1)/fs,MXL); end
        ZZ=aminterp((0:(L-1))/fs,log(freq),ZZ,XX,log(YY),'max'); XX=XX(:); YY=YY(:);
    end
    TL=length(XX); FL=length(YY);
    pc=pcolor(YY,XX,ZZ'); set(pc,'EdgeColor','none'); title(ZZname);
    
    xlabel('Frequency (Hz)'); ylabel('Time (s)');
    ylim([0,(L-1)/fs]); xlim([freq(1),freq(end)]);
    
    coib1(coib1==0)=NaN; coib2(coib2==0)=NaN;
    ib=find(L-coib1-coib2>0,1,'first'); if isempty(ib), ib=1; end
    plot(freq(end:-1:ib),coib1(end:-1:ib)/fs,'-k','LineWidth',2);
    plot(freq(ib:end),(L-coib2(ib:end)+1)/fs,'-k','LineWidth',2);
    plot([freq(ib),freq(ib)],[coib1(ib)/fs,(L-coib2(ib)+1)/fs],'-k','LineWidth',2);
    coib1(isnan(coib1))=0; coib2(isnan(coib2))=0; if strcmpi(CutEdges,'off'), coib1(:)=0; coib2(:)=0; end
    
    if ~isempty(strfind(lower(PlotMode),'+'))
        axes('Position',[0.15,0.7,0.8,0.25],'Layer','top','XLim',[freq(1),freq(end)],'XScale','log','XTickLabel',{},'Box','on','FontSize',16);
        hold all;
        mx=zeros(FL,1); for fn=1:FL, mx(fn)=mean(ZZ(fn,~isnan(ZZ(fn,:))),2); end
        plot(YY,mx,'-k','LineWidth',2); ylabel({'Time-averaged',ZZname});
        if max(mx)>0, ylim([0,1.1*max(mx)]); end
        
        if ~isempty(strfind(lower(PlotMode),'++'))
            sZZ=sort(ZZ,2); ZL=zeros(FL,1)*NaN;
            for fn=1:FL, uid=find(~isnan(sZZ(fn,:)),1,'last'); if ~isempty(uid), ZL(fn)=uid; end, end
            lx=zeros(FL,1)*NaN; ux=zeros(FL,1)*NaN;
            for fn=1:FL
                if ~isnan(ZL(fn))
                    lx(fn)=sZZ(fn,max([1,round(0.025*ZL(fn))]));
                    ux(fn)=sZZ(fn,round(0.975*ZL(fn)));
                end
            end
            idnn=find(~isnan(lx)); % not-NaN indices
            hold all; fill([YY(idnn);flipud(YY(idnn))],[lx(idnn);flipud(ux(idnn))],[1,1,0],'FaceAlpha',0.5);
            if max(ux)>0, ylim([0,1.1*max(ux)]); end
        end
    end
end

if nargout>2
    wopt=struct; %simulation parameters
    wopt.wp=wp; %parameters of the wavelet
    wopt.TFRname='WT'; wopt.fs=fs;
    wopt.Wavelet=Wavelet;
    wopt.f0=f0;
    wopt.fmin=fmin;
    wopt.fmax=fmax;
    wopt.nv=nv; wopt.nvsim=nvsim;
    wopt.Padding=PadMode;
    wopt.RelTol=RelTol;
    wopt.Preprocess=Preprocess;
    wopt.Plot=PlotMode;
    wopt.Display=DispMode;
    wopt.CutEdges=CutEdges;
    
    varargout{1}=wopt;
end










%----------------------------------------------------------------------------------------------------------------------
%======================================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================================================================
%----------------------------------------------------------------------------------------------------------------------
    % Based on wavelet function in time [twf] and frequency [fwt],
    % determines wavelet parameters such as constant [Cpsi],
    % peaks in time [tpeak] and frequency [ompeak],
    % epsilon-support in frequency [xi1e,xi2e] and in time [t1e,t2e],
    % 50%-support in frequency [xi1h,xi2h] and time [t1h,t2h],
    % time-frequency resolution [tfres] (inverse of multiplication of the latter),
    % and number of voices [nv] (if 'auto') for specified relative accuracy [racc] (=epsilon);
    % assigns all these values into the wavelet parameters structure [wp].
    function parcalc(racc)
        racc=min(racc,1-10^(-6)); %current \epsilon
        ctol=max([racc/1000,10^(-12)]); %parameter of numerical accuracy
        MIC=max([100000,10*L]); %maximum interval count for one-time calculations

        %==================================================================
        %Determine values for known frequency and/or time-domain forms
        if ~isempty(fwt) %if the frequency-domain form is known
            wp.fwt=fwt;
            if isempty(wp.ompeak) %peak frequency
                wp.ompeak=1; if strcmpi(Wavelet,'Morlet'), wp.ompeak=2*pi*f0; end
                if wp.xi1>0 && isfinite(wp.xi2), wp.ompeak=sqrt(wp.xi1*wp.xi2);
                elseif isfinite(wp.xi2), wp.ompeak=wp.xi2/2; end
                if fwt(wp.ompeak)==0 || isnan(fwt(wp.ompeak)) || ~isfinite(fwt(wp.ompeak))
                    cp1=wp.ompeak*exp(-10^(-14)); cp2=wp.ompeak*exp(10^(-14)); kk=1;
                    while kk<10^(28)
                        cv1=abs(fwt(cp1)); cv2=abs(fwt(cp2)); kk=kk*2;
                        if isfinite(cv1) && cv1>0, wp.ompeak=cp1; break; end
                        if isfinite(cv2) && cv2>0, wp.ompeak=cp2; break; end
                        cp1=cp1*exp(-kk*10^(-14)); if cp1<=max([wp.xi1,0]), cp1=(cp1*exp(kk*10^(-14))+max([wp.xi1,0]))/2; end
                        cp2=cp2*exp(kk*10^(-14)); if cp2>=wp.xi2, cp2=(cp2*exp(-kk*10^(-14))+wp.xi2)/2; end
                    end
                    cv=abs(fwt(wp.ompeak));
                    while isnan(cv) || cv==0 %if search failed
                        if isfinite(wp.xi2), pp=max([wp.xi1,0])+(wp.xi2-max([wp.xi1,0]))*rand(MIC,1);
                        else pp=exp(atan(pi*(rand(MIC,1)-1/2))); end
                        [cv,ipeak]=max(abs(fwt(pp))); wp.ompeak=pp(ipeak);
                    end
                end
                wp.ompeak=fminsearch(@(x)-abs(fwt(exp(x))),log(wp.ompeak),optimset('TolX',10^(-14),'Display','off'));
                wp.ompeak=exp(wp.ompeak);
            end
            if isempty(wp.fwtmax)
                wp.fwtmax=fwt(wp.ompeak);
                if isnan(wp.fwtmax), wp.fwtmax=fwt(wp.ompeak+10^(-14)); end
            end
            vfun=@(u)conj(fwt(exp(u))); xp=log(wp.ompeak); lim1=log(max([wp.xi1,0])); lim2=log(wp.xi2);
            
            %Test admissibility %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            wstate=warning('off','all');
            if wp.xi1<=0, AC=fwt(0); else AC=0; end
            if isnan(AC)
                cx0=10^(-14);
                while fwt(cx0)>10^(-14), cx0=cx0/2; end
                while isnan(fwt(cx0)), cx0=cx0*2; end
                AC=fwt(cx0);
            end
            if AC>10^(-12) && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                fprintf(2,'Wavelet does not seem to be admissible (its Fourier transform does not vanish at zero frequency)!\n');
                fprintf(2,'Parameters estimated from its frequency domain form, e.g. integration constant Cpsi (which is \n');
                fprintf(2,'infinite for non-admissible wavelets), cannot be estimated appropriately (the same concerns the \n');
                fprintf(2,'number-of-voices ''nv'', when set to ''auto'', so frequency discretization might be also not appropriate).\n');
                fprintf(2,'It is recommended to use only admissible wavelets.\n');
                fprintf(2,'----------------------------------------------------------------------------------------------------\n');
            end
            warning(wstate);
            
            [QQ,wflag,xx,ss]=sqeps(vfun,xp,[lim1,lim2],racc,MIC,...
                [log((wp.ompeak/fmax)*fs/L/8),log(8*(wp.ompeak/(fs/L))*fs)]); %¬¬¬¬¬¬¬¬
            wp.xi1e=exp(ss(1,1)); wp.xi2e=exp(ss(1,2)); wp.xi1h=exp(ss(2,1)); wp.xi2h=exp(ss(2,2));
            if isempty(wp.C), wp.C=(QQ(1,1)+QQ(1,2))/2; end %¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            if isempty(wp.D) %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                wstate=warning('off','all');
                [D1,errD1]=quadgk(@(u)conj(fwt(1./u)),1/wp.ompeak,exp(-xx(1,1)),'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                [D2,errD2]=quadgk(@(u)-conj(fwt(1./u)),1/wp.ompeak,exp(-xx(1,2)),'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                [D3,errD3]=quadgk(@(u)conj(fwt(1./u)),exp(-xx(1,1)),exp(-xx(4,1)),'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                [D4,errD4]=quadgk(@(u)-conj(fwt(1./u)),exp(-xx(1,2)),exp(-xx(4,2)),'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                if abs((errD1+errD2+errD3+errD4)/(D1+D2+D3+D4))<10^(-4), wp.D=(wp.ompeak/2)*(D1+D2+D3+D4); else wp.D=Inf; end
                warning(wstate);
            end %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            if wflag==1 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                fprintf('The frequency-domain wavelet function is not well-behaved (e.g. decays very slowly as frequency tends to zero\n');
                fprintf('or infinity). The integration might be not accurate (and therefore e.g. the calculated number-of-voices ''nv'',\n');
                fprintf('if set to ''auto'', so frequency discretization might be also not appropriate).\n');
                fprintf(2,'----------------------------------------------------------------------------------------------------\n');
            end
            
            if isempty(twf) %if time domain form is not known
                [PP,wflag,xx,ss]=sqeps(@(x)abs(fwt(x)).^2,wp.ompeak,[max([wp.xi1,0]),wp.xi2],racc,MIC,...
                    [0,8*(wp.ompeak/(fs/L))*fs]); %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                Etot=sum(PP(1,:))/2/pi;
                
                CL=2^nextpow2(MIC/8); CT=CL/(2*abs(ss(1,2)-ss(1,1)));
                CNq=ceil((CL+1)/2); cxi=(2*pi/CT)*(CNq-CL:CNq-1)'; idm=find(cxi<=max([wp.xi1,0])); idc=find(cxi>max([wp.xi1,0]) & cxi<wp.xi2); idp=find(cxi>=wp.xi2);
                Cfwt=[zeros(length(idm),1);fwt(cxi(idc));zeros(length(idp),1)]; idnan=find(isnan(Cfwt));
                if ~isempty(idnan), idnorm=find(~isnan(Cfwt)); Cfwt(idnan)=interp1(idnorm,Cfwt(idnorm),idnan,'spline','extrap'); end
                Ctwf=ifft((CL/CT)*Cfwt([CL-CNq+1:CL,1:CL-CNq])); Ctwf=Ctwf([CNq+1:CL,1:CNq]);
                Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                Eest=(CT/CL)*sum(Etwf); perr=Inf;
                while (abs(Etot-Eest)+Iest1+Iest2)/Etot<=perr
                    CT=CT/2; perr=(abs(Etot-Eest)+Iest1+Iest2)/Etot;
                    CNq=ceil((CL+1)/2); cxi=(2*pi/CT)*(CNq-CL:CNq-1)'; idm=find(cxi<=max([wp.xi1,0])); idc=find(cxi>max([wp.xi1,0]) & cxi<wp.xi2); idp=find(cxi>=wp.xi2);
                    Cfwt=[zeros(length(idm),1);fwt(cxi(idc));zeros(length(idp),1)]; idnan=find(isnan(Cfwt));
                    if ~isempty(idnan), idnorm=find(~isnan(Cfwt)); Cfwt(idnan)=interp1(idnorm,Cfwt(idnorm),idnan,'spline','extrap'); end
                    Ctwf=ifft((CL/CT)*Cfwt([CL-CNq+1:CL,1:CL-CNq])); Ctwf=Ctwf([CNq+1:CL,1:CNq]);
                    Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                    Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                    Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                    Eest=(CT/CL)*sum(Etwf);
                end
                CL=16*CL; CT=CT*2;
                CNq=ceil((CL+1)/2); cxi=(2*pi/CT)*(CNq-CL:CNq-1)'; idm=find(cxi<=max([wp.xi1,0])); idc=find(cxi>max([wp.xi1,0]) & cxi<wp.xi2); idp=find(cxi>=wp.xi2);
                Cfwt=[zeros(length(idm),1);fwt(cxi(idc));zeros(length(idp),1)]; idnan=find(isnan(Cfwt));
                if ~isempty(idnan), idnorm=find(~isnan(Cfwt)); Cfwt(idnan)=interp1(idnorm,Cfwt(idnorm),idnan,'spline','extrap'); end
                Ctwf=ifft((CL/CT)*Cfwt([CL-CNq+1:CL,1:CL-CNq])); Ctwf=Ctwf([CNq+1:CL,1:CNq]);
                Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                Eest=(CT/CL)*sum(Etwf);
                
                if (abs(Etot-Eest)+Iest1+Iest2)/Etot>0.01 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                    fprintf(['Cannot accurately invert the specified frequency-domain form of the wavelet function to find its\n',...
                        'time domain form and corresponding characteristics (e.g. cone-of-influence borders).\n',...
                        'This might be because the wavelet function decays too slowly in time or frequency.\n']);
                    fprintf(2,'----------------------------------------------------------------------------------------------------\n');
                end
                
                Ctwf=Ctwf(1:2*CNq-3); ct=(CT/CL)*(-(CNq-2):CNq-2)'; %make symmetric
                wp.twf={Ctwf,ct};
                Ctwf=Ctwf.*exp(-1i*wp.ompeak*ct); %demodulate (wavelet only) %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                
                %Estimate general parameters
                if isempty(wp.tpeak) %peak time
                    ipeak=find(abs(Ctwf)==max(abs(Ctwf)));
                    if length(ipeak)==1
                        a1=abs(Ctwf(ipeak-1)); a2=abs(Ctwf(ipeak)); a3=abs(Ctwf(ipeak+1));
                        wp.tpeak=ct(ipeak);
                        if abs(a1-2*a2+a3)>2*eps, %use queadratic interpolation to find exact peak location
                            wp.tpeak=wp.tpeak+(1/2)*(a1-a3)/(a1-2*a2+a3)*(CT/CL);
                        end
                    else
                        wp.tpeak=mean(ct(ipeak));
                    end
                end
                if isempty(wp.twfmax)
                    [~,ipeak]=min(abs(ct-wp.tpeak));
                    wp.twfmax=interp1(ct(ipeak-1:ipeak+1),abs(Ctwf(ipeak-1:ipeak+1)),wp.tpeak,'spline');
                end
                
                %Calculate the cumulative integrals
                ct=[ct-CT/CL/2;ct(end)+CT/CL/2]; %to use midpoint rule
                CS=(CT/CL)*cumsum(Ctwf); CS=[0;CS(:)]/CS(end); CS=abs(CS);
                ICS=(CT/CL)*cumsum(Ctwf(end:-1:1)); ICS=ICS(end:-1:1); ICS=[ICS(:);0]/ICS(1); ICS=abs(ICS);
                
                %Estimate epsilon-supports
                xid=find(CS(1:end-1)<racc/2 & CS(2:end)>=racc/2,1,'first');
                if isempty(xid), wp.t1e=ct(1);
                else
                    a1=CS(xid)-racc/2; a2=CS(xid+1)-racc/2;
                    wp.t1e=ct(xid)-a1*(ct(xid+1)-ct(xid))/(a2-a1);
                end
                xid=find(ICS(1:end-1)>=racc/2 & ICS(2:end)<racc/2,1,'last');
                if isempty(xid), wp.t2e=ct(end);
                else
                    a1=ICS(xid)-racc/2; a2=ICS(xid+1)-racc/2;
                    wp.t2e=ct(xid)-a1*(ct(xid+1)-ct(xid))/(a2-a1);
                end
                xid=find(CS(1:end-1)<0.25 & CS(2:end)>=0.25,1,'first');
                if isempty(xid), wp.t1h=ct(1);
                else
                    a1=CS(xid)-0.25; a2=CS(xid+1)-0.25;
                    wp.t1h=ct(xid)-a1*(ct(xid+1)-ct(xid))/(a2-a1);
                end
                xid=find(ICS(1:end-1)>=0.25 & ICS(2:end)<0.25,1,'last');
                if isempty(xid), wp.t2h=ct(end);
                else
                    a1=ICS(xid)-0.25; a2=ICS(xid+1)-0.25;
                    wp.t2h=ct(xid)-a1*(ct(xid+1)-ct(xid))/(a2-a1);
                end
                
            end
            
        end
        %------------------------------------------------------------------
        if ~isempty(twf) %if the time-domain form is known
            wp.twf=twf;
            if isempty(wp.tpeak) %peak time
                wp.tpeak=max([min([0,wp.t2-abs(wp.t2)/2]),wp.t1+abs(wp.t1)/2]);
                if isfinite(wp.t1) && isfinite(wp.t2), wp.tpeak=(wp.t1+wp.t2)/2; end
                if twf(wp.tpeak)==0 || isnan(twf(wp.tpeak)) || ~isfinite(twf(wp.tpeak))
                    cp1=wp.tpeak-10^(-14); cp2=wp.tpeak+10^(-14); kk=1;
                    while kk<10^(28)
                        cv1=abs(twf(cp1)); cv2=abs(twf(cp2)); kk=kk*2;
                        if isfinite(cv1) && cv1>0, wp.tpeak=cp1; break; end
                        if isfinite(cv2) && cv2>0, wp.tpeak=cp2; break; end
                        cp1=cp1-kk*10^(-14); if cp1<=wp.t1, cp1=(cp1+kk*10^(-14)+wp.t1)/2; end
                        cp2=cp2+kk*10^(-14); if cp2>=wp.t2, cp2=(cp2-kk*10^(-14)+wp.t2)/2; end
                    end
                    cv=abs(twf(wp.tpeak));
                    while isnan(cv) || cv==0 %if search failed
                        if isfinite(wp.t1) && isfinite(wp.t2), pp=wp.t1+(wp.t2-wp.t1)*rand(MIC,1);
                        elseif isfinite(wp.t1), pp=wp.t1+atan((pi/2)*rand(MIC,1));
                        elseif isfinite(wp.t2), pp=wp.t2-atan((pi/2)*rand(MIC,1));
                        else pp=atan(pi*(rand(MIC,1)-1/2));
                        end
                        [cv,ipeak]=max(abs(twf(pp))); wp.tpeak=pp(ipeak);
                    end
                end
                wp.tpeak=fminsearch(@(x)-abs(twf(x)),wp.tpeak,optimset('TolX',10^(-14),'Display','off'));
            end
            if isempty(wp.twfmax)
                wp.twfmax=twf(wp.tpeak);
                if isnan(wp.twfmax), wp.twfmax=twf(wp.tpeak+10^(-14)); end
            end
            
            %Test admissibility %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            wstate=warning('off','all');
            AC=quadgk(@(u)-twf(u),wp.tpeak,xx(1,1),'MaxIntervalCount',MIC,'AbsTol',10^(-16),'RelTol',0)+...
                quadgk(@(u)-twf(u),xx(1,1),xx(4,1),'MaxIntervalCount',MIC,'AbsTol',10^(-16),'RelTol',0)+...
                quadgk(@(u)twf(u),wp.tpeak,xx(1,2),'MaxIntervalCount',MIC,'AbsTol',10^(-16),'RelTol',0)+...
                quadgk(@(u)twf(u),xx(1,2),xx(4,2),'MaxIntervalCount',MIC,'AbsTol',10^(-16),'RelTol',0);
            if AC>10^(-8) && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                fprintf(2,'Wavelet does not seem to be admissible (its Fourier transform does not vanish at zero frequency)!\n');
                fprintf(2,'Parameters estimated from its frequency domain form, e.g. integration constant Cpsi (which is \n');
                fprintf(2,'infinite for non-admissible wavelets), cannot be estimated appropriately (the same concerns the \n');
                fprintf(2,'number-of-voices ''nv'', when set to ''auto'', so frequency discretization might be also not appropriate).\n');
                fprintf(2,'It is recommended to use only admissible wavelets.\n');
                fprintf(2,'----------------------------------------------------------------------------------------------------\n');
            end
            warning(wstate);
            
            %Calculate the frequency domain characteristics first, if not
            %known (will be needed afterwards, mainly [wp.ompeak])
            if isempty(fwt) %if frequency domain form is not known
                compeak=wp.ompeak;
                if isempty(compeak) %if not known, roughly estimate the peak frequency
                    [~,~,~,bss]=sqeps(@(u)abs(twf(u)).^2,wp.tpeak,[wp.t1,wp.t2],0.01,MIC,[wp.t1,wp.t2]);
                    BL=2^(nextpow2(MIC)); BNq=ceil((BL+1)/2); BT=bss(1,2)-bss(1,1);
                    bt=linspace(bss(1,1),bss(1,2),BL)'; bxi=(2*pi/BT)*[0:BNq-1,BNq-BL:-1]';
                    Bfwt=fft(twf(bt)); ix=find(bxi>max([wp.xi1,0]) & bxi<wp.xi2);
                    [~,imax]=max(abs(Bfwt(ix))); compeak=bxi(ix(imax));
                end
                [PP,wflag,xx,ss]=sqeps(@(x)abs(twf(x)).^2,wp.tpeak,[wp.t1,wp.t2],racc,MIC,...
                    [-8*(2*pi*fmax/compeak)*L/fs,8*(2*pi*fmax/compeak)*L/fs]); %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                Etot=sum(PP(1,:));
                
                CL=2^nextpow2(MIC/8); CT=2*abs(ss(1,2)-ss(1,1));
                CNq=ceil((CL+1)/2); ct=(CT/CL)*(CNq-CL:CNq-1)'; idm=find(ct<=wp.t1); idc=find(ct>wp.t1 & ct<wp.t2); idp=find(ct>=wp.t2);
                Ctwf=[zeros(length(idm),1);twf(ct(idc));zeros(length(idp),1)]; idnan=find(isnan(Ctwf));
                if ~isempty(idnan), idnorm=find(~isnan(Ctwf)); Ctwf(idnan)=interp1(idnorm,Ctwf(idnorm),idnan,'spline','extrap'); end
                Cfwt=(CT/CL)*fft(Ctwf([CL-CNq+1:CL,1:CL-CNq])); Cfwt=Cfwt([CNq+1:CL,1:CNq]);
                Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                Eest=(1/CT)*sum(Efwt); perr=Inf;
                while (abs(Etot-Eest)+Iest1+Iest2)/Etot<=perr
                    CT=CT*2; perr=(abs(Etot-Eest)+Iest1+Iest2)/Etot;
                    CNq=ceil((CL+1)/2); ct=(CT/CL)*(CNq-CL:CNq-1)'; idm=find(ct<=wp.t1); idc=find(ct>wp.t1 & ct<wp.t2); idp=find(ct>=wp.t2);
                    Ctwf=[zeros(length(idm),1);twf(ct(idc));zeros(length(idp),1)]; idnan=find(isnan(Ctwf));
                    if ~isempty(idnan), idnorm=find(~isnan(Ctwf)); Ctwf(idnan)=interp1(idnorm,Ctwf(idnorm),idnan,'spline','extrap'); end
                    Cfwt=(CT/CL)*fft(Ctwf([CL-CNq+1:CL,1:CL-CNq])); Cfwt=Cfwt([CNq+1:CL,1:CNq]);
                    Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                    Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                    Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                    Eest=(1/CT)*sum(Efwt);
                end
                CL=16*CL; CT=CT*2;
                CNq=ceil((CL+1)/2); ct=(CT/CL)*(CNq-CL:CNq-1)'; idm=find(ct<=wp.t1); idc=find(ct>wp.t1 & ct<wp.t2); idp=find(ct>=wp.t2);
                Ctwf=[zeros(length(idm),1);twf(ct(idc));zeros(length(idp),1)]; idnan=find(isnan(Ctwf));
                if ~isempty(idnan), idnorm=find(~isnan(Ctwf)); Ctwf(idnan)=interp1(idnorm,Ctwf(idnorm),idnan,'spline','extrap'); end
                Cfwt=(CT/CL)*fft(Ctwf([CL-CNq+1:CL,1:CL-CNq]));
                Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt([CNq+1:CL,1:CNq])).^2;
                Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                Eest=(1/CT)*sum(Efwt);
                
                if (abs(Etot-Eest)+Iest1+Iest2)/Etot>0.01 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                    fprintf(['Cannot accurately invert the specified time-domain form of the wavelet function to find its\n',...
                        'frequency-domain form and corresponding characteristics (e.g. optimal number-of-voices ''nv'').\n',...
                        'This might be because the wavelet function decays too slowly in time or frequency.\n']);
                    fprintf(2,'----------------------------------------------------------------------------------------------------\n');
                end
                
                Cfwt=Cfwt(2:CNq); cxi=(2*pi/CT)*(1:CNq-1)'; %take only positive frequencies
                if 2*abs((CT/CL)*sum(Ctwf))<ctol %represent low-frequency structure
                    Atot=abs((2*pi/CT)*sum(Cfwt./cxi));
                    cxi0=cxi(1); Cfwt0=(CT/CL)*sum(Ctwf.*exp(-1i*cxi0*ct));
                    axi=NaN*zeros(CL,1); Afwt=NaN*zeros(CL,1); kn=1;
                    while 2*abs(Cfwt0)/Atot>ctol
                        cxi0=cxi0/2; Cfwt0=(CT/CL)*sum(Ctwf.*exp(-1i*cxi0*ct));
                        axi(kn)=cxi0; Afwt(kn)=Cfwt0; kn=kn+1;
                    end
                    axi=axi(1:kn-1); Afwt=Afwt(1:kn-1);
                else
                    ix=min([1+find(diff(abs(Cfwt(2:end)))<=0,1,'first'),length(Cfwt)]);
                    cxi0=interp1([0;abs(Cfwt(1:ix))],[0;cxi(1:ix)],ctol/2,'spline');
                    axi=[]; Afwt=[];
                end
                CS0=interp1([0;cxi],[0;Cfwt],cxi0/4,'spline'); CS0=2*CS0; %initial cumulative sum
                
                %Move to the logarithmic frequency scale
                [~,imxi]=max(abs(Cfwt)); %peak position
                bxi1=linspace(log(cxi0),log(cxi(1)),ceil(2*CL/3))'; zxi1=(bxi1(1:end-1)+bxi1(2:end))/2;
                bxi2=linspace(log(cxi(1)),log(cxi(imxi)),ceil(2*CL/3))'; zxi2=(bxi2(1:end-1)+bxi2(2:end))/2;
                bxi3=linspace(log(cxi(imxi)),log(cxi(end)),ceil(2*CL/3))'; zxi3=(bxi3(1:end-1)+bxi3(2:end))/2;
                zxi=[zxi1;zxi2;zxi3]; bxi=[bxi1(1:end-1);bxi2(1:end-1);bxi3]; dbxi=diff(bxi);
                Zfwt=interp1([0;axi;cxi],[0;Afwt;Cfwt],exp(zxi),'spline');
                wp.fwt={Zfwt,exp(zxi)};
                
                %Estimate general parameters
                if isempty(wp.ompeak) %peak frequency
                    ipeak=find(abs(Zfwt)==max(abs(Zfwt)));
                    if length(ipeak)==1
                        a1=abs(Zfwt(ipeak-1)); a2=abs(Zfwt(ipeak)); a3=abs(Zfwt(ipeak+1));
                        wp.ompeak=zxi(ipeak);
                        if abs(a1-2*a2+a3)>2*eps, %use queadratic interpolation to find exact peak location
                            wp.ompeak=wp.ompeak+(1/2)*(a1-a3)/(a1-2*a2+a3)*dbxi(ipeak);
                        end
                    else
                        wp.ompeak=mean(zxi(ipeak));
                    end
                    wp.ompeak=exp(wp.ompeak);
                end
                if isempty(wp.fwtmax)
                    [~,ipeak]=min(abs(cxi-wp.ompeak));
                    wp.fwtmax=interp1(cxi(ipeak-1:ipeak+1),abs(Cfwt(ipeak-1:ipeak+1)),wp.ompeak,'spline');
                end
                if isempty(wp.C), wp.C=(1/2)*sum(conj(Zfwt).*dbxi); end %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                if isempty(wp.D) %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    wp.D=Inf;
                    if abs(Zfwt(2)/Zfwt(1))>exp(zxi(2)-zxi(1)) %determine if Dpsi is finite, i.e. fwt\sim\xi^(1+a>0) when xi->0
                        wp.D=(wp.ompeak/2)*sum(exp(-zxi).*conj(Zfwt).*dbxi);
                    end
                end %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                
                %Calculate the cumulative integrals
                CS=CS0+cumsum(Zfwt.*dbxi); CS=[CS0;CS(:)]/CS(end); CS=abs(CS);
                ICS=cumsum(flipud(Zfwt.*dbxi)); ICS=ICS(end:-1:1); ICS=[ICS(:);0]/(ICS(1)+CS0); ICS=abs(ICS);
                
                %Estimate epsilon-supports
                xid=find(CS(1:end-1)<racc/2 & CS(2:end)>=racc/2,1,'first');
                if isempty(xid), wp.xi1e=exp(bxi(1));
                else
                    a1=CS(xid)-racc/2; a2=CS(xid+1)-racc/2;
                    wp.xi1e=exp(bxi(xid)-a1*(bxi(xid+1)-bxi(xid))/(a2-a1));
                end
                xid=find(ICS(1:end-1)>=racc/2 & ICS(2:end)<racc/2,1,'last');
                if isempty(xid), wp.xi2e=exp(bxi(end));
                else
                    a1=ICS(xid)-racc/2; a2=ICS(xid+1)-racc/2;
                    wp.xi2e=exp(bxi(xid)-a1*(bxi(xid+1)-bxi(xid))/(a2-a1));
                end
                xid=find(CS(1:end-1)<0.25 & CS(2:end)>=0.25,1,'first');
                if isempty(xid), wp.xi1h=exp(bxi(1));
                else
                    a1=CS(xid)-0.25; a2=CS(xid+1)-0.25;
                    wp.xi1h=exp(bxi(xid)-a1*(bxi(xid+1)-bxi(xid))/(a2-a1));
                end
                xid=find(ICS(1:end-1)>=0.25 & ICS(2:end)<0.25,1,'last');
                if isempty(xid), wp.xi2h=exp(bxi(end));
                else
                    a1=ICS(xid)-0.25; a2=ICS(xid+1)-0.25;
                    wp.xi2h=exp(bxi(xid)-a1*(bxi(xid+1)-bxi(xid))/(a2-a1));
                end
                
            end
            
            %Return to the time-domain form (we needed to estimate [wp.ompeak] if not known)
            vfun=@(u)conj(twf(u).*exp(-1i*wp.ompeak*u)); %demodulated wavelet in the time domain
            xp=wp.tpeak; lim1=wp.t1; lim2=wp.t2;
            
            [QQ,wflag,xx,ss]=sqeps(vfun,xp,[lim1,lim2],racc,MIC,...
                [-8*(2*pi*fmax/wp.ompeak)*L/fs,8*(2*pi*fmax/wp.ompeak)*L/fs]); %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            wp.t1e=ss(1,1); wp.t2e=ss(1,2); wp.t1h=ss(2,1); wp.t2h=ss(2,2);
            if wflag==1 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                fprintf('The time-domain wavelet function is not well-behaved (e.g. decays very slowly as time tends to +/- infinity).\n');
                fprintf('The integration might be not accurate (and therefore e.g. cone-of-influence borders).\n');
                fprintf(2,'----------------------------------------------------------------------------------------------------\n');
            end
            
        end
        
    end
%======================================================================================================================


end


%============ Function for accurate integration and estimation ============
%================= of epsilon-supports (needed in parcalc) ================
%INPUT:
%[vfun] - function @(x) to be integrated
%[xp] - position of the peak
%[lims]=[lim1,lim2] - limits where [vfun] is defined
%[racc] - accuracy (i.e. epsilon)
%[MIC] - maximum number of intervals in quadgk(...)
%[nlims] - the maximum limits dictated by the signal
%OUTPUT:
%[QQ]=[Q1,Q2;q1m,q2m] - integrals of the function:
%     [Q1,Q2] - over [lim1,xp] and [xp,lim2]
%     [q1m,q2m] - over [x1m,xp] and [xp,x2m]
%[wflag]=0|1 - was integration fully successful (0) or not (1)
%[xx]=[x1h,x2h;x1e,x2e;x1m,x2m;x1n,x2n] - different important points:
%     [x1e,x2e] - interval where [vfun] is larger than 10^(-8) of maximum
%     [x1h,x2h] - where [vfun] first drops below 0.5 of its maximum
%     [x1m,x2m] - interval over which integration with quadgk is precise
%                 (just =[x1e,x2e] if all is precise, i.e. wflag=0)
%     [x1n,x2n] - limits of overflow (i.e. when function returns NaN due to
%                 numerics, as for exp(-x^2)*exp(x) for large x); if no
%                 overflow, then [x1n,x2n]=[lim1,lim2]
%optional:
%[ss]=[s1h,s2h;s1e,s2e] - supports of the function over its variable:
%     [s1e,s2e] - epsilon-support
%     [s1h,s2h] - 0.5-support
function [QQ,wflag,xx,ss]=sqeps(vfun,xp,lims,racc,MIC,nlims)

wflag=0; %indicates are there any problems with integration
ctol=max([racc/1000,10^(-12)]); %numerical accuracy
lim1=lims(1); lim2=lims(2);
nlim1=nlims(1); nlim2=nlims(2);

kk=1; shp=0; %check for function to be well-defined on the peak, otherwise introduce some shift [shp]
while (~isfinite(vfun(xp+shp)) || isnan(vfun(xp+shp))), shp=kk*10^(-14); kk=kk*(-2); end
vmax=vfun(xp+shp);

%determine x after for which function first crosses 0.5 of its maximal value (if such point exists)
if isfinite(lim1), tx1=lim1-0.01*(lim1-xp); qv1=abs(vfun(tx1)/vmax); else qv1=NaN; end
if qv1<0.5, x1h=fzero(@(u)abs(vfun(u)/vmax)-0.5,[xp+shp,tx1],optimset('Display','off'));
elseif isnan(qv1), x1h=fzero(@(u)abs(vfun(xp-abs(u))/vmax)-0.5,shp,optimset('Display','off')); x1h=xp-abs(x1h);
else x1h=xp+(lim1-xp)/100;
end
if isfinite(lim2), tx2=lim2-0.01*(lim2-xp); qv2=abs(vfun(tx2)/vmax); else qv2=NaN; end
if qv2<0.5, x2h=fzero(@(u)abs(vfun(u)/vmax)-0.5,[xp+shp,tx2],optimset('Display','off'));
elseif isnan(qv2), x2h=fzero(@(u)abs(vfun(xp+abs(u))/vmax)-0.5,shp,optimset('Display','off')); x2h=xp+abs(x2h);
else x2h=xp+(lim2-xp)/100;
end
if isnan(x1h)
    x1h=fminsearch(@(u)abs(abs(vfun(xp-abs(u))/vmax)-0.5),shp,optimset('TolFun',10^(-12)));
    x1h=xp-abs(x1h)/100;
end
if isnan(x2h)
    x2h=fminsearch(@(u)abs(abs(vfun(xp+abs(u))/vmax)-0.5),shp,optimset('TolFun',10^(-12)));
    x2h=xp+abs(x2h)/100;
end
    

%determine x after which function is always below 10^(-8) of its maximal value (if such point exists), and not only the first crossing
if isfinite(lim1), tx1=lim1-0.01*(lim1-xp); qv1=(abs(vfun(tx1))+abs(vfun((tx1+lim1)/2))+abs(vfun((tx1+3*lim1)/4)))/abs(vmax); else qv1=NaN; end
if qv1<10^(-8)/3
    x1e=fzero(@(u)abs(vfun(u)/vmax)+abs(vfun((u+lim1)/2)/vmax)+abs(vfun((u+3*lim1)/4)/vmax)-10^(-8),[xp+shp,tx1],optimset('Display','off'));
elseif isnan(qv1)
    x1e=fzero(@(u)abs(vfun(xp-abs(u))/vmax)+abs(vfun(xp-sqrt(3)*abs(u))/vmax)+...
        abs(vfun(xp-sqrt(5)*abs(u))/vmax)-10^(-8),shp,optimset('Display','off')); %to avoid selecting first minima for oscillating functions
    x1e=xp-abs(x1e);
else x1e=xp+(lim1-xp)/2;
end
if isfinite(lim2), tx2=lim2-0.01*(lim2-xp); qv2=(abs(vfun(tx2))+abs(vfun((tx2+lim2)/2))+abs(vfun((tx2+3*lim2)/4)))/abs(vmax); else qv2=NaN; end
if qv2<10^(-8)
    x2e=fzero(@(u)abs(vfun(u)/vmax)+abs(vfun((u+lim2)/2)/vmax)+abs(vfun((u+3*lim2)/4)/vmax)-10^(-8),[xp+shp,tx2],optimset('Display','off'));
elseif isnan(qv2)
    x2e=fzero(@(u)abs(vfun(xp+abs(u))/vmax)+abs(vfun(xp+sqrt(3)*abs(u))/vmax)+...
        abs(vfun(xp+sqrt(5)*abs(u))/vmax)-10^(-8),shp,optimset('Display','off'));
    x2e=xp+abs(x2e);
else x2e=xp+(lim2-xp)/2;
end
if isnan(x1e)
    x1e=fminsearch(@(u)abs(abs(vfun(x1h-abs(u))/vmax)-10^(-8)),0,optimset('TolFun',10^(-12)));
    x1e=x1h-abs(x1e); lim1=x1e; wflag=1;
end
if isnan(x2e)
    x2e=fminsearch(@(u)abs(abs(vfun(x2h+abs(u))/vmax)-10^(-8)),0,optimset('TolFun',10^(-12)));
    x2e=x2h+abs(x2e); lim2=x2e; wflag=1;
end


wstate=warning('off','all'); %disable warnings


%=============== Integrate given function to find Q1 and Q2 ===============
Q1=0; Q2=0;
%--------------------------------------------------------------------------
[qv,eb]=quadgk(@(u)vfun(u),xp,x1e,'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
[qv,eb]=quadgk(@(u)vfun(u),xp,x1e,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q1+qv)),'RelTol',0);
x1m=x1e; q1m=qv;
if abs(eb/(Q1+qv))>ctol
    wflag=1;
    [qv,eb]=quadgk(@(u)vfun(u),xp,x1h,'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
    [qv,eb]=quadgk(@(u)vfun(u),xp,x1h,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q1+qv)),'RelTol',0);
    x1m=x1h; q1m=qv;
    if abs(eb/(Q1+qv))<=ctol
        while 1==1
            [qv,eb]=quadgk(@(u)vfun(u),xp,max([xp+(x1m-xp)*2,x1e]),'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
            [qv,eb]=quadgk(@(u)vfun(u),xp,max([xp+(x1m-xp)*2,x1e]),'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q1+qv)),'RelTol',0);
            if abs(eb/(Q1+qv))>ctol, break; end
            x1m=max([xp+(x1m-xp)*2,x1e]); q1m=qv;
        end
    else
        while abs(eb/(Q1+qv))>ctol
            x1m=xp+(x1m-xp)/2;
            [qv,eb]=quadgk(@(u)vfun(u),xp,x1m,'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
            [qv,eb]=quadgk(@(u)vfun(u),xp,x1m,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q1+qv)),'RelTol',0);
            q1m=qv;
        end
    end
    [qv,eb]=quadgk(@(u)vfun(u),x1m,x1e,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q1+qv)),'RelTol',0);
    if abs(eb)<0.5*abs(qv), Q1=Q1+qv; end
end
Q1=Q1+q1m;
if wflag==0
    [qv,eb]=quadgk(@(u)vfun(u),x1e,lim1,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q1),'RelTol',0);
    if ~isfinite(lim1) && isnan(qv) %to avoid overflow
        lim1=x1e; while ~isnan(vfun(2*lim1)), lim1=2*lim1; end
        [qv,eb]=quadgk(@(u)vfun(u),x1e,lim1,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q1),'RelTol',0);
    end
    if abs(eb/Q1)>ctol
        wflag=1;
        [qv,eb]=quadgk(@(u)vfun(u),x1e,max([min([8*x1e,nlim1]),lim1]),'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q1),'RelTol',0);
    end
    if abs(eb)<0.5*abs(qv), Q1=Q1+qv; end
end
Q1=-Q1;
%--------------------------------------------------------------------------
[qv,eb]=quadgk(@(u)vfun(u),xp,x2e,'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
[qv,eb]=quadgk(@(u)vfun(u),xp,x2e,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q2+qv)),'RelTol',0);
x2m=x2e; q2m=qv;
if abs(eb/(Q2+qv))>ctol
    wflag=1;
    [qv,eb]=quadgk(@(u)vfun(u),xp,x2h,'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
    [qv,eb]=quadgk(@(u)vfun(u),xp,x2h,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q2+qv)),'RelTol',0);
    x2m=x2h; q2m=qv;
    if abs(eb/(Q2+qv))<=ctol
        while 1==1
            [qv,eb]=quadgk(@(u)vfun(u),xp,min([xp+(x2m-xp)*2,x2e]),'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
            [qv,eb]=quadgk(@(u)vfun(u),xp,min([xp+(x2m-xp)*2,x2e]),'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q2+qv)),'RelTol',0);
            if abs(eb/(Q2+qv))>ctol, break; end
            x2m=min([xp+(x2m-xp)*2,x2e]); q2m=qv;
        end
    else
        while abs(eb/(Q2+qv))>ctol
            x2m=xp+(x2m-xp)/2;
            [qv,eb]=quadgk(@(u)vfun(u),xp,x2m,'MaxIntervalCount',MIC,'AbsTol',0,'RelTol',0.1*ctol);
            [qv,eb]=quadgk(@(u)vfun(u),xp,x2m,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q2+qv)),'RelTol',0);
            q2m=qv;
        end
    end
    [qv,eb]=quadgk(@(u)vfun(u),x2m,x2e,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*(Q2+qv)),'RelTol',0);
    if abs(eb)<0.5*abs(qv), Q2=Q2+qv; end
end
Q2=Q2+q2m;
if wflag==0
    [qv,eb]=quadgk(@(u)vfun(u),x2e,lim2,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q2),'RelTol',0);
    if ~isfinite(lim2) && isnan(qv) %to avoid overflow
        lim2=x2e; while ~isnan(vfun(2*lim2)), lim2=2*lim2; end
        [qv,eb]=quadgk(@(u)vfun(u),x2e,lim2,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q2),'RelTol',0);
    end
    if abs(eb/Q2)>ctol
        wflag=1;
        [qv,eb]=quadgk(@(u)vfun(u),x2e,min([max([8*x2e,nlim2]),lim2]),'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q2),'RelTol',0);
    end
    if abs(eb)<0.5*abs(qv), Q2=Q2+qv; end
end

QQ=[Q1,Q2;-q1m,q2m]; xx=[x1e,x2e;x1h,x2h;x1m,x2m;lim1,lim2];
%--------------------------------------------------------------------------

%======== Find 0.5- and \epsilon-supports of the specified function =======
if nargout>3
    Q=Q1+Q2;
    s1e=fz(racc/2); s2e=fz(1-racc/2); s1h=fz(0.5/2); s2h=fz(1-0.5/2);
    ss=[s1e,s2e;s1h,s2h];
end
%--------------------------------------------------------------------------

warning(wstate); %restore the warning settings


    %function for finding the epsilon-supports
    function x0=fz(zv)
        if zv<abs(Q1/Q), cx1=x1m; cq1=Q1+q1m; ra=exp(-1/2); rb=exp(1/2);
        else cx1=x2m; cq1=Q1+q2m; ra=exp(1/2); rb=exp(-1/2); end
        if abs(1-abs((Q-cq1)/Q))<zv
            while 1==1
                nx=xp+ra*(cx1-xp);
                if nx<lim1, nx=(cx1+lim1)/2; end
                if nx>lim2, nx=(cx1+lim2)/2; end
                [pv,~]=quadgk(@(u)vfun(u),cx1,nx,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q),'RelTol',0);
                if abs(1-abs((Q-cq1-pv)/Q))>=zv, cx2=nx; cq2=cq1+pv; break; end
                cx1=nx; cq1=cq1+pv;
            end
        else
            while 1==1
                nx=xp+rb*(cx1-xp);
                if nx<lim1, nx=(cx1+lim1)/2; end
                if nx>lim2, nx=(cx1+lim2)/2; end
                [pv,~]=quadgk(@(u)vfun(u),cx1,nx,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q),'RelTol',0);
                if abs(1-abs((Q-cq1-pv)/Q))<zv, cx2=nx; cq2=cq1+pv; break; end
                cx1=nx; cq1=cq1+pv;
            end
        end
        [pv,~]=quadgk(@(u)vfun(u),cx1,(cx1+cx2)/2,'MaxIntervalCount',MIC,'AbsTol',0.1*abs(ctol*Q),'RelTol',0);
        qfun=@(x)1-abs((Q-(cq1+pv+quadgk(@(u)vfun(u),(cx1+cx2)/2,x,'MaxIntervalCount',round(MIC/10),'AbsTol',0.5*abs(ctol*Q),'RelTol',0)))/Q);
        x0=fzero(@(x)abs(qfun(x))-zv,[cx1,cx2],optimset('TolX',10^(-14)));
    end


end



%================== Predictive padding function ===========================

% fsig=fcast(sig,fs,NP,fint,Optional:MaxOrder,w)
% - given signal [sig] sampled at [fs] Hz, uses DFT and weighted least
% squares (with weights specified by [w], which should be of the same
% length as [sig]), finds main sinusoidal components present in the signal
% and uses them to predict signal for NP consequtive time-steps, which is
% returned in [fsig]; [fint]=[fmin,fmax] specifies the allowable frequency
% range for sinusoids: tones with frequencies out of it are not continued
% for prediction. The number of sinusoids is determined using Bayesian
% (Schwarz) information criterion, but it cannot exceed [MaxOrder].

function fsig = fcast(sig,fs,NP,fint,varargin)

L=length(sig); T=L/fs; t=0:1/fs:T-1/fs;
MaxOrder=L; if nargin>3 && ~isempty(varargin{1}), MaxOrder=varargin{1}; end
w=[]; if nargin>4, w=varargin{2}(:); end, rw=sqrt(w);
FTol=1/T/100; %accuracy of frequencies determination
Nq=ceil((L+1)/2); ftfr=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L;
Y=sig(:); if ~isempty(rw), Y=rw.*Y; end
orstd=std(Y); %original signal's standard deviation
rr=(1+sqrt(5))/2; %multiplier to easily use golden section search afterwards

v=zeros(L,1); ic=v; frq=v; amp=v; phi=v; itn=0;
while itn<MaxOrder
    itn=itn+1;
    aftsig=abs(fft(Y)); [~,imax]=max(aftsig(2:Nq-1)); imax=imax+1;
    
    %========================= Forward search =============================
    %Find two points between which minima occurs
    nf=ftfr(imax);
    FM=[ones(L,1),cos(2*pi*nf*t(:)),sin(2*pi*nf*t(:))];
    if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
    nb=FM\Y; nerr=std(Y-FM*nb);
    
    df=FTol; perr=Inf; %starting frequency step and 'previous error'
    while nerr<perr
        if abs(nf-fs/2+FTol)<eps, break; end
        pf=nf; perr=nerr; pb=nb;
        nf=min([pf+df,fs/2-FTol]);
        FM=[ones(L,1),cos(2*pi*nf*t(:)),sin(2*pi*nf*t(:))];
        if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
        nb=FM\Y; nerr=std(Y-FM*nb);
        df=df*rr;
    end
    
    %Use golden section search to find exact minimum
    if nerr<perr
        cf=[nf,nf,nf]; cerr=[nerr,nerr,nerr]; cb=[nb(:),nb(:),nb(:)];
    elseif abs(nf-pf-FTol)<eps
        cf=[pf,pf,pf]; cerr=[perr,perr,perr]; cb=[pb(:),pb(:),pb(:)];
    else
        cf=[0,pf,nf]; cerr=[0,perr,nerr]; cb=[zeros(length(pb),1),pb(:),nb(:)];
        cf(1)=pf-df/rr/rr;
        FM=[ones(L,1),cos(2*pi*cf(1)*t(:)),sin(2*pi*cf(1)*t(:))];
        if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
        cb(:,1)=FM\Y; cerr(1)=std(Y-FM*cb(:,1));
    end
    while (cf(2)-cf(1)>FTol && cf(3)-cf(2)>FTol)
        tf=cf(1)+cf(3)-cf(2);
        FM=[ones(L,1),cos(2*pi*tf*t(:)),sin(2*pi*tf*t(:))];
        if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
        tb=FM\Y; terr=std(Y-FM*tb);
        if terr<cerr(2) && tf<cf(2), cf=[cf(1),tf,cf(2)]; cerr=[cerr(1),terr,cerr(2)]; cb=[cb(:,1),tb(:),cb(:,2)]; end
        if terr<cerr(2) && tf>cf(2), cf=[cf(2),tf,cf(3)]; cerr=[cerr(2),terr,cerr(3)]; cb=[cb(:,2),tb(:),cb(:,3)]; end
        if terr>cerr(2) && tf<cf(2), cf=[tf,cf(2),cf(3)]; cerr=[terr,cerr(2),cerr(3)]; cb=[tb(:),cb(:,2),cb(:,3)]; end
        if terr>cerr(2) && tf>cf(2), cf=[cf(1),cf(2),tf]; cerr=[cerr(1),cerr(2),terr]; cb=[cb(:,1),cb(:,2),tb(:)]; end
    end
    fcf=cf(2); fcb=cb(:,2); fcerr=cerr(2); %forward values
    
    %========================= Backward search ============================
    %Find two points between which minima occurs
    nf=ftfr(imax);
    FM=[ones(L,1),cos(2*pi*nf*t(:)),sin(2*pi*nf*t(:))];
    if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
    nb=FM\Y; nerr=std(Y-FM*nb);
    
    df=FTol; perr=Inf; %starting frequency step and 'previous error'
    while nerr<perr
        if abs(nf-FTol)<eps, break; end
        pf=nf; perr=nerr; pb=nb;
        nf=max([pf-df,FTol]);
        FM=[ones(L,1),cos(2*pi*nf*t(:)),sin(2*pi*nf*t(:))];
        if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
        nb=FM\Y; nerr=std(Y-FM*nb);
        df=df*rr;
    end
    
    %Use golden section search to find exact minimum
    if nerr<perr
        cf=[nf,nf,nf]; cerr=[nerr,nerr,nerr]; cb=[nb(:),nb(:),nb(:)];
    elseif abs(nf-pf-FTol)<eps
        cf=[pf,pf,pf]; cerr=[perr,perr,perr]; cb=[pb(:),pb(:),pb(:)];
    else
        cf=[nf,pf,0]; cerr=[nerr,perr,0]; cb=[nb(:),pb(:),zeros(length(pb),1)];
        cf(3)=pf+df/rr/rr;
        FM=[ones(L,1),cos(2*pi*cf(1)*t(:)),sin(2*pi*cf(1)*t(:))];
        if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
        cb(:,3)=FM\Y; cerr(3)=std(Y-FM*cb(:,3));
    end
    while (cf(2)-cf(1)>FTol && cf(3)-cf(2)>FTol)
        tf=cf(1)+cf(3)-cf(2);
        FM=[ones(L,1),cos(2*pi*tf*t(:)),sin(2*pi*tf*t(:))];
        if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
        tb=FM\Y; terr=std(Y-FM*tb);
        if terr<cerr(2) && tf<cf(2), cf=[cf(1),tf,cf(2)]; cerr=[cerr(1),terr,cerr(2)]; cb=[cb(:,1),tb(:),cb(:,2)]; end
        if terr<cerr(2) && tf>cf(2), cf=[cf(2),tf,cf(3)]; cerr=[cerr(2),terr,cerr(3)]; cb=[cb(:,2),tb(:),cb(:,3)]; end
        if terr>cerr(2) && tf<cf(2), cf=[tf,cf(2),cf(3)]; cerr=[terr,cerr(2),cerr(3)]; cb=[tb(:),cb(:,2),cb(:,3)]; end
        if terr>cerr(2) && tf>cf(2), cf=[cf(1),cf(2),tf]; cerr=[cerr(1),cerr(2),terr]; cb=[cb(:,1),cb(:,2),tb(:)]; end
    end
    bcf=cf(2); bcb=cb(:,2); bcerr=cerr(2); %backward values
    
    %================ Assigning values and subtracting ====================
    if fcerr<bcerr, cf=fcf; cb=fcb; cerr=fcerr;
    else cf=bcf; cb=bcb; cerr=bcerr; end
    
    frq(itn+1)=cf; amp(itn+1)=sqrt(cb(2)^2+cb(3)^2); phi(itn+1)=atan2(-cb(3),cb(2));
    amp(1)=amp(1)+cb(1); v(itn)=cerr;
    
    FM=[ones(L,1),cos(2*pi*cf*t(:)),sin(2*pi*cf*t(:))];
    if ~isempty(rw), FM=FM.*(rw*ones(1,3)); end
    Y=Y-FM*cb;
    
    %===================== Information criterion ==========================
    CK=3*itn+1;
    cic=L*log(cerr)+CK*log(L); %Bayesian (Schwarz) IC (the best)
    %cic=L*log(cerr)+2*CK*L/(L-CK-1); %Akaike IC (corrected)
    %cic=cerr*(L+CK)/(L-CK); %Forward Prediction Error (FPE)
    
    ic(itn)=cic;
    if v(itn)/orstd<2*eps, break; end
    if itn>2 && cic>ic(itn-1) && ic(itn-1)>ic(itn-2), break; end
end

frq=frq(1:itn+1); amp=amp(1:itn+1); phi=phi(1:itn+1); v=v(1:itn); ic=ic(1:itn);
fsig=zeros(NP,1); nt=(T:1/fs:T+(NP-1)/fs)';
if size(sig,2)>size(sig,1), fsig=fsig'; nt=nt'; end
for k=1:length(frq)
    if frq(k)>fint(1) && frq(k)<fint(2)
        fsig=fsig+amp(k)*cos(2*pi*frq(k)*nt+phi(k));
    else
        fsig=fsig+amp(k)*cos(2*pi*frq(k)*(T-1/fs)+phi(k));
    end
end

end


%=============== Interpolation function for plotting ======================

% nZ=aminterp(X,Y,Z,XI,YI,method)
% - the same as interp2, but uses different interpolation types,
% maximum-based ([method]='max') or average-based ([method]='avg'), where
% [nZ] at each point [XI,YI] represents maximum or average among values of
% [Z] corresponding to the respective quadrant in [X,Y]-space which
% includes point [XI,YI];
% X,Y,XI,YI are all linearly spaced vectors, and the lengths of XI,YI
% should be the same or smaller than that of X,Y; Z should be real,
% ideally positive; X and Y correspond to the 2 and 1 dimension of Z, as
% always. 

function nZ = aminterp(X,Y,Z,XI,YI,method)

%Interpolation over X
nZ=zeros(size(Z,1),length(XI))*NaN;
xstep=mean(diff(XI));
xind=1+floor((1/2)+(X-XI(1))/xstep); xind=xind(:);
xpnt=[0;find(xind(2:end)>xind(1:end-1));length(xind)];
if strcmpi(method,'max')
    for xn=1:length(xpnt)-1
        xid1=xpnt(xn)+1; xid2=xpnt(xn+1);
        nZ(:,xind(xid1))=max(Z(:,xid1:xid2),[],2);
    end
else
    for xn=1:length(xpnt)-1
        xid1=xpnt(xn)+1; xid2=xpnt(xn+1);
        nZ(:,xind(xid1))=mean(Z(:,xid1:xid2),2);
    end
end

Z=nZ;

%Interpolation over Y
nZ=zeros(length(YI),size(Z,2))*NaN;
ystep=mean(diff(YI));
yind=1+floor((1/2)+(Y-YI(1))/ystep); yind=yind(:);
ypnt=[0;find(yind(2:end)>yind(1:end-1));length(yind)];
if strcmpi(method,'max')
    for yn=1:length(ypnt)-1
        yid1=ypnt(yn)+1; yid2=ypnt(yn+1);
        nZ(yind(yid1),:)=max(Z(yid1:yid2,:),[],1);
    end
else
    for yn=1:length(ypnt)-1
        yid1=ypnt(yn)+1; yid2=ypnt(yn+1);
        nZ(yind(yid1),:)=mean(Z(yid1:yid2,:),1);
    end
end


end
