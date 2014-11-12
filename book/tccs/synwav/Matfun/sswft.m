%Version 1.01 stable
%**************************************************************************
%*************** Synchrosqueezed Windowed Fourier Transform ***************
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
% [SWFT,freq,Optional:wopt,WFT,wfreq,IFR]=sswft(sig,fs,Optional:'PropertyName',PropertyValue)
%
% INPUT:
% sig - signal for which to calculate SWFT
% fs  - sampling frequency of the signal
% 
% Properties: ({{...}} denotes default)
% /same as for wft:/
% 'fmin':value (default = 0)
%            minimal frequency for which to calculate SWFT
% 'fmax':value (default = fs/2, i.e. the Nyquist frequency)
%            maximal frequency for which to calculate SWFT
% 'fstep':{{'auto'}}|'auto-NB'|value
%            frequency step, which determines frequency discretization,
%            so that the next frequency equals previous one plus [fstep];
%            when set to 'auto-NB' (e.g. 'auto-20') determines [fstep]
%            automatically as described in [1], so that it equals
%            1/NB of the frequency region containing 50% of the
%            window function; the default 'auto' is equivalent to 'auto-10'
% 'Window':{{'Gaussian'}}|'Hann'|'Blackman'|'Exp'|'Kaiser-a'|{@(xi)fwt(xi),[xi1,xi2],@(t)twf(t),[t1,t2]}
%            window used in SWFT calculation, for a list of all supported
%            names and their properties see Appendix E in [1]. However, you
%            can use any window by specifying its frequency domain form,
%            i.e. window FT (defined by function [fwt], the argument of
%            which is cyclic frequency) and/or time-domain form [twf] (the
%            argument is time) together with corresponding full supports;
%            only [fwt] or [twf] is enough, so if one of them is not known
%            just put empty field [] for it and its support (but it is
%            better to specify both [fwt] and [twf] when available). Thus,
%            the Gaussian window with [f0=1] can be alternatively defined as
%            {@(xi)exp(-(1/2)*(xi.^2)),[-Inf,Inf],@(t)(1/sqrt(2*pi))*exp(-(1/2)*(t.^2)),[-Inf,Inf]},
%            while Hann window with [f0=1] -- e.g. as {[],[],@(t)(1+cos(2*pi*t/4.4))/2,[-2.2,2.2]}
% 'f0':value (default = 1)
%            window resolution parameter, which determines the tradeoff
%            between the time and frequency resolutions: the higher it is,
%            the closer in frequency components can be resolved in (S)WFT,
%            but the slower time-variations, e.g. amplitude/frequency
%            modulation, can be reliably represented; for the way it is
%            introduced for each window see Appendix E in [1], while if the
%            window is user-defined in terms of its function in frequency
%            and/or time (see 'Window' property), then 'f0' obviously does
%            not influence anything.
% 'Preprocess':{{'on'}}|'off'
%            perform or not an initial signal preprocessing, which consists
%            of subtracting 3rd order polynomial fit and then bandpassing
%            the signal in the band of interest [fmin,fmax], for which SWFT
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
%            sswft(signal(n1:n2),fs,'Padding',{signal(1:n1-1),signal(n2+1:end)});
%            the same applies if you can predict signal out of the time-limits
% 'RelTol':value (default = 0.01) (in [1] commonly referred as \epsilon)
%            relative tolerance (e.g. 0.01 means 1%), which specifies cone
%            of influence for (S)WFT (i.e. range of (S)WFT coefficients which are
%            determined up to this accuracy in respect of boundary errors);
%            determines also the minimal number of values to pad signal
%            with, so that relative contribution of effects of implicit
%            periodic signal continuation due to convolution in frequency
%            domain is smaller [RelTol], see [1] for details; additonally,
%            determines for how wider frequency region should be calculated
%            WFT so that SWFT coefficients in the [fmin,fmax] are
%            determined with a specified accuracy in relation to boundary
%            effects in frequency (i.e. at least (1-RelTol) of the full
%            area under the WFT peaks is acocunted for, see [1]).
% 'Plot':{{'off'}} - do not plot anything
%          'amp'   - plots SWFT amplitude (i.e. its absolute value) together
%                    with line denoting the cone of influence
%          'amp+'  - additionally shows time-averaged SWFT amplitude
%          'amp++' - additionally shows 95% range of SWFT amplitude
%          'pow','pow+','pow++' - the same but for the SWFT power (i.e. its 
%                                 squared modulus)
%           IMPORTANT: to avoid plotting huge data (in which case it might
%           be very slow to render and modify the figure, and MatLab can
%           even crash), the plotted SWFT is resampled to have no more than
%           few data points displayed per pixel (for the current screen
%           resolution). Unless 'Display' is 'off', it will always notify if
%           the SWFT size exceeds the number of pixels on the plot and the
%           resamling is therefore performed. Note, that in this case considerable
%           zooming in of the parts of displayed plots might not show the full
%           structure of the original SWFT. If one wants to investigate the
%           resultant plot in fine details (and not only see how it looks),
%           then the original, full SWFT might be displayed without resampling
%           by adding 'wr' to the end of this option, e.g. 'ampwr' or 'pow++wr'.
% 'Display':{{'on'}} - displays all relevant information about progress etc.
%           'notify' - displays only information if something went wrong
%            'off'   - does not display anything
% 'CutEdges':{{'off'}}|'on'
%           determine should (S)WFT coefficients be set to NaNs out of the cone
%           of influence (determined for the current 'RelTol'); if you wish
%           to analyze only (S)WFT within the cone of influence, which is
%           recommended if you want to estimate e.g. only the time-averaged
%           quantities, set it to 'on'.
% /additional to the properties of wft:/
% 'CutFreq':{{'on'}}|'off'
%           as discussed in [1], WFT used for synchrosqueezing should
%           be calculated in a broader frequency range to assure the
%           specified relative error 'RelTol' in SWFT in a given range
%           [fmin,fmax]; if WFT is requested as output, 'CutFreq' specifies
%           should it be truncated to the frequency band of SWFT
% 'Threshold':value|'MAD' (default = 0, i.e. no denoising)
%           before performing synchrosqueezing, one can additionally
%           denoise the calculated WFT by hard thresholding, i.e. set
%           WFT coefficients lower than the specified 'Threshold' to zero;
%           when set to 'MAD', the value of the threshold is determined as
%           the median between median deviations of the WFT amplitude
%           (around its median value) at frequencies from [(fmax+fmin)/2]
%           to [fmax], multiplied on 1.4826*sqrt(2*log(L)) (L is the signal
%           length in samples), as proposed in arXiv:1105.0010 for WT
%           (for 'MAD' it is better to calculate WFT for a full range,
%           i.e. set 'fmax' to fs/2 or just leave as default; but even in
%           this case the calculated threshold might be not appropriate in
%           some cases, e.g. when there is some AM/FM component at
%           frequencies between [(fmax+fmin)/2] and [fmax])
%
%
% OUTPUT:
% SWFT  - synchrosqueezed windowed Fourier transform of the signal [sig],
%         with rows corresponding to frequencies and columns - to time;
%         represents FNxL matrix, where FN is the number of frequencies
%         and L is the length of the signal in samples
% freq  - frequencies corresponding to rows of SWFT
% wopt  - structure with all parameters of the window and simulation
% WFT   - windowed Fourier transform, from which SWFT was constructed,
%         with rows corresponding to frequencies and columns - to time;
%         represents WNxL matrix, where WN is the number of frequencies for
%         which it was computed (which span broader frequency range than
%         [freq] for SWFT, unless 'CutFreq' is set to 'on')
% wfreq - frequencies corresponding to rows of WFT
% IFR   - the Instantaneous Frequency Representation corresponding to the
%         current WFT ([IFR] is \nu_G/2\pi in notations of [1]), i.e.
%         the rate of phase growth of WFT coefficients at each time and 
%         frequency, i.e. in fact WFT instantaneous frequency;
%         obviously, has the same size as returned WFT.
%
%-------------------------------Examples-----------------------------------
%
% [SWFT,freq]=sswft(sig,100,'fmin',0.1,'fmax',10,'f0',2)
% given a signal [sig] sampled at 100 Hz, returns its SWFT [SWFT] based on a
% Gaussian window (default) with [f0=2] and calculated for frequencies from
% 0.1 to 10 Hz, returned in [freq].
%
%-----------------------Additional possibilities---------------------------
%
% One can also pass the structure with the properties as a third input
% argument instead of specifying them by pairs, e.g.
% opts=struct; opts.f0=2; opts.fmax=12; opts.Padding='predictive';
% [SWFT,freq]=sswft(sig,fs,opts);
% You can also add further parameters in the usual way, e.g.
% [SWFT,freq]=sswft(sig,fs,opts,'Display','off');
% If you add a parameter contained in [opts], it will be ovewritten, e.g.
% [SWFT,freq]=sswft(sig,fs,opts,'fmax',10);
% will change 'fmax' property from 12 specified in [opts.fmax] to 10.
%
% When one needs to calculate SWFT using the same window and simulation
% parameters, but for different signals, it is a good idea to do it like
% [SWFT1,freq1,wopt]=sswft(sig1,fs1,...(parameters));
% [SWFT2,freq2]=sswft(sig2,fs2,wopt); [SWFT3,freq3]=sswft(sig3,fs3,wopt); ...
% This will also avoid recalculating window parameters each time, thus
% slightly speeding up the computations.
%
% The same can be done if you calculated SWFT of a signal but want to obtain
% it using some different parameter (e.g. Padding), for example
% [SWFT1,freq1,wopt]=sswft(sig,fs,'Padding',0,...(other parameters));
% [SWFT2,freq2]=sswft(sig,fs,wopt,'Padding','predictive');
% Unless you overwrite 'f0', 'Window' or 'RelTol' properties, the window
% parameters will not be recalculated for the second time.
%
%------------------------------Changelog-----------------------------------
%
% v1.01:
% - changed estimation of WFT instantaneous frequencies [IFR] from direct
%   numerical differentiation to differentiation in the frequency domain
%   (see [1]); however, although this is in principle more accurate with
%   respect to time discretization, I have not noticed any differences...
%
%--------------------------------------------------------------------------


function [SWFT,freq,varargout] = sswft(signal,fs,varargin)

L=length(signal); signal=signal(:);

%Default parameters
Window='Gaussian'; f0=1;
fmin=0; fmax=fs/2;
fstep='auto';
PadMode='predictive';
RelTol=0.01;
Preprocess='on';
DispMode='on';
PlotMode='off';
CutEdges='off';
CutFreq='on';
Threshold=0;
%Update if user defined
vst=1; recflag=1;
if nargin>2 && isstruct(varargin{1})
    opts=varargin{1}; vst=2;
    if isfield(opts,'Window'), Window=opts.Window; end
    if isfield(opts,'f0'), f0=opts.f0; end
    if isfield(opts,'fmin'), fmin=opts.fmin; end
    if isfield(opts,'fmax'), fmax=opts.fmax; end
    if isfield(opts,'fstep'), fstep=opts.fstep; end
    if isfield(opts,'Padding'), PadMode=opts.Padding; end
    if isfield(opts,'RelTol'), RelTol=opts.RelTol; end
    if isfield(opts,'Preprocess'), Preprocess=opts.Preprocess; end
    if isfield(opts,'Plot'), PlotMode=opts.Plot; end
    if isfield(opts,'Display'), DispMode=opts.Display; end
    if isfield(opts,'CutEdges'), CutEdges=opts.CutEdges; end
    if isfield(opts,'CutFreq'), CutFreq=opts.CutFreq; end
    if isfield(opts,'Threshold'), Threshold=opts.Threshold; end
    if isfield(opts,'wp')
        recflag=0;
        for vn=vst:nargin-2
            if strcmpi(varargin{vn},'Window'), recflag=1; end
            if strcmpi(varargin{vn},'f0'), recflag=1; end
            if strcmpi(varargin{vn},'RelTol'), recflag=1; end
        end
        wp=opts.wp;
    end
    if isfield(opts,'fstepsim') && recflag==1, fstep=opts.fstepsim; end
end
for vn=vst:2:nargin-2
    if strcmpi(varargin{vn},'Window'), Window=varargin{vn+1};
    elseif strcmpi(varargin{vn},'f0'), f0=varargin{vn+1};
    elseif strcmpi(varargin{vn},'fmin'), fmin=varargin{vn+1};
    elseif strcmpi(varargin{vn},'fmax'), fmax=varargin{vn+1};
    elseif strcmpi(varargin{vn},'fstep'), fstep=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Padding'), PadMode=varargin{vn+1};
    elseif strcmpi(varargin{vn},'RelTol'), RelTol=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Preprocess'), Preprocess=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Plot'), PlotMode=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Display'), DispMode=varargin{vn+1};
    elseif strcmpi(varargin{vn},'CutEdges'), CutEdges=varargin{vn+1};
    elseif strcmpi(varargin{vn},'CutFreq'), CutFreq=varargin{vn+1};
    elseif strcmpi(varargin{vn},'Threshold'), Threshold=varargin{vn+1};
    else
        error(['There is no Property "',varargin{vn},'" (which is ',num2str(1+(vn-vst)/2,'%d'),...
            '''th out of ',num2str(ceil((nargin-1-vst)/2),'%d'),' specified)']);
    end
end

%===========================Window function================================
%[fwt] and [twf] - window function in frequency and time
%wp - structure with window parameters containing fields:
%     xi1,xi2 - window full support in the frequency domain;
%     ompeak,tpeak - window peak frequency and peak time
%     t1,t2 - window full support in the time domain;
%     fwtmax,twfmax - maximum value of abs(fwt) and abs(twf);
%     C,omg - coefficients needed for reconstruction
%--------------------------------------------------------------------------
fwt=[]; twf=[];
if recflag==1
    wp=struct;
    wp.fwtmax=[]; wp.twfmax=[]; wp.C=[]; wp.omg=[];
    wp.xi1=-Inf; wp.xi2=Inf; wp.ompeak=[];
    wp.t1=-Inf; wp.t2=Inf; wp.tpeak=[];
end
if iscell(Window)
    fwt=Window{1};
    if length(Window{2})==2, wp.xi1=Window{2}(1); wp.xi2=Window{2}(2); end
    twf=Window{3};
    if length(Window{4})==2, wp.t1=Window{4}(1); wp.t2=Window{4}(2); end
elseif strcmpi(Window,'Gaussian')
    fwt=@(xi)exp(-(f0^2/2)*xi.^2);
    twf=@(t)(1/sqrt(2*pi)/f0)*exp(-t.^2/(2*f0^2));
    wp.ompeak=0; wp.C=pi*twf(0); wp.omg=0;
    wp.tpeak=0;
elseif strcmpi(Window,'Exp')
    q=6.5*f0;
    twf=@(t)exp(-abs(t)/q);
    fwt=@(xi)2*q./(1+(q^2)*xi.^2);
    wp.ompeak=0; wp.C=pi*twf(0); wp.omg=0; wp.tpeak=0;
elseif strcmpi(Window,'Hann')
    q=4.4*f0;
    twf=@(t)(1+cos(2*pi*t/q))/2; wp.t1=-q/2; wp.t2=q/2;
    fwt=@(xi)(-(2*pi/q)^2)*sin(xi*q/2)./(xi.*(xi.^2-(2*pi/q).^2));
    wp.ompeak=0; wp.C=pi*twf(0); wp.omg=0; wp.tpeak=0;
elseif strcmpi(Window,'Blackman')
    q=5.6*f0; alpha=0.16;
    twf=@(t)(1+cos(2*pi*t/q))/2-alpha*(1+cos(4*pi*t/q))/2; wp.t1=-q/2; wp.t2=q/2;
    fwt=@(xi)((-(2*pi/q)^2)*sin(xi*q/2)./xi).*(1./(xi.^2-(2*pi/q).^2)-4*alpha./(xi.^2-(4*pi/q).^2));
    wp.ompeak=0; wp.C=pi*twf(0); wp.omg=0; wp.tpeak=0;
elseif ~isempty(strfind(lower(Window),'kaiser'))
    a=3; if length(Window)>6, a=str2double(Window(8:length(Window))); end
    q=3*sqqrt(1+abs(a-1/a))*f0;
    B=besseli(0,pi*a);
    twf=@(t)besseli(0,pi*a*sqrt(1-(2*t/q).^2))/B; wp.t1=-q/2; wp.t2=q/2;
    wp.C=pi*twf(0); wp.omg=0; wp.tpeak=0;
else
    error('Invalid window name');
end
%==========================================================================


%Determine parameters of the window function (and frequency bin width if needed)
if recflag==1
    if strcmpi(DispMode,'on'), fprintf('Estimating window parameters...\n'); end
    parcalc(RelTol); %calculate/update window parameters
end
coib1=ceil(abs(wp.t1e*fs)); coib2=ceil(abs(wp.t2e*fs)); %cone-of-influence edges


%Produce warning if there are no (S)WFT in the cone-of-influence
if wp.t2e-wp.t1e>L/fs
    if ~strcmpi(DispMode,'off')
        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
        wstr='';
        wstr=[wstr,sprintf('For used window function parameters and signal time-length there are no (S)WFT coefficients determined with specified accuracy %1.2e.\n',RelTol)];
        wstr=[wstr,sprintf('Transform is inaccurate under the specified precision. Either allow lower accuracy by increasing ''RelTol'', or consider changing window resolution parameter(s).\n')];
        if strcmpi(CutEdges,'on')
            wstr=[wstr,sprintf('Continuing without cone-of-influence.\n')];
        end
        fprintf(wstr);
        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
    end
    if strcmpi(CutEdges,'on'), CutEdges='off'; end
end


%Define WFT frequencies
wfmin=fmin+wp.xi1e/2/pi; wfmax=fmax+wp.xi2e/2/pi; %minimal and maximal frequency for WFT so that SWFT satisfy accuracy RelTol in range of [fmin,fmax]
fstepsim=fstep; wp.fstep=fstep;
if ~isempty(strfind(fstep,'auto')) %determine frequency step [fstep] if needed
    Nb=10; if length(fstep)>4, Nb=str2double(fstep(6:length(fstep))); end
    wp.fstep=(wp.xi2h-wp.xi1h)/(2*pi*Nb);
    c10=floor(log10(wp.fstep)); fdig=floor(wp.fstep/(10^c10));
    fstep=floor(wp.fstep/(10^c10))*(10^c10);
    if strcmpi(DispMode,'on')
        fprintf('Optimal frequency bin width was determined to be %f Hz (rounded to %d x 10^%d)\n',wp.fstep,fdig,c10);
    end
end
wfreq=(floor(wfmin/fstep)*fstep:fstep:ceil(wfmax/fstep)*fstep)'; % frequencies (determined solely by fstep, which by default depends only on window parameters)
SN=length(wfreq); % number of frequencies


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
        w=2.^(-(L/fs-(1:L)/fs)/(wp.t2h-wp.t1h));
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
        w=2.^(-(L/fs-(1:L)/fs)/(wp.t2h-wp.t1h));
        padright=fcast(signal,fs,n2,[max([fmin,fs/L]),fmax],min([ceil(SN/2)+5,round(L/3)]),w);
    elseif isnumeric(p2), padright=p2*ones(n2,1);
    elseif strcmpi(p2,'symmetric'), padright=[flipud(signal(L-min([n2,L])+1:L));zeros(n2-L,1)];
    elseif strcmpi(p2,'periodic'), padright=[signal(1:min([n2,L]));zeros(n2-L,1)];
    else error('Bad ''Padding'' property');
    end
    signal=[padleft;signal;padright];
elseif strcmpi(PadMode,'predictive')
    if strcmpi(DispMode,'on'), fprintf('Applying predictive padding (may take some time)...\n'); end
    w=2.^(-(L/fs-(1:L)/fs)/(wp.t2h-wp.t1h));
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


%Windowed Fourier Transform by itself
WFT=zeros(SN,L)*NaN; IFR=zeros(SN,L)*NaN; ouflag=0; if wp.t2e-wp.t1e>L/fs, coib1=0; coib2=0; end
if strcmpi(DispMode,'on'), pos=0; fprintf('Calculating Windowed Fourier Transform and its phase growth rate (%d frequencies from %3.3f to %3.3f): ',SN,wfreq(1),wfreq(end)); end
for sn=1:SN
    freqwf=wfreq(sn)-ff; %frequencies for the window function
    ii=find(freqwf>wp.xi1/2/pi & freqwf<wp.xi2/2/pi); %take into account only frequencies within the window support
    if ~isempty(fwt)
        fw=fwt(2*pi*freqwf(ii)); nid=find(isnan(fw) | ~isfinite(fw));
        if ~isempty(nid) %to avoid NaNs due to numerics, e.g. sin(0)/0
            fw(nid)=fwt(2*pi*freqwf(ii(nid))+10^(-14));
            nid=find(isnan(fw) | ~isfinite(fw)); fw(nid)=0;
            if ~isempty(nid), ouflag=1; ouval=2*pi*freqwf(nid(1)); end
        end
    else
        timewf=(1/fs)*[-(1:ceil((NL-1)/2))+1,NL+1-(ceil((NL-1)/2)+1:NL)]';
        jj=find(timewf>wp.t1 & timewf<wp.t2); tw=zeros(NL,1); %take into account only times within the window support
        tw(jj)=twf(timewf(jj)).*exp(-1i*2*pi*wfreq(sn)*timewf(jj)); nid=find(isnan(tw) | ~isfinite(tw));
        if ~isempty(nid) %to avoid NaNs due to numerics, e.g. sin(0)/0
            tw(nid)=twf(timewf(nid)+10^(-14));
            nid=find(isnan(tw) | ~isfinite(tw)); tw(nid)=0;
            if ~isempty(nid), ouflag=1; ouval=timewf(nid(1)); end
        end
        fw=(1/fs)*fft(tw); fw=fw(ii);
    end
    cc=zeros(NL,1); cc(ii)=fx(ii).*fw(:); %convolution in the frequency domain
    out=ifft(cc,NL); % calculate WFT at each time
    WFT(sn,1:L)=out(1+n1:NL-n2);
    
    dcc=zeros(NL,1); dcc(ii)=fx(ii).*(1i*2*pi*ff(ii)).*fw(:); %convolution in the frequency domain
    dout=ifft(dcc,NL); % calculate WFT time derivative at each time
    IFR(sn,1:L)=(1/2/pi)*imag(dout(1+n1:NL-n2)./out(1+n1:NL-n2));
    
    if strcmpi(DispMode,'on') && floor(100*sn/SN)>floor(100*(sn-1)/SN)
        cstr=num2str(floor(100*sn/SN)); fprintf([repmat('\b',1,pos),cstr,'%%']); pos=length(cstr)+1;
    end
end
if strcmpi(DispMode,'on'), fprintf('\n'); end
if ouflag==1
    if ~isempty(fwt)
        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
        fprintf('Possibly overflow/underflow (e.g. Inf/Inf): specified frequency-domain form of window function\n');
        fprintf('returns NaN or Inf, e.g. when its argument is %e. In all such cases it is set to zero.\n',ouval);
        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
    else
        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
        fprintf('Possibly overflow/underflow (e.g. Inf/Inf): specified time-domain form of window function\n');
        fprintf('returns NaN or Inf, e.g. when its argument is %e. In all such cases it is set to zero.\n',ouval);
        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
    end
end

%Set to NaN all WFT coefficients outside the cone of influence if specified
%(the SWFT will be then cut to the cone of influence automatically)
if strcmpi(CutEdges,'on')
    WFT(:,1:coib1)=NaN; WFT(:,1+L-coib2:L)=NaN;
    IFR(:,1:coib1)=NaN; IFR(:,1+L-coib2:L)=NaN;
end
if nargout>5
    if strcmpi(CutFreq,'on'), varargout{4}=IFR(wfreq>=fmin & wfreq<=fmax,:);
    else varargout{4}=IFR; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Synchrosqueezing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate borders (boms) and widths (doms) of frequency bins
freq=wfreq(wfreq>=fmin & wfreq<=fmax); %frequencies for SWFT (same as for WFT in the corresponding range)
FN=length(freq); % number of frequencies

%Filter out white noise (if needed) using hard-thresholding
if strcmpi(Threshold,'MAD')
    iwfreq=find(wfreq>(fmax+fmin)/2);
    MAD=zeros(length(iwfreq),1); %median absolute deviation
    for sn=1:length(iwfreq)
        cwn=iwfreq(sn);
        MAD(sn)=median(abs(abs(WFT(cwn,~isnan(WFT(cwn,:))))-median(abs(WFT(cwn,~isnan(WFT(cwn,:)))))));
    end
    cthr=1.4826*sqrt(2*log(L))*median(MAD);
else
    cthr=Threshold;
end
if cthr>0
    WFT(abs(WFT)<cthr)=0; %set WFT coefficients lower than threshold to zero
end

%Transform phase growth rate [IFR] to numbers of frequency bins
IFR=1+floor((1/2)+(IFR-freq(1))/fstep); %for equally-spaced frequencies
IFR(IFR<1 | IFR>FN | isnan(IFR))=0;

%Very fast implementation of synchrosqueezing
if strcmpi(DispMode,'on'), fprintf('Synchrosqueezing...\n'); end
[allsc,alltt,allfr]=find(IFR); clear IFR;
allsv=(2*pi*fstep/wp.C)*WFT(sub2ind([SN,L],allsc,alltt));
SWFT=sparse(allfr,alltt,allsv,FN,L); clear allsc allfr alltt allsv;

%Cut or delete WFT to free more memory
if nargout<4, clear WFT;
elseif strcmpi(CutFreq,'on')
    WFT=WFT(wfreq>=fmin & wfreq<=fmax,:);
    wfreq=freq;
end

SWFT=full(SWFT); %transform SWFT to a full form

if ~strcmpi(PlotMode,'off')
    if strcmpi(DispMode,'on'), fprintf('Plotting...\n'); end
    scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,scrsz(3)/2,6*scrsz(4)/8]);
    axes('Position',[0.15,0.1,0.8,0.5333],'Layer','top','Box','on','FontSize',16);
    hold all;
    
    YY=freq; XX=(0:(L-1))/fs; ZZ=abs(SWFT); ZZname='SWFT amplitude';
    if ~isempty(strfind(lower(PlotMode),'pow')), ZZ=ZZ.^2; ZZname='SWFT power'; end
    MYL=round(scrsz(3)); MXL=round(scrsz(4)); %maximum number of points seen in plots
    if isempty(strfind(lower(PlotMode),'wr')) && (size(ZZ,1)>MYL || size(ZZ,2)>MXL)
        if ~strcmpi(DispMode,'off')
            fprintf('Warning: SWFT contains more data points (%d x %d) than pixels in the plot, so for a\n',size(ZZ,1),size(ZZ,2));
            fprintf('         better performance its resampled version (%d x %d) will be displayed instead.\n',min([MYL,size(ZZ,1)]),min([MXL,size(ZZ,2)]));
        end
        if size(ZZ,1)>MYL, YY=linspace(freq(1),freq(end),MYL); end
        if size(ZZ,2)>MXL, XX=linspace(0,(L-1)/fs,MXL); end
        ZZ=aminterp((0:(L-1))/fs,freq,ZZ,XX,YY,'max'); XX=XX(:); YY=YY(:);
    end
    TL=length(XX); FL=length(YY);
    pc=pcolor(YY,XX,ZZ'); title(ZZname);
    
    xlabel('Frequency (Hz)'); ylabel('Time (s)');
    ylim([0,(L-1)/fs]); xlim([freq(1),freq(end)]); set(pc,'EdgeColor','none');
    
    coib1(coib1==0)=NaN; coib2(coib2==0)=NaN;
    plot([freq(1);freq(end)],[coib1;coib1]/fs,'-k','LineWidth',2); plot([freq(1);freq(end)],[L-coib2+1;L-coib2+1]/fs,'-k','LineWidth',2);
    coib1(isnan(coib1))=0; coib2(isnan(coib2))=0; if strcmpi(CutEdges,'off'), coib1(:)=0; coib2(:)=0; end
    
    if ~isempty(strfind(lower(PlotMode),'+'))
        axes('Position',[0.15,0.7,0.8,0.25],'XLim',[freq(1),freq(end)],'XTickLabel',{},'FontSize',16);
        hold all;
        ni1=find(~isnan(ZZ(1,:)),1,'first'); ni2=find(~isnan(ZZ(1,:)),1,'last'); mx=mean(ZZ(:,ni1:ni2),2);
        plot(YY,mx,'-k','LineWidth',2); ylabel({'Time-averaged',ZZname});
        if max(mx)>0, ylim([0,1.1*max(mx)]); end
        
        if ~isempty(strfind(lower(PlotMode),'++'))
            sZZ=sort(ZZ(:,ni1:ni2),2); ZL=size(sZZ,2);
            lx=sZZ(:,max([1,round(0.025*ZL)])); ux=sZZ(:,round(0.975*ZL));
            idnn=find(~isnan(lx)); % not-NaN indices
            fill([YY(idnn);flipud(YY(idnn))],[lx(idnn);flipud(ux(idnn))],[1,1,0],'FaceAlpha',0.5);
            if max(ux)>0, ylim([0,1.1*max(ux)]); end
        end
    end
end

if nargout>3, varargout{2}=WFT; end
if nargout>4, varargout{3}=wfreq; end

if nargout>2
    wopt=struct; %simulation parameters
    wopt.wp=wp; %parameters of the window
    wopt.TFRname='SWFT'; wopt.fs=fs;
    wopt.Window=Window;
    wopt.f0=f0;
    wopt.fmin=fmin;
    wopt.fmax=fmax;
    wopt.fstep=fstep; wopt.fstepsim=fstepsim;
    wopt.Padding=PadMode;
    wopt.RelTol=RelTol;
    wopt.Preprocess=Preprocess;
    wopt.Plot=PlotMode;
    wopt.Display=DispMode;
    wopt.CutEdges=CutEdges;
    wopt.Threshold=Threshold;
    wopt.CutFreq=CutFreq;
    
    varargout{1}=wopt;
end







%----------------------------------------------------------------------------------------------------------------------
%======================================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================================================================
%----------------------------------------------------------------------------------------------------------------------
    % Based on window function in time [twf] and frequency [fwt],
    % determines window parameters such as constant [Cg],
    % peaks in time [tpeak] and frequency [ompeak],
    % epsilon-support in frequency [xi1e,xi2e] and in time [t1e,t2e],
    % 50%-support in frequency [xi1h,xi2h] and time [t1h,t2h],
    % time-frequency resolution [tfres] (inverse of multiplication of the latter),
    % and frequency step [fstep] (if 'auto') for specified relative accuracy [racc] (=epsilon);
    % assigns all these values into the window parameters structure [wp].
    function parcalc(racc)
        racc=min(racc,0.5-10^(-10)); %current \epsilon
        ctol=max([racc/1000,10^(-12)]); %parameter of numerical accuracy
        MIC=max([100000,10*L]); %maximum interval count for one-time calculations
        
        %Times and frequencies for testing (with 4 times higher time duration and sampling frequency than the signal)
        nt=(1/(4*fs))*(-8*L+1:8*L-1)'; nt=nt(nt>wp.t1 & nt<wp.t2);
        nxi=(2*pi*(4*fs)/(16*L-1))*(-8*L+1:8*L-1)'; nxi=nxi(nxi>wp.xi1 & nxi<wp.xi2);
        
        %==================================================================
        %Determine values for known frequency and/or time-domain forms
        if ~isempty(fwt) %if the frequency-domain form is known
            wp.fwt=fwt;
            if isempty(wp.ompeak) %peak frequency
                ipeak=find(abs(fwt(nxi))==max(abs(fwt(nxi)))); wp.ompeak=mean(nxi(ipeak));
                wp.ompeak=fminsearch(@(x)-abs(fwt(x)),wp.ompeak,optimset('TolX',10^(-14),'Display','off'));
            end
            if isempty(wp.fwtmax)
                wp.fwtmax=fwt(wp.ompeak);
                if isnan(wp.fwtmax), wp.fwtmax=fwt(wp.ompeak+10^(-14)); end
            end
            if abs(wp.ompeak)>10^(-12) %center around the peak if needed
                if ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                    fprintf('The Fourier transform (FT) of the specified window function is not centered around its maximum, which occurs at %f.\n',wp.ompeak)
                    fprintf('Therefore, we consider its modified version shifted in frequency so that the window FT amplitude is now peaked at zero.\n');
                    fprintf(2,'----------------------------------------------------------------------------------------------------\n');
                end
                fwt=@(xi)fwt(xi+wp.ompeak);
                if ~isempty(twf)
                    twf=@(t)twf(t).*exp(-1i*wp.ompeak*t);
                end
                wp.xi1=wp.xi1-wp.ompeak; wp.xi2=wp.xi2-wp.ompeak;
                wp.fwt=fwt; wp.ompeak=0;
            end
            vfun=@(u)fwt(u); xp=wp.ompeak; lim1=wp.xi1; lim2=wp.xi2;
            
            [QQ,wflag,xx,ss]=sqeps(vfun,xp,[lim1,lim2],racc,MIC,...
                [-8*(2*pi*fs),8*(2*pi*fs)]); %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            wp.xi1e=ss(1,1); wp.xi2e=ss(1,2); wp.xi1h=ss(2,1); wp.xi2h=ss(2,2);
            if isempty(wp.C)
                if exist('twf','var') && ~isempty(twf)
                    wp.C=pi*twf(0);
                    if isnan(wp.C), wp.C=pi*twf(10^(-14)); end
                else
                    wp.C=(QQ(1,1)+QQ(1,2))/2;
                end
            end
            if wflag==1 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                fprintf('The frequency-domain window function is not well-behaved (e.g. decays very slowly as frequency tends\n');
                fprintf('to +/- infinity). The integration might be not accurate (and therefore e.g. the calculated frequency \n');
                fprintf('step ''fstep'', if set to ''auto'', so frequency discretization might be also not appropriate).\n');
                fprintf(2,'----------------------------------------------------------------------------------------------------\n');
            end
            
            if isempty(wp.omg) %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                wstate=warning('off','all');
                px1=min([wp.ompeak-xx(1,1),xx(1,2)-wp.ompeak]); px2=min([wp.ompeak-xx(4,1),xx(4,2)-wp.ompeak]);
                [Y1,errY1]=quadgk(@(u)(u.*fwt(wp.ompeak+u)-u.*fwt(wp.ompeak-u)),0,px1,'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                [Y2,errY2]=quadgk(@(u)(u.*fwt(wp.ompeak+u)-u.*fwt(wp.ompeak-u)),px1,px2,'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                [Y3,errY3]=quadgk(@(u)-u.*fwt(wp.ompeak-u),px2,wp.ompeak-xx(4,1),'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                [Y4,errY4]=quadgk(@(u)u.*fwt(wp.ompeak+u),px2,xx(4,2)-wp.ompeak,'MaxIntervalCount',2*MIC,'AbsTol',0,'RelTol',10^(-12));
                if abs((errY1+errY2+errY3+errY4)/(Y1+Y2+Y3+Y4))<10^(-4), wp.omg=wp.ompeak+(Y1+Y2+Y3+Y4)/(2*wp.C); else wp.omg=Inf; end
                warning(wstate);
            end
            
            
            if isempty(twf) %if time domain form is not known
                [PP,wflag,xx,ss]=sqeps(@(x)abs(fwt(x)).^2,wp.ompeak,[wp.xi1,wp.xi2],racc,MIC,...
                    [-8*(2*pi*fs),8*(2*pi*fs)]); %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                Etot=sum(PP(1,:))/2/pi;
                
                CL=2^nextpow2(MIC/8); CT=CL/(2*abs(ss(1,2)-ss(1,1)));
                CNq=ceil((CL+1)/2); cxi=(2*pi/CT)*(CNq-CL:CNq-1)'; idm=find(cxi<=wp.xi1); idc=find(cxi>wp.xi1 & cxi<wp.xi2); idp=find(cxi>=wp.xi2);
                Cfwt=[zeros(length(idm),1);fwt(cxi(idc));zeros(length(idp),1)]; idnan=find(isnan(Cfwt));
                if ~isempty(idnan), idnorm=find(~isnan(Cfwt)); Cfwt(idnan)=interp1(idnorm,Cfwt(idnorm),idnan,'spline','extrap'); end
                Ctwf=ifft((CL/CT)*Cfwt([CL-CNq+1:CL,1:CL-CNq])); Ctwf=Ctwf([CNq+1:CL,1:CNq]);
                Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                Eest=(CT/CL)*sum(Etwf); perr=Inf;
                while (abs(Etot-Eest)+Iest1+Iest2)/Etot<=perr
                    CT=CT/2; perr=(abs(Etot-Eest)+Iest1+Iest2)/Etot;
                    CNq=ceil((CL+1)/2); cxi=(2*pi/CT)*(CNq-CL:CNq-1)'; idm=find(cxi<=wp.xi1); idc=find(cxi>wp.xi1 & cxi<wp.xi2); idp=find(cxi>=wp.xi2);
                    Cfwt=[zeros(length(idm),1);fwt(cxi(idc));zeros(length(idp),1)]; idnan=find(isnan(Cfwt));
                    if ~isempty(idnan), idnorm=find(~isnan(Cfwt)); Cfwt(idnan)=interp1(idnorm,Cfwt(idnorm),idnan,'spline','extrap'); end
                    Ctwf=ifft((CL/CT)*Cfwt([CL-CNq+1:CL,1:CL-CNq])); Ctwf=Ctwf([CNq+1:CL,1:CNq]);
                    Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                    Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                    Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                    Eest=(CT/CL)*sum(Etwf);
                end
                CL=16*CL; CT=CT*2;
                CNq=ceil((CL+1)/2); cxi=(2*pi/CT)*(CNq-CL:CNq-1)'; idm=find(cxi<=wp.xi1); idc=find(cxi>wp.xi1 & cxi<wp.xi2); idp=find(cxi>=wp.xi2);
                Cfwt=[zeros(length(idm),1);fwt(cxi(idc));zeros(length(idp),1)]; idnan=find(isnan(Cfwt));
                if ~isempty(idnan), idnorm=find(~isnan(Cfwt)); Cfwt(idnan)=interp1(idnorm,Cfwt(idnorm),idnan,'spline','extrap'); end
                Ctwf=ifft((CL/CT)*Cfwt([CL-CNq+1:CL,1:CL-CNq])); Ctwf=Ctwf([CNq+1:CL,1:CNq]);
                Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                Eest=(CT/CL)*sum(Etwf);
                
                if (abs(Etot-Eest)+Iest1+Iest2)/Etot>0.01 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                    fprintf(['Cannot accurately invert the specified frequency-domain form of the window function to find its\n',...
                        'time domain form and corresponding characteristics (e.g. cone-of-influence borders).\n',...
                        'This might be because the window function decays too slowly in time or frequency.\n']);
                    fprintf(2,'----------------------------------------------------------------------------------------------------\n');
                end
                
                
                Ctwf=Ctwf(1:2*CNq-3); ct=(CT/CL)*(-(CNq-2):CNq-2)'; %make symmetric
                wp.twf={Ctwf,ct};
                
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
                ipeak=find(abs(twf(nt))==max(abs(twf(nt)))); wp.tpeak=mean(nt(ipeak));
                wp.tpeak=fminsearch(@(x)-abs(twf(x)),wp.tpeak,optimset('TolX',10^(-14),'Display','off'));
            end
            
            %Calculate the frequency domain characteristics first, if not
            %known (will be needed afterwards, mainly [wp.ompeak])
            if isempty(fwt) %if frequency domain form is not known
                [PP,wflag,xx,ss]=sqeps(@(x)abs(twf(x)).^2,wp.tpeak,[wp.t1,wp.t2],racc,MIC,...
                    [-8*L/fs,8*L/fs]); %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
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
                Cfwt=(CT/CL)*fft(Ctwf([CL-CNq+1:CL,1:CL-CNq])); Cfwt=Cfwt([CNq+1:CL,1:CNq]);
                Etwf=abs(Ctwf).^2; Efwt=abs(Cfwt).^2;
                Iest1=(CT/CL)*sum(abs(Etwf(3:end)-2*Etwf(2:end-1)+Etwf(1:end-2)))/24; %error of integration in time
                Iest2=(1/CT)*sum(abs(Efwt(3:end)-2*Efwt(2:end-1)+Efwt(1:end-2)))/24; %error of integration in frequency
                Eest=(1/CT)*sum(Efwt);
                
                if (abs(Etot-Eest)+Iest1+Iest2)/Etot>0.01 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                    fprintf(['Cannot accurately invert the specified time-domain form of the window function to find its\n',...
                        'frequency-domain form and corresponding characteristics (e.g. optimal frequency step ''fstep'').\n',...
                        'This might be because the window function decays too slowly in time or frequency.\n']);
                    fprintf(2,'----------------------------------------------------------------------------------------------------\n');
                end
                
                Cfwt=Cfwt(1:2*CNq-3); cxi=(2*pi/CT)*(-(CNq-2):CNq-2)'; %make symmetric
                wp.fwt={Cfwt,cxi};
                
                %Estimate general parameters
                if isempty(wp.ompeak) %peak frequency
                    ipeak=find(abs(Cfwt)==max(abs(Cfwt)));
                    if length(ipeak)==1
                        a1=abs(Cfwt(ipeak-1)); a2=abs(Cfwt(ipeak)); a3=abs(Cfwt(ipeak+1));
                        wp.ompeak=cxi(ipeak);
                        if abs(a1-2*a2+a3)>2*eps, %use queadratic interpolation to find exact peak location
                            wp.ompeak=wp.ompeak+(1/2)*(a1-a3)/(a1-2*a2+a3)*(2*pi/CT);
                        end
                    else
                        wp.ompeak=mean(cxi(ipeak));
                    end
                end
                if isempty(wp.fwtmax)
                    [~,ipeak]=min(abs(cxi-wp.ompeak));
                    wp.fwtmax=interp1(cxi(ipeak-1:ipeak+1),abs(Cfwt(ipeak-1:ipeak+1)),wp.ompeak,'spline');
                end
                if isempty(wp.C) %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    wp.C=pi*twf(0);
                    if isnan(wp.C), wp.C=pi*twf(10^(-14)); end
                end
                if isempty(wp.omg) %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                    wp.omg=sum((2*pi/CT)*cxi.*Cfwt)/(2*wp.C);
                end
                
                %Calculate the cumulative integrals
                cxi=[cxi-pi/CT;cxi(end)+pi/CT]; %to use midpoint rule
                CS=(2*pi/CT)*cumsum(Cfwt); CS=[0;CS(:)]/CS(end); CS=abs(CS);
                ICS=(2*pi/CT)*cumsum(Cfwt(end:-1:1)); ICS=ICS(end:-1:1); ICS=[ICS(:);0]/ICS(1); ICS=abs(ICS);
                
                %Estimate epsilon-supports
                xid=find(CS(1:end-1)<racc/2 & CS(2:end)>=racc/2,1,'first');
                if isempty(xid), wp.xi1e=cxi(1);
                else
                    a1=CS(xid)-racc/2; a2=CS(xid+1)-racc/2;
                    wp.xi1e=cxi(xid)-a1*(cxi(xid+1)-cxi(xid))/(a2-a1);
                end
                xid=find(ICS(1:end-1)>=racc/2 & ICS(2:end)<racc/2,1,'last');
                if isempty(xid), wp.xi2e=cxi(end);
                else
                    a1=ICS(xid)-racc/2; a2=ICS(xid+1)-racc/2;
                    wp.xi2e=cxi(xid)-a1*(cxi(xid+1)-cxi(xid))/(a2-a1);
                end
                xid=find(CS(1:end-1)<0.25 & CS(2:end)>=0.25,1,'first');
                if isempty(xid), wp.xi1h=cxi(1);
                else
                    a1=CS(xid)-0.25; a2=CS(xid+1)-0.25;
                    wp.xi1h=cxi(xid)-a1*(cxi(xid+1)-cxi(xid))/(a2-a1);
                end
                xid=find(ICS(1:end-1)>=0.25 & ICS(2:end)<0.25,1,'last');
                if isempty(xid), wp.xi2h=cxi(end);
                else
                    a1=ICS(xid)-0.25; a2=ICS(xid+1)-0.25;
                    wp.xi2h=cxi(xid)-a1*(cxi(xid+1)-cxi(xid))/(a2-a1);
                end
                
                
                if abs(wp.ompeak)>10^(-12) %center around the peak
                    if ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                        fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                        fprintf('The Fourier transform (FT) of the specified window function is not centered around its maximum, which occurs at %f.\n',wp.ompeak)
                        fprintf('Therefore, we consider its modified version shifted in frequency so that the window FT amplitude is now peaked at zero.\n');
                        fprintf(2,'----------------------------------------------------------------------------------------------------\n');
                    end
                    wp.xi1=wp.xi1-wp.ompeak; wp.xi2=wp.xi2-wp.ompeak;
                    wp.xi1e=wp.xi1e-wp.ompeak; wp.xi2e=wp.xi2e-wp.ompeak;
                    wp.xi1h=wp.xi1h-wp.ompeak; wp.xi2h=wp.xi2h-wp.ompeak;
                    wp.omg=wp.omg-wp.ompeak;
                    
                    twf=@(t)twf(t).*exp(-1i*wp.ompeak*t); wp.twf=twf;
                    Ctwf=[zeros(length(idm),1);twf(ct(idc));zeros(length(idp),1)]; idnan=find(isnan(Ctwf));
                    if ~isempty(idnan), idnorm=find(~isnan(Ctwf)); Ctwf(idnan)=interp1(idnorm,Ctwf(idnorm),idnan,'spline','extrap'); end
                    Cfwt=(CT/CL)*fft(Ctwf([CL-CNq+1:CL,1:CL-CNq])); Cfwt=Cfwt([CNq+1:CL,1:CNq]);
                    Cfwt=Cfwt(1:2*CNq-3); cxi=(2*pi/CT)*(-(CNq-2):CNq-2)'; %make symmetric
                    
                    wp.fwt={Cfwt,cxi}; wp.ompeak=0;
                end
                
            end
            
            if isempty(wp.twfmax)
                wp.twfmax=twf(wp.tpeak);
                if isnan(wp.twfmax), wp.twfmax=twf(wp.tpeak+10^(-14)); end
            end
            vfun=@(u)twf(u); xp=wp.tpeak; lim1=wp.t1; lim2=wp.t2;
            
            [QQ,wflag,xx,ss]=sqeps(vfun,xp,[lim1,lim2],racc,MIC,...
                [-8*L/fs,8*L/fs]); %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            wp.t1e=ss(1,1); wp.t2e=ss(1,2); wp.t1h=ss(2,1); wp.t2h=ss(2,2);
            if wflag==1 && ~strcmpi(DispMode,'off') %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                fprintf(2,'--------------------------------------------- Warning! ---------------------------------------------\n');
                fprintf('The time-domain window function is not well-behaved (e.g. decays very slowly as time tends to +/- infinity).\n');
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
