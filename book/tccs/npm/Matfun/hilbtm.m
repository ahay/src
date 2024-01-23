%function hret = hilbtm(IMF,chkplot)
%
%   INPUT:
%         IMF: an IMF 2D matrix,with dimension(Npt,Nimf)
%     chkplot: when chkplot=1,plot result on the screen,default=0
%  OUTPUT:
%        hret: the improved hilbert transform result
%
%  NOTE: 
%     1.This code is an improved version of Hilbert Transform of matlab.
%       The Matlab version uses '2 times fft method' for Hilbert-Transform calculation.
%     2.For FFT computation, a signal waveform been repeated and repeated from (-oo ~ +oo).
%       So if the wave have different value ,slope and curvature  at beginning and endding point,
%       that makes a jump appearance after the wave is repeated again and again.This discontinuous 
%       makes the FFT result a problem,and is called 'Gibbs phenomenon'.
%     3.This hilbtm.m, it uses perfect cos waves instead of the previous hilbt_m used. 
%       Also,instead of just replacing an end section of the original wave with the new wave it uses
%       a weight to gradually introduce the new wave. This results in more subler jumps of frequency 
%       at the joint from the original. Data appears to be "better" in this version after the hilbert
%       transform - frequency is smoother, less jitter from the Gibbs effect.
%
% References:
%   
%
%
% code writer: Karin Blank (karin.blank@nasa.gov)
%              modified to fit into fah and deal with multiple IMFs KBB
% code writer: S.C.Su-add a figure for checking result
% footnote:S.C.Su 2009/10/13
%
% This code is hard to express the structure clearly.Brief it as follows:
%  1.have 2 parts--fixing the head,fixing the end
%  2.maybe 2 conditions---maximun near the end,minimum near the end
%  3.by the above 1.and 2.,there are 4 different scenario :
%    head-max,head-min,end-max,end-min
% So there are some sentences for finding extrema,zero-crossing ,
%    and make judgement for which extrema near the head or end.
%
% After 4 different scenario is decided,this code extend the ends out by cosine wave
%  4.the cosine have the same frequency for the first extrema
%  5.Adding the extra cosine from the first extrema,to a zero value,zero slope,zero curvature position
%    After adding those extra cosine wave in the ends,the result waveform start and end with  
%    the same 'zero value,zero slope,zero curvature 'position at beginning and endding points. 
% So there are some sentences for finding frequency,extending adequate distance,arranging the cosine connected position   
%
% Because the Gibbs phenomenon is fixed,we use matlab function 'hilbert.m' to process Hilbert Transformation
%
% Association: hilbert.m
% this function is a improved edition used for calculating the Hilbert-Transform. 
% matlab original function disadvantage is been solved.
%
% Concerned function: hilbert.m
%                     above mentioned m file must be put together
%
%

function hret = hilbtm(IMF,chkplot)

%when you need to check Hilbtm.m result,let chkplot=1
%else chkplot=0
if nargin<2
  chkplot = 0;
end
  
%check input for its dimension arrangement
orient = 0;
if(size(IMF, 1) > size(IMF,2))
    IMF = IMF.';
    orient = 1; %indicates needs to be flipped on return
end

% initial the result with zero 
hret = zeros(size(IMF,1), size(IMF,2));

%add warning sentence
 disp('When maximum less then zero,the result should be checked')
 disp('When minimum bigger then zero,the result should be checked')
    
% start to do the calculation IMF by IMF------loop A
for k=1:size(IMF,1)

% assign this IMF data as imf matrix    
    imf = IMF(k,:);
%Checking the necessarity of Gibbs phenonemon,the result :
% when bx=1 means beginning point needs no rectification ,ex=1 means endding point needs no rectification
% when bx=0 means beginning point needs rectification ,ex=0 means endding point needs rectification
    [bx, ex] = skiphilbt_m(imf);
    
    if(k==size(IMF,1)) %test if last element is a residual
        [max_x, max_y] = local_max(1:length(imf), imf);
        
        if(length(max_x) <= 2)
            bx = 1;
            ex = 1;
        end
    end
        
    thestart = 0;
    wavestart = [];
    wavend = [];
    y = imf;

% when bx=0 means beginning point needs rectification----loop B start
% start of extending,when bx=0 means the data need extending extra segemnt
    if (bx == 0)
        flagmax = 0;
        flagmin = 0;

        zeroc =0;

        n = length(imf);
        n_mn=length(imf);
        n_mx=length(imf);

%finds first wave from the beginning segment- look for and extrema and a zero crossing
% finding starts here------------------------------------loop C
        for i=2:n-1
        
            %finds local maximum
            if ((imf(i) > imf(i-1)) & (flagmax == 0)) %might be an extrema
                flagmax = 1;
                Xmax = i;
            end
            if((imf(i) < imf(i+1)) & (flagmax == 1)) %not an extrema
                flagmax = 0;
                Xmax = 0; %%
            end
            if((imf(i) > imf(i+1)) & (flagmax == 1)) %is an extrema
                flagmax = 0;
                n_mx = Xmax;
            end
            %the result 
            % if maximum is found its position is assigned as n_mx
            % if maximum is not found the n_mx value remains Npt

            %finds local minimum
            if((imf(i) < imf(i-1)) & (flagmin == 0)) %might be an extrema
                flagmin = 1;
                Xmin = i;
            end
            if((imf(i) > imf(i+1)) & (flagmin == 1)) %not an extrema
                flagmin = 0;
                Xmin = 0;
            end
            if((imf(i) < imf(i+1)) & (flagmin == 1)) %is an extrema
                flagmin = 0;

                n_mn = Xmin;
            end
            %the result 
            % if minimum is found its position is assigned as n_mn
            % if minimum is not found the n_mx value remains Npt

            %finds zero crossing
            if (((imf(i) == 0) | (imf(i) > 0 & imf(i+1) < 0) | (imf(i) < 0 & imf(i+1) > 0)) & zeroc <= 0)
                zeroc = i;
            end
            %the result
            %break out of loop if conditions are met (since matlab doesn't support do/while loops
            %set those result to plot the figure.
            if(((zeroc > 0) & (n_mn < length(imf) | n_mx < length(imf))))
                
                break; %after zero-crossing pt is founded ,jump out loop C,stop finding from beginning segment
            end
        end     %-----------------------------------------the end of loop C
           
        % this is for the last one IMF ,the almost trend component
        if(((zeroc > 0) + (n_mn > 0) + (n_mx > 0)) < 2)
            h=[];
            if orient == 1
                hret = hret.';
            end
            return;
        end

        %make sure there is one valid extrema and one valid zero crossing
        if(zeroc > 0 && (n_mn < length(imf) || n_mx < length(imf)))
        %start to process the beginning segemnt,produce extending part
       
        
            % after finding from beginning segment,we have a maximun,a minimum,a zero-crossing
            % Now we need tofind which extreme is more closer to the beginning end.
            % the following 'if' sentence is to process for maximum or minumum conditions        
            
            %if minimum closer to border--use first minimum value 
            if(n_mn < n_mx)
                %     wavstart = cos(-pi*3.5:pi/(abs(zeroc-n_mn)*2):pi*3)*abs(imf(n_mn));
                %estimate where zero crossing actually occurs
                slope = (imf(zeroc+1) - imf(zeroc)); %(the slope of the zero-crossing segment,delta x = 1)
                b = imf(zeroc+1) - (slope * (zeroc+1));%extend slope to y axis,b is the distence on y and imf(i+1)
                x = -b/slope; %accurate zero-crossing x coordinate is obtained.
                
                %find out frequency interval
                interval = abs(0.5*pi/(x-n_mn));
                interval_a = interval;
                
                %find the length of a correct wave
                wavestart_n = ceil((pi * 2.5)/interval);%for 2.5pi how many points
                
                %extend wave length if the extended part segment not long enough
                %extend to the sarting point
                if(wavestart_n <= n_mn)
                    i=ceil((((n_mn+1)*interval)/pi-2.5)/2);
                    wavestart_n = ceil((pi * (2.5 + i * 2))/interval);
                end

                %make the cosine wave 
                %(within the same period in the beginning-extend for 2.5 pi:0.5pi~3.0pi)
                %(within the first minimum value of the imf as the cosine amplitude)
                x2 = 0:interval:((wavestart_n-1)*interval);%the x-axis value of new coordinate
                wavestart = cos(pi*.5+x2) * abs(imf(n_mn));%the y-axis value of cosine function(0.5pi~3.0pi)

                %calculate weights for interpolation
                weight = 1/n_mn;
                weight_intr = weight;
                
                curr = wavestart_n;
                
                %interpolate between actual wave and artificial wave
                %it means extend out gradually between original wave and artificial cosine wave
                for i=n_mn+1:-1:1
                    wavestart(curr) = (imf(i) * (1-weight)) + (wavestart(curr) * weight);
                    curr = curr-1;
                    weight = weight + weight_intr;%weighting is always increasing step by step
                end

                y = [wavestart, imf(n_mn+1:end)];%add to the original in the heading position
                thestart = n_mn;
                
            else %maximum is closer to border--use first maximum value 
            
                %estimate where zero crossing actually occurs
                slope = (imf(zeroc+1) - imf(zeroc)); %(delta x = 1)%(the slope of the zero-crossing segment,delta x = 1)
                b = imf(zeroc+1) - (slope * (zeroc+1));%extend slope to y axis,b is the distence on y and imf(i+1)
                x = -b/slope;%accurate zero-crossing x coordinate is obtained.
                
                %find out frequency interval 
                interval = abs(pi/((x-n_mx)*2));
                interval_a = interval;
                
                wavestart_n = ceil((pi * 1.5)/interval);
                
                %extend wave length if the extended part segment not long enough
                %extend to the sarting point
                if(wavestart_n <= n_mx)
                   %calculate to the beginning point
                    i = ceil((((n_mx+1)*interval)/pi - 1.5) / 2);
                    wavestart_n = ceil((pi * (1.5 + i * 2))/interval);
                end

                %make the cosine wave 
                %(within the same period in the beginning-extend for 1.5 pi:2.5pi~4.0pi)
                %(within the first maximum value of the imf as the cosine amplitude)    
                x2 = 0:interval:((wavestart_n-1)*interval);%the x-axis value of new coordinate
                wavestart = cos(pi*2.5+x2) * abs(imf(n_mx));%the y-axis value of cosine function(2.5pi~4.0pi)
                
                %calculate weights for interpolation
                weight = 1/n_mx;
                weight_intr = weight;
                
                curr = wavestart_n;
                
                %interpolate between actual wave and artificial wave
                %means extend out gradually between original wave and artificiaal cosine wave
                for i=n_mx+1:-1:1
                    wavestart(curr) = (imf(i) * (1-weight)) + (wavestart(curr) * weight);
                    curr=curr-1;
                    weight = weight + weight_intr;%weighting is always increasing step by step
                end

                y = [wavestart, imf(n_mx+1:end)];%add to the original in the heading position
                thestart = n_mx;
            end % end of maximum and maximum condition processing 
        end     % end of extending
    end         %end of processing about the begginning part of imf  ----loop B end
    
    
    % when ex=0 means endding point needs rectification----loop B start
    % start of extending,when ex=0 means the data need extending extra segemnt
    if (ex == 0)
        %search again, this time flipped
        %find out maximum,minimum ,zero-crossing points in the endding part
        imf = fliplr(imf);
        zeroc2 = 0;

        n_mn = length(y);
        n_mx = length(y);
        flagmax = 0;
        flagmin = 0;

      %the conceptual explantion of the following sentance is almost the same with the former
      %so footnote will not go on,except for extra new syntax ,new algorithm
      %modification is still go on within the following sentences 

        for i=2:n-1
            %finds local maximum
            if ((imf(i) > imf(i-1)) & (flagmax == 0)) %might be an extrema
                flagmax = 1;
                Xmax = i;
            end
            if((imf(i) < imf(i+1)) & (flagmax == 1)) %not an extrema
                flagmax = 0;
            end
            if((imf(i) > imf(i+1)) & (flagmax == 1)) %is an extrema
                flagmax = 0;
                n_mx = Xmax;
            end

            %finds local minimum
            if((imf(i) < imf(i-1)) & (flagmin == 0)) %might be an extrema
                flagmin = 1;
                Xmin = i;
            end
            if((imf(i) > imf(i+1)) & (flagmin == 1)) %not an extrema
                flagmin = 0;
            end
            if((imf(i) < imf(i+1)) & (flagmin == 1)) %is an extrema
                flagmin = 0;

                n_mn = Xmin;
            end

            %finds zero crossing
            if (((imf(i) == 0) | (imf(i) > 0 & imf(i+1) < 0) | (imf(i) < 0 & imf(i+1) > 0)) & zeroc2 <= 0)
                zeroc2 = i;
            end

            if((zeroc2 > 0) & (n_mn < length(y) | n_mx < length(y)))
                     break;
            end
        end

        imf = fliplr(imf);

        zeroc2 = length(y)-zeroc2;
        n_mn = length(y)-n_mn;
        n_mx = length(y)-n_mx;
        

        %makes sure zeroc2+1 will function
        if(zeroc2 == length(y))
            zeroc2 = zeroc2-1;
        end
        
        if(((zeroc2 > 0) + (n_mn > 0) + (n_mx > 0)) < 2)
            h = [];
            if orient == 1
                hret = hret.';
            end

            return;
        end

        if(n_mn > n_mx) %if min is closer to end
            if(n_mn == zeroc2)
                zeroc2 = length(y)-1;
            end

            slope = (y(zeroc2+1) - y(zeroc2)); %(delta x = 1)
            b = y(zeroc2+1) - (slope * (zeroc2+1));
            x = -b/slope;
            
            %find out frequency interval 
            interval = abs(pi/((x-n_mn)*2));
            interval_b = interval;
            
            %find the length of a correct wave
            wavend_n = ceil((pi *1.5)/interval);
            

         
            if(wavend_n <= length(y)-n_mn)
                i = ceil((((length(y)-n_mn+1)*interval)/pi - 1.5) / 2);
                wavend_n = ceil((pi * (1.5 + i * 2))/interval);
            end

            x2 = 0:interval:(wavend_n-1)*interval;
            wavend = cos(-pi * 3 + x2) * abs(y(n_mn));

            weight = 1/(length(y)-n_mn);
            weight_intr = weight;
            weight = 1-weight;  %flip weight values
            
            for i=0:(length(y)-n_mn)
                wavend(i+1) = (y(i+n_mn) * weight) + (wavend(i+1)* (1 - weight));
                weight = weight - weight_intr;
            end

            thend = length(y)-n_mn+1;
            y = [y(1:n_mn-1), wavend];
            

        else    %if max is closer to end
            %        wavend = cos(-pi*4:pi/(abs(n_mx-zeroc2)*2):pi*4.5)*abs(y(n_mx));
            if(n_mx == zeroc2)
                zeroc2 = length(y)-1;
            end

            slope = (y(zeroc2+1) - y(zeroc2)); %(delta x = 1)
            b = y(zeroc2+1) - (slope * (zeroc2+1));
            x = -b/slope;
            
            %find out frequency interval   
            interval = abs(pi/((x-n_mx)*2));
            interval_b = interval;
            
            %find the length of a correct wave
            wavend_n = ceil((pi *2.5)/interval);

            if(wavend_n <= length(y)-n_mx)
                i = ceil((((length(y)-n_mx+1)*interval)/pi - 2.5) / 2);
                wavend_n = ceil((pi * (2.5 + i * 2))/interval);
            end
            
            x2 = 0:interval:(wavend_n-1)*interval;
            wavend = cos(-pi * 4 + x2) * abs(y(n_mx));
            
            weight = 1/(length(y)-n_mx);
            weight_intr = weight;
            weight = 1-weight;  %flip weight values
            
            for i=0:(length(y)-n_mx)
                wavend(i+1) = (y(n_mx+i) * weight) + (wavend(i+1)* (1 - weight));
                weight = weight - weight_intr;
            end
            
            thend = length(y)-n_mx+1;
            y = [y(1:n_mx-1), wavend];
        end
    end    %end of processing about the endding part of imf  ----loop B end

       
    %the checking figure plotting syntax
    if chkplot==1
        
        imforiginal=[zeros(1,(length(wavestart)-thestart)),imf,zeros(1,(length(wavend)-thend))];
        numberori=numel(imforiginal);
        numberext=numel(y);
        if numberori==numberext
            gg=1:numberori;
            figure(91)
            plot(gg,y,'r',gg,imforiginal,'b');title('After end-process of original IMF,the result');legend('Modified','Original');
            disp('program paused for figure checking,press any key to continue') 
            pause
        else
            disp('number not match,BAD,no figure,Check please!')
        end
    end   
    
    %DO Hilbert-Transform with matlab hilbert.m
    h = hilbert(y);
    
    %remove those additional segments out (in head and tail that we add  in)
    h = h(length(wavestart)-thestart+1:length(imf)+length(wavestart)-thestart);
    %use extra head length and original imf numbers can easily remove those
    %extra parts in the head and in the end
    
   %copy data into return array
    if (bx==1 & ex ==1)
        hret(k,:) = 0;
    else   
        hret(k,:) = h;
    end   
    clear h y
    %clear for correct calculation
end 

if orient == 1
    hret = hret.';
end

return;        
