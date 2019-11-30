%This function fits the data in decay to single exponential decay,
%described by the following model:

% F(tau, t) = A*exp(-t/tau)

%The function returns a decay time, a color shift, an offset,
%the final reduced chi squared, and an exit flag as to whether the
%algorithm converged (positive values are good, negative values are bad).

%The fit guess is reconvolved with the irf after application of a color
%shift. The color shift start point is given by the cShift parameter. The
%boolean shiftFixed dictates whether the shift will be included as a free
%parameter in the fit (0 = shift is a free parameter, 1 = shift is fixed to
%starting value). 

%The offsetMode parameter dictates whether the offset is fit (offFixed = 0) or
%whether it is fixed during the analysis (offsetMode = 1,2). If
%offsetMode = 1, the offset is fixed to 0. If offsetMode = 2, the offset is
%determined from the time bins before the start of the fitting and after
%the end of the fitting.

%The function uses the fmincon built-in Matlab function for optimization
%and minimizes the reduced chi squared of the fit, specifically calculated
%as the following (where st and fi are the time bins of the start and
%finish of the fitting, respectively). Note: This is a LOCAL NOT A GLOBAL
%optimization, selected to decrease run time.

% chiSq = sum((guess(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);  
%(Each value is weighted by the number of photon counts)

%Default parameters that are perhaps relevant to the astute user:
%Start point: all coefficients equal, decay constants: 0.5, 1.5, 2.5 ns
%256 time bins of ADC resolution, fitting only to time bins 23-240, laser
%period 12.5 ns
%This script assumes the IRF has already been trimmed to the desired time
%bins (such as 26-36 of the ADC range).

%last edits 1/24/2019 Julia Lazzari-Dean


function [tF, cF, offset, chiSq, residTrace, SSE, exitFlag] = floptimize3_1exp(decay,irf,cShift,shiftFixed,offsetMode,startParam,stFi,period,viewDecay)

%initialize parameters correctly for deckard usage
p = period; %time range
n = length(irf);
t = (1:n)'; %time scale in time bin units
dt = p/n;
tp = (0:dt:(p-dt))'; %time scale 
st = floor(stFi(1));
fi = floor(stFi(2));
tStart = startParam; %the start point by default for lifetime analysis
residTrace = zeros(size(tp));

%find index of first zero in the decay (after the peak, use this as fi if it is less than
%fi
zeroBins = decay == 0;
zeroInd = -1;
[~, peakInd] = max(decay);
for i = peakInd:n
    if zeroBins(i,1) == 1
        zeroInd = i;
        break;
    end
end
if (zeroInd > 0 && zeroInd < fi)
    fi = zeroInd;
end

%calculate the offset from the values before the start and after the finish
%if offset subtraction mode is 2
if offsetMode == 2
    bkgd1 = decay(1:st);
    offset = mean(bkgd1);
    decay = decay - offset;
    decay(decay<0) = 0;
end

%decide whether shift should be fixed or free based on the value of
%shiftFixed boolean
if (shiftFixed == 0)
    %shift is a free parameter for the function now
    param = [tStart cShift]; 
    nParam = 2; %number of free params we are fitting to
elseif (shiftFixed == 1)
    c = cShift;
    param = tStart; %the only parameter is starting lifetime
    nParam = 1;
else
    error('shiftFixed boolean should be either 0 or 1');
end

%setup for fmincon
options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
    'MaxFunctionEvaluations',1e5,'StepTolerance',1e-4,'OptimalityTolerance',1e-4,...
    'ConstraintTolerance',1e-6);
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(1,nParam);
ub = zeros(1,nParam);
if (shiftFixed)
    for i=1:nParam
        ub(1,i) = 10;
        lb(1,i) = 0;
    end
else
    for i=1:(nParam-1)
        ub(1,i) = 10;
        lb(1,i) = 0;
    end
    %bound shift between +/- 5 time bin units
    ub(1,nParam) = 5;
    lb(1,nParam) = -5;
end
nonlcon = [];
[pOpt,chiSq,exitFlag] = fmincon(@expScoring,param,A,b,Aeq,beq,lb,ub,nonlcon,options);

    function [resid] = expScoring(param)
        %this nested function calculates the scoring function        
        %this function has access to all of the values within floptimize.
                     
        if (shiftFixed == 1)
            %evalulate the function again with the current guess
            tGuess = param;
            x = exp(-tp/tGuess);
            %shift the IRF by the designated amount and reconvolve with the
            %result
            irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
            z = convol(irs, x);
        elseif (shiftFixed == 0)
            %evalulate the function again with the current guess
            tGuess = param(1);
            cG = param(2);
            x = exp(-tp/tGuess);
            %shift the IRF by the current guess for the shift amount and 
            %reconvolve with the current exponential guess
            irs = (1-cG+floor(cG))*irf(rem(rem(t-floor(cG)-1, n)+n,n)+1) + (cG-floor(cG))*irf(rem(rem(t-ceil(cG)-1, n)+n,n)+1);
            z = convol(irs, x);            
        else
            error('Boolean for shift is not set correctly.');
        end
        
        %fix the offset coefficient or leave it floating
        if offsetMode == 1 || offsetMode == 2
            z = [zeros(size(z,1),1) z];
        else
            z = [ones(size(z,1),1) z];
        end
        
        %scale the hypothesized fit to the decay
        scale = lsqnonneg(z,decay);
        z = z*scale;
             
        %recalculate the scoring function
        resid = sum((z(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);        
    end

%extract the final values from the output
if (shiftFixed == 1)
    %evalulate the function again with the current guess
    tF = pOpt;
    cF = c;
elseif (shiftFixed == 0)
    %evalulate the function again with the current guess
    tF = pOpt(1);
    cF = pOpt(2);
else
    error('Boolean for shift is not set correctly.');
end

%recalculate the final state of the function
x = exp(-tp/tF);
%shift the IRF by the current guess for the shift amount and
%reconvolve with the current exponential guess
irs = (1-cF+floor(cF))*irf(rem(rem(t-floor(cF)-1, n)+n,n)+1) + (cF-floor(cF))*irf(rem(rem(t-ceil(cF)-1, n)+n,n)+1);
gC = convol(irs, x);
%fix the offset coefficient or leave it floating
if(offsetMode == 1 || offsetMode == 2)
    gC = [zeros(size(gC,1),1) gC];
else
    gC = [ones(size(gC,1),1) gC];
end
coeff = lsqnonneg(gC,decay);
gCW = gC*coeff;

%from lsqnonneg - coeff(1) is the offset and coeff(2)is the scale factor
%between the two - I THINK.
if offsetMode == 0
    offset = coeff(1);
elseif offsetMode == 1
    offset = 0;
end

SSE = 0; %calculate the sum of squared errors
for i=st:fi
    SSE = SSE + (decay(i,1)-gCW(i,1)).^2;
end

%calculate the weighted residuals
residTrace(st:fi) = (decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi)));

if viewDecay
    %this is just here for troubleshooting
    %make a plot of the guesses versus the data
    hold off
    subplot('position',[0.1 0.4 0.8 0.5])
    plot(t,log10(decay),t,log10(irs./max(irs).*max(decay)),t,log10(gCW));
    v = axis;
    v(1) = min(t);
    v(2) = max(t);
    axis(v);
    xlabel('Time in ns');
    ylabel('Log Count');
    s = sprintf('Shift, Offset = %3.3f  %3.3f',cF,offset);
    text(max(t)/4*3,v(4)-0.05*(v(4)-v(3)),s);
    s = 'Lifetime = ';
    for i=1:length(tF)
        s = [s sprintf('%3.3f',tF(i)) '   '];
    end
    text(max(t)/4*3,v(4)-0.19*(v(4)-v(3)),s);
    subplot('position',[0.1 0.1 0.8 0.2])
    plot(t(st:fi),residTrace(st:fi));
    v = axis;
    v(1) = min(t);
    v(2) = max(t);
    
    axis(v);
    xlabel('Time in ns');
    ylabel('Residue');
    s = sprintf('%3.3f', chiSq);
    text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
    set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])
end

end

