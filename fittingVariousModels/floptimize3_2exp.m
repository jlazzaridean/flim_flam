%This function fits the data in decay to an unconstrained 2 exponential
%decay model described by:

% F(tau1, tau2, a1, a2, t) = a1*exp(-t/tau1) + a2*exp(-t/tau2)

%Dependencies: convol.m, sortATs.m, assignParams_2exp.m
%Note: convol.m function is used from open source code
%orginally documented in:
%   Enderlein, J; Erdmann, R, Optics Communications, 1997, 134, 371-378.

%Input parameters:
% decay - time resolved fluorescence signal (one decay at a time), as an N x
%  1 array (N is the ADC resolution)
% IRF - measured instrument response, also as an Nx1 array (I recommend
%  trimming only to time bins of interest, e.g. 26-36 out of 256 for the
%  current MaiTai setup). The fit guess is reconvolved with the IRF
%  after application of a color shift (below) before chiSq is evaluated.
% cShift - starting value of the color shift parameter.
% shiftFixed - if 0, will optimize cShift. If 1, will fix the shift to
%  whatever was passed as cShift.
% dcOffset - time-independent background in the FLIM decay. The way this
%  data was taken, that is effectively 0 (and the code doesn't currently
%  support any other values). Signal before the laser pulse will be
%  interpreted as arising from the previous laser cycle.

%The function returns the following:
% tm - amplitude weighted average decay time (tm) (a1*t1 + a2*t2, where a1 +
% a2 = 1
% aF - weights of each time constant, in the order of the time constants
% tF - fitted decay times, in ascending order.
% cF - color shift, either the fixed value or the fitted value
% offset - should be 0 in this implementation.
% chiSq - goodness of fit, see below.
% exitFlag - the stopping condition from the Matlab optimization routine.

% chiSq = sum((guess(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);
%(Difference between the measured and the fit value for all time bins used,
% weighted by the number of photon counts in that time bin, as is typical
% for Poisson processes.

%The fit guess is reconvolved with the IRF after application of a color
%shift. The color shift start point is given by the cShift parameter. The
%boolean shiftFixed dictates whether the shift will be included as a free
%parameter in the fit (0 = shift is a free parameter, 1 = shift is fixed to
%starting value). 

%The offsetMode value determines how the offset is found.
    %offsetMode 0: offset is floating (determined during the fit)
    %offsetMode 1: offset is fixed to 0 during analysis (assumes no dark counts)
    %offsetMode 2: offset is determined from the time bins before the laser
        %pulse. This value is subtracted from all data, and then the offset
        %is fixed to 0 during fitting. If there is an incomplete decay (i.e.
        %residual signal from the previous laser cycle) in the early bins, 
        %offsetMode 2 will give you artefactual results.

%The function uses the fmincon built-in Matlab function for optimization
%and minimizes the reduced chi squared of the fit, specifically calculated
%as the following (where st and fi are the time bins of the start and
%finish of the fitting, respectively). Note: This is a LOCAL NOT A GLOBAL
%optimization, selected to decrease run time.

%Last edits: 1/30/2020, Julia Lazzari-Dean

function [tm, aFs, tFs, cF, offset, chiSq, residTrace, SSE, exitFlag] = floptimize3_2exp(decay,IRF,configS)

%extract relevant parameters from the configuration structure
p = configS(1,1).nsPeriod; %time range
n = length(IRF);
t = (1:n)'; %time scale in time bin units
dt = p/n;
tp = (0:dt:(p-dt))'; %time scale 
st = configS(1,1).stFi(1)/configS(1,1).adcBin;
fi = configS(1,1).stFi(2)/configS(1,1).adcBin;
residTrace = zeros(size(tp));
cShift = configS(1,1).cShift;
startParam = configS(1,1).startParam;
fixedParam = configS(1,1).fixedParam;
offsetMode = configS(1,1).offsetMode;
shiftFixed = configS(1,1).shiftFixed;

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

%fixedParam should be all 0s or 1s
for i=1:size(fixedParam,2)
    if(fixedParam(1,i) ~= 0 && fixedParam(1,i) ~= 1)
        error('All values in the fixed param array must be either zero or one.');
    end
end

%startParam will be a vector with the desired starting values for each
%parameter. It is paired with the vector fixedParam, which will specify
%which of the parameters is fixed. 
nParam = length(startParam) - nnz(fixedParam);
paramToUse = ~fixedParam;
param = startParam .* paramToUse;
param(param == 0) = [];
%if the shift will be free, add it to the end of the param vector
if ~shiftFixed
    param(1,end+1) = cShift;
    nParam = nParam + 1;
end

%setup for fmincon - extract relevant values from the config file
maxFunEval = configS(1,1).maxFunEval;
stepTol = configS(1,1).stepTol;
optimTol = configS(1,1).optimTol;
constraintTol = configS(1,1).constraintTol;

%create a structure with fit options and bounds for the function call
options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
    'MaxFunctionEvaluations',maxFunEval,'StepTolerance',stepTol,...
    'OptimalityTolerance',optimTol,'ConstraintTolerance',constraintTol);
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

    function [chiSquared] = expScoring(param)
        %	This nested function calculates the scoring function
        %this function has access to all of the values within floptimize.
        
        %assign the current guesses and fixed parameters
        [a1G, a2G, t1G, t2G, uP] = assignParams_2exp(param,startParam,fixedParam);
        
        %fix the shift as is appropriate
        if(~shiftFixed)
            cG = param(nParam); %either as the last parameter in the fitted values
            uP = uP + 1;
        else
            cG = cShift; %...or as the original input
        end
        
        if (uP ~= nParam)
            error('Incorrect matching of parameters in fitting function.');
        end
        
        %recalculate the guesss
        x = a1G*exp(-tp/t1G) + a2G*exp(-tp/t2G);
        %shift the IRF by the current guess for the shift amount and
        %reconvolve with the current exponential guess
        irs = (1-cG+floor(cG))*IRF(rem(rem(t-floor(cG)-1, n)+n,n)+1) + (cG-floor(cG))*IRF(rem(rem(t-ceil(cG)-1, n)+n,n)+1);
        z = convol(irs, x);
        
        %fix the offset coefficient or leave it floating
        if offsetMode == 1 || offsetMode == 2
            z = [zeros(size(z,1),1) z];
        else
            z = [ones(size(z,1),1) z];
        end
        
        %scale the hypothesized fit to the decay
        scale = lsqnonneg(z,decay);
        z = z*scale;
        
        %%recalculate the scoring function
        chiSquared = sum((z(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);
    end

%extract the final values from the output
[a1F, a2F, t1F, t2F, uPF] = assignParams_2exp(pOpt,startParam,fixedParam);

%assign the shift
if(~shiftFixed)
    cF = pOpt(nParam); %either as the last parameter in the fitted values
    uPF = uPF + 1;
else
    cF = cShift; %...or as the original input
end

%sanity check for parameter assignment
if (uPF ~= nParam)
    error('Incorrect matching of parameters in final assignment.');
end

%recalculate the final state of the function
x = a1F*exp(-tp/t1F) + a2F*exp(-tp/t2F);

%shift the IRF by the current guess for the shift amount and
%reconvolve with the current exponential guess
irs = (1-cF+floor(cF))*IRF(rem(rem(t-floor(cF)-1, n)+n,n)+1) + (cF-floor(cF))*IRF(rem(rem(t-ceil(cF)-1, n)+n,n)+1);
gC = convol(irs, x);
%fix the offset coefficient or leave it floating
if(offsetMode == 1 || offsetMode == 2)
    gC = [zeros(size(gC,1),1) gC];
else
    gC = [ones(size(gC,1),1) gC];
end
coeff = lsqnonneg(gC,decay);
gCW = gC*coeff;

%from lsqnonneg - i believe coeff(1) is the offset and coeff(2)is the scale factor
%between the two
if offsetMode == 0
    offset = coeff(1);
elseif offsetMode == 1
    offset = 0;
end

%calculate tm at  the final state
tm = a1F*t1F + a2F*t2F;

%put the taus and coefficients into arrays so they can be returned easily
aF = [a1F a2F];
tF = [t1F t2F];
%sort them such that the taus are in ascending order for ease of processing
%later.
[aFs,tFs] = sortATs(aF, tF);

SSE=0;
for i=st:fi
    SSE = SSE + (decay(i,1)-gCW(i,1)).^2;
end

%calculate the weighted residuals
residTrace(st:fi) = (decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi)));      

if configS(1,1).viewDecay
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
    
    %program will wait for the user key press to acknowledge having viewed
    %the decay
    pause;
end

end

