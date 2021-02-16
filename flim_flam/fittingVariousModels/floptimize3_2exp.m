%This function fits the data in decay to an unconstrained 2 exponential
%decay model described by:

% F(tau1, tau2, a1, a2, t) = a1*exp(-t/tau1) + a2*exp(-t/tau2)

%Code was originally adapted (although heavily modified) from open source code
%orginally documented in:
%   Enderlein, J; Erdmann, R, Optics Communications, 1997, 134, 371-378.
%All edits were made by Julia Lazzari-Dean.
%Last edits 10/11/2019

%Dependencies: convol.m, sortATs.m, assignParams_2exp.m

%Input parameters:
% decay - time resolved fluorescence signal (one decay at a time), as an N x
%  1 array (N is the ADC resolution)
% irf - measured instrument response, also as an Nx1 array (I recommend
%  trimming only to time bins of interest, e.g. 26-36 out of 256 for the
%  current MaiTai setup). The fit guess is reconvolved with the irf
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

%The fit guess is reconvolved with the irf after application of a color
%shift. The color shift start point is given by the cShift parameter. The
%boolean shiftFixed dictates whether the shift will be included as a free
%parameter in the fit (0 = shift is a free parameter, 1 = shift is fixed to
%starting value). 

%The offFixed boolean dictates whether the offset is fit (offFixed = 0) or
%whether it is fixed to 0 during the analysis (offFixed = 1). There is not
%currently support for setting offset to a particular value here - instead
%subtract the offset from each decay in the upstream parsing code.

%The function uses the fmincon built-in Matlab function for optimization
%and minimizes the reduced chi squared of the fit, specifically calculated
%as the following (where st and fi are the time bins of the start and
%finish of the fitting, respectively). Note: This is a LOCAL NOT A GLOBAL
%optimization, selected to decrease run time.

function [tm, aFs, tFs, cF, offset, chiSq, residTrace, SSE, exitFlag] = floptimize3_2exp(decay,irf,cShift,shiftFixed,offFixed,startParams,fixedParams,stFi,period)

%initialize parameters correctly for deckard usage
p = period; %time range
n = length(irf);
t = (1:n)'; %time scale in time bin units
dt = p/n;
tp = (0:dt:(p-dt))'; %time scale
st = stFi(1);
fi = stFi(2);
residTrace = zeros(size(tp));

%fixedParam should be all 0s or 1s
for i=1:size(fixedParams,2)
    if(fixedParams(1,i) ~= 0 && fixedParams(1,i) ~= 1)
        error('All values in the fixed param array must be either zero or one.');
    end
end

%startParams will be a vector with the desired starting values for each
%parameter. It is paired with the vector fixedParams, which will specify
%which of the parameters is fixed. 
nParam = length(startParams) - nnz(fixedParams);
paramToUse = ~fixedParams;
param = startParams .* paramToUse;
param(param == 0) = [];
%if the shift will be free, add it to the end of the param vector
if ~shiftFixed
    param(1,end+1) = cShift;
    nParam = nParam + 1;
end

%setup for fmincon based optimization
options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
    'MaxFunctionEvaluations',1e5,'StepTolerance',1e-6,'ConstraintTolerance',1e-6);
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
        [a1G, a2G, t1G, t2G, uP] = assignParams_2exp(param,startParams,fixedParams);
        
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
        irs = (1-cG+floor(cG))*irf(rem(rem(t-floor(cG)-1, n)+n,n)+1) + (cG-floor(cG))*irf(rem(rem(t-ceil(cG)-1, n)+n,n)+1);
        z = convol(irs, x);
        
        %fix the offset coefficient or leave it floating
        if(offFixed)
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
[a1F, a2F, t1F, t2F, uPF] = assignParams_2exp(pOpt,startParams,fixedParams);

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
irs = (1-cF+floor(cF))*irf(rem(rem(t-floor(cF)-1, n)+n,n)+1) + (cF-floor(cF))*irf(rem(rem(t-ceil(cF)-1, n)+n,n)+1);
gC = convol(irs, x);
%fix the offset coefficient or leave it floating
if(offFixed)
    gC = [zeros(size(gC,1),1) gC];
else
    gC = [ones(size(gC,1),1) gC];
end
coeff = lsqnonneg(gC,decay);
gCW = gC*coeff;

%from lsqnonneg - coeff(1) is the offset and coeff(2)is the scale factor
%between the two - I THINK.
offset = coeff(1);

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
	%SSE = SSE + (decay(i,1)-gCW(i,1)).^2/abs(gCW(i,1));
    SSE = SSE + (decay(i,1)-gCW(i,1)).^2;
end

%calculate the weighted residuals
residTrace(st:fi) = (decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi)));      

end

