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

%The offFixed boolean dictates whether the offset is fit (offFixed = 0) or
%whether it is fixed to 0 during the analysis (offFixed = 1). There is not
%currently support for setting offset to a particular value here - instead
%subtract the offset from each decay in the upstream parsing code.

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


function [tF, cF, offset, chiSq, residTrace, SSE, exitFlag] = floptimize3_1exp(decay,irf,cShift,shiftFixed,offFixed,startParam,stFi,period)

%initialize parameters correctly for deckard usage
p = period; %time range
n = length(irf);
t = (1:n)'; %time scale in time bin units
dt = p/n;
tp = (0:dt:(p-dt))'; %time scale 
st = stFi(1);
fi = stFi(2);
tStart = startParam; %the start point by default for lifetime analysis
residTrace = zeros(size(tp));

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
        if(offFixed)
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

SSE = 0; %calculate the sum of squared errors
for i=st:fi
    SSE = SSE + (decay(i,1)-gCW(i,1)).^2;
end

%calculate the weighted residuals
residTrace(st:fi) = (decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi)));      

end

