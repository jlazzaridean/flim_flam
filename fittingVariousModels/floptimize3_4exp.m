%This function fits the data in decay to an unconstrained 3 exponential
%decay model described by:

% F(tau1, tau2, tau3, tau4, a1, a2, a3, a4, t)
%   = a1*exp(-t/tau1) + a2*exp(-t/tau2) + a3*exp(-t/tau3) + a4*exp(-t/tau4);

%The function returns a weighted average decay time (tm), an array of
%weights (aF), an array of decay constants (tF), a color shift, an offset,
%the final reduced chi squared, and an exit flag as to whether the
%algorithm converged (positive values are good, negative values are bad).

%The fit guess is reconvolved with the IRF after application of a color
%shift. The color shift start point is given by the cShift parameter. The
%boolean shiftFixed dictates whether the shift will be included as a free
%parameter in the fit (0 = shift is a free parameter, 1 = shift is fixed to
%starting value). 

%The offFixed boolean dictates whether the offset is fit (offsetMode = 0) or
%whether it is fixed to 0 during the analysis (offsetMode = 1). To subtract
%the offset as the average of a few time bins before the start of the
%decay, use offsetMode = 2.

%The function uses the fmincon built-in Matlab function for optimization
%and minimizes the reduced chi squared of the fit, specifically calculated
%as the following (where st and fi are the time bins of the start and
%finish of the fitting, respectively). Note: This is a LOCAL NOT A GLOBAL
%optimization, selected to decrease run time.

% chiSq = sum((guess(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);  
%(Each value is weighted by the number of photon counts)

%This script assumes the IRF has already been trimmed to the desired time
%bins (such as 26-36 of the ADC range).

%last edits 1/30/2020 Julia Lazzari-Dean

function [tm, aFs, tFs, cF, offset, chiSq, residTrace, SSE, exitFlag] = floptimize3_4exp(decay,IRF,configS)

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

    function [resid] = expScoring(param)
        %	This nested function calculates the scoring function  
        %this function has access to all of the values within floptimize.
         
        %match up the input parameters to model parameters based on the
        %pattern of fixed and free values.
        [a1G, a2G, a3G, a4G, t1G, t2G, t3G, t4G, uP] = assignParams_4exp(param,...
            startParam,fixedParam);
        
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
        x = a1G*exp(-tp/t1G) + a2G*exp(-tp/t2G) + a3G*exp(-tp/t3G) + a4G*exp(-tp/t4G);
        %shift the IRF by the current guess for the shift amount and
        %reconvolve with the current exponential guess
        irs = (1-cG+floor(cG))*IRF(rem(rem(t-floor(cG)-1, n)+n,n)+1) + (cG-floor(cG))*IRF(rem(rem(t-ceil(cG)-1, n)+n,n)+1);
        z = convol(irs, x);
        
        if offsetMode == 1 || offsetMode == 2
            %offset is fixed to zero at this stage (it either is zero or it
            %has already been subtracted based on a particular range)
            z = [zeros(size(z,1),1) z];
        else
            %determine the offset based on the "gap" between the fit and
            %the experimentally measured decay
            z = [ones(size(z,1),1) z];
        end
        
        %scale the hypothesized fit to the decay
        scale = lsqnonneg(z,decay);
        z = z*scale;       
     
        %recalculate the scoring function
        resid = sum((z(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);            
    end

%extract the final values from the output
[a1F, a2F, a3F, a4F, t1F, t2F, t3F, t4F,uPF] = assignParams_4exp(pOpt,startParam,fixedParam);

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
x = a1F*exp(-tp/t1F) + a2F*exp(-tp/t2F) + a3F*exp(-tp/t3F) + a4F*exp(-tp/t4F);

%shift the IRF by the current guess for the shift amount and
%reconvolve with the current exponential guess
irs = (1-cF+floor(cF))*IRF(rem(rem(t-floor(cF)-1, n)+n,n)+1) + (cF-floor(cF))*IRF(rem(rem(t-ceil(cF)-1, n)+n,n)+1);
gC = convol(irs, x);

%fix the offset coefficient or leave it floating
if offsetMode == 1 || offsetMode == 2
    gC = [zeros(size(gC,1),1) gC];
else
    gC = [ones(size(gC,1),1) gC];
end
coeff = lsqnonneg(gC,decay);
gCW = gC*coeff;

%from lsqnonneg - i believe coeff(1) is the offset and coeff(2)is the scale factor
%between the two
%from lsqnonneg - i believe coeff(1) is the offset and coeff(2)is the scale factor
%between the two
if offsetMode == 0
    offset = coeff(1);
elseif offsetMode == 1
    offset = 0;
end

%calculate tm at  the final state
tm = a1F*t1F + a2F*t2F + a3F*t3F + a4F*t4F;

%put the taus and coefficients into arrays so they can be returned easily
aF = [a1F a2F a3F a4F];
tF = [t1F t2F t3F t4F];
%sort them such that the taus are in ascending order for ease of processing
%later.
[aFs,tFs] = sortATs(aF, tF);

%calculate the weighted residuals
residTrace(st:fi) = (decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi)));

SSE=0; %calculate the sum of squared errors
for i=st:fi
    SSE = SSE + (decay(i,1)-gCW(i,1)).^2;
end

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

