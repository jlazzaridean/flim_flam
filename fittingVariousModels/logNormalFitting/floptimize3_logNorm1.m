%order of start params is tau bar, then tau sigma.

function [tauBF, tauSF, geoMean, geoStd, chiSq, cF, offset, residTrace, SSE, exitFlag] = floptimize3_logNorm1(decay,irf,cShift,shiftFixed,dcOffset,startParams,fixedParams)

%initialize parameters correctly for deckard usage
p = 12.5; %time range
n = length(irf);
t = (1:n)'; %time scale in time bin units
dt = p/n;
tp = (0:dt:(p-dt))'; %time scale 
st = 20;
fi = 240;
residTrace = zeros(size(tp));

%startParams will be a vector with the desired starting values for each
%parameter. It is paired with the vector fixedParams, which will specify
%which of the parameters is fixed. If a value is listed as "1" it is fixed,
%if it is listed as "0" it is free.
if (shiftFixed == 0)
    %shift is a free parameter for the function now - start it at 0
    nParam = length(startParams) - nnz(fixedParams) + 1;
    param = -1;
    for i = 1:length(startParams)
        if(fixedParams(1,i) == 0)
            %if we will be fitting this value, put the start point in the
            %vector to go into the fitting function.
            if(param == -1)
                param(1,1) = startParams(1,i);
            else
                param(1,end+1) = startParams(1,i);
            end
        elseif (fixedParams(1,i) ==1)
            %do nothing for now - will find the set values later.
        else
            error('Values in the fixedParams vector should be 0 or 10.');
        end
    end
    param(1,end+1) = cShift;
    
elseif (shiftFixed == 1)
    %shift is a fixed parameter for the function now - don't add it to
    %parameter vector
    nParam = length(startParams) - nnz(fixedParams);
    param = -1;
    for i = 1:length(startParams)
        if(fixedParams(1,i) == 0)
            %if we will be fitting this value, put the start point in the
            %vector to go into the fitting function.
            if(param == -1)
                param(1,1) = startParams(1,i);
            else
                param(1,end+1) = startParams(1,i);
            end
        elseif (fixedParams(1,i) ==1)
            %do nothing for now - will find the set values later.
        else
            error('Values in the fixedParams vector should be 0 or 10.');
        end
    end
else
    error('shiftFixed boolean should be either 0 or 1');
end

%parse the input parameters to initialize offset
%if dcOffset = -1, let the shift float
%if dcOffset is any other value, it will fix the offset to zero (HACKED RN)
if (dcOffset == -1)
    offFixed = 0; %allow the offset to float during the fitting 
else
    offFixed = 1; %fix the shift to the input value during the fitting
    %HACKED - THIS DOESN'T WORK -- WILL FIX IT TO ZERO BY DEFAULT
end

%setup for fmincon
options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
    'MaxFunctionEvaluations',1e5,'StepTolerance',1e-10,'ConstraintTolerance',1e-10);
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(1,nParam);
ub = zeros(1,nParam);
if (shiftFixed)
    for i=1:nParam
        ub(1,i) = 10;
        lb(1,i) = -10;
    end
else
    for i=1:(nParam-1)
        ub(1,i) = 10;
        lb(1,i) = -10;
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
        [tauBG, tauSG, uP] = assignParams_logNorm1(param,startParams,fixedParams);
        
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
        x = evalGaussFunc_noIntegral(tp,tauBG,tauSG);
        %shift the IRF by the current guess for the shift amount and
        %reconvolve with the current exponential guess
        %evaulate the result of the function with the current settings
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
     
        %recalculate the scoring function
        resid = sum((z(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);            
    end

%extract the final values from the output
[tauBF, tauSF, uPF] = assignParams_logNorm1(pOpt,startParams,fixedParams);

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
x = evalGaussFunc_noIntegral(tp,tauBF,tauSF);

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

%calculate geoMean and geoStd (more meaningful parameters)
geoMean = exp(tauBF);
geoStd = exp(tauSF);

%calculate the weighted residuals
residTrace(st:fi) = (decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi)));

SSE=0;
for i=st:fi
	%SSE = SSE + (decay(i,1)-gCW(i,1)).^2/abs(gCW(i,1));
    SSE = SSE + (decay(i,1)-gCW(i,1)).^2;
end

end

