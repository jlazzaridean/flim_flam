function [modAIC] = calcAIC(decay,guess,nParams,shiftFixed)

st = 20;
fi = 240;

if(~shiftFixed)
    error('Not implemented for variable shift yet!');
end

logLik = 0;
for i=st:fi
    %dPD = decay(i,1)*log(decay(i,1)/guess(i,1)) - (decay(i,1) -
    %guess(i,1)); Prendergast - this is the likelihood ratio chi squared
    dPD = decay(i,1)*log(guess(i,1)) - guess(i,1) - log(factorial(round(decay(i,1))));
    logLik = logLik + dPD;
end

%using nParams+1 because the variance of the residuals has to be a
%parameter? also implementing the "low number of samples AIC"
smallSampleCorr = (2*(nParams + 1)^2 + 2*(nParams + 1))/(fi-st-nParams-2);
modAIC = 2*(nParams+1) + (-2)*logLik + smallSampleCorr; 

modAIC = -1;

end

