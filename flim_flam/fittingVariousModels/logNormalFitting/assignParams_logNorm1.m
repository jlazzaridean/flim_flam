function [tauB, tauS, usedParams] = assignParams_logNorm1(pFit, pStart, pFixed)

usedParams = 0;

if(pFixed(1,1) == 0) % tau bar
    tauB = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    tauB = pStart(1);
end
if(pFixed(1,2) == 0) % now tau sigma
    tauS = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    tauS = pStart(2);
end

end