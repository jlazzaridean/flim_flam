function [a1, a2, tauB1, tauB2, tauS1, tauS2, usedParams] = assignParams_logNorm2(pFit, pStart, pFixed)

usedParams = 0;

if(pFixed(1,1) == 0) %weight of first log norm dist
    a1 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    a1 = pStart(1);
end
if(pFixed(1,2) == 0) %weight of second log norm dist
    a2 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    a2 = pStart(2);
end

%normalize the coefficients
aTot = a1 + a2;
a1 = a1/aTot;
a2 = a2/aTot;

%now pull out tau bar and tau sigma for each distribution
if(pFixed(1,3) == 0) % tau bar 1
    tauB1 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    tauB1 = pStart(1);
end
if(pFixed(1,4) == 0) % tau bar 2
    tauB2 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    tauB2 = pStart(2);
end

if(pFixed(1,5) == 0) % tau sigma 1
    tauS1 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    tauS1 = pStart(1);
end
if(pFixed(1,6) == 0) % tau sigma 2
    tauS2 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    tauS2 = pStart(2);
end

end