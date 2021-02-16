%Helper function compares starting and fixed parameters for 4 exponential fits.

function [a1, a2, t1, t2, usedParams] = assignParams_2exp(pFit, pStart, pFixed)

usedParams = 0;
%pull out the coefficients first
if(pFixed(1,1) == 0) %a1
    a1 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    a1 = pStart(1);
end
if(pFixed(1,2) == 0) %a2
    a2 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    a2 = pStart(2);
end

%renormalize coefficients
aTot = a1 + a2;
a1 = a1/aTot;
a2 = a2/aTot;

%now pull out the taus
if(pFixed(1,3) == 0) %tau1
    t1 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    t1 = pStart(3);
end
if(pFixed(1,4) == 0) %tau2
    t2 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    t2 = pStart(4);
end

end

