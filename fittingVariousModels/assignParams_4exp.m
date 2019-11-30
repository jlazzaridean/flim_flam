%Helper function compares starting and fixed parameters for 4 exponential fits.

function [a1, a2, a3, a4, t1, t2, t3, t4, usedParams] = assignParams_4exp(pFit, pStart, pFixed)

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
if(pFixed(1,3) == 0) %a3
    a3 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    a3 = pStart(3);
end
if(pFixed(1,4) == 0) %a4
    a4 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    a4 = pStart(4);
end

%renormalize coefficients
aTot = a1 + a2 + a3 + a4;
a1 = a1/aTot;
a2 = a2/aTot;
a3 = a3/aTot;
a4 = a4/aTot;

%now pull out the taus
if(pFixed(1,5) == 0) %tau1
    t1 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    t1 = pStart(5);
end
if(pFixed(1,6) == 0) %tau2
    t2 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    t2 = pStart(6);
end
if(pFixed(1,7) == 0) %tau3
    t3 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    t3 = pStart(7);
end
if(pFixed(1,8) == 0) %tau4
    t4 = pFit(usedParams + 1);
    usedParams = usedParams + 1;
else
    t4 = pStart(8);
end

end

