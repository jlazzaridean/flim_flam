%This function allows batch analysis of multiple decays using the
%exponential series method. It also combines the assorted fit outputs into
%one useful output structure to avoid variable clutter.

function [results] = batchESM(decays,IRF,tMax,nSteps,cShift,shiftFixed,offFixed)

results = struct('decay',0,'IRF',0,'taus',0,'aFinal',0,'shift',0,'offset',0,...
    'chiSq',0,'AIC',0,'exitFlag',0,'residTrace',0);

tSize = tMax/nSteps;
taus = (tSize:tSize:tMax)';

for i=1:size(decays,2)
    results(i,1).decay = decays(:,i);
    results(i,1).IRF = IRF;
    results(i,1).taus = taus;
    
    [results(i,1).aFinal, results(i,1).chiSq, results(i,1).shift,...
        results(i,1).offset, results(i,1).residTrace, results(i,1).AIC,...
        results(i,1).exitFlag] = floptimize3_ESM(decays(:,i),...
        IRF,tMax,nSteps,cShift,shiftFixed,offFixed);
end

end

