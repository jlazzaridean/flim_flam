%This function takes a set of decays and an IRF (as well as a variety of
%other fit parameters) and performs a single fit per decay (i.e. global
%analysis). It returns a structure of the results.

%need to comment input data format more thoroughly

function [results] = batchFloptimizeGeneral(decays,IRF,model,cShift,shiftFixed,dcOffset,startParams,fixedParams,stFi,period_ns)

if(strcmp(model,'LogNorm1'))
    %set up a new structure
    results = struct('decay',0,'tauBar',0,'tauSig',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'exitFlag',0,'geoMean',0,'geoStd',0,'residTrace',0,...
        'SSE',0,'startParams',0,'fixedParams',0);
    
    for i = 1:size(decays,2)
        results(i,1).decay = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParams = startParams;
        results(i,1).fixedParams = fixedParams;
        [results(i,1).tauBar, results(i,1).tauSig, results(i,1).geoMean,...
            results(i,1).geoStd, results(i,1).chiSq, results(i,1).shift,...
            results(i,1).offset, results(i,1).residTrace,results(i,1).SSE,...
            results(i,1).exitFlag] = floptimize3_logNorm1(results(i,1).decays,...
            IRF,cShift,shiftFixed,dcOffset,startParams,fixedParams,stFi,period_ns);
    end
    %note that the order for start params in lognorm1 is taubar, tau sigma
elseif(strcmp(model,'LogNorm2'))
    %set up a new structure
    results = struct('decays',0,'tm',0,'a',0,'tauBar',0,'tauSig',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'exitFlag',0,'geoMean',0,'geoStd',0,'residTrace',0,...
        'SSE',0,'startParams',0,'fixedParams',0);
    
    for i = 1:size(decays,2)
        results(i,1).decay = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParams = startParams;
        results(i,1).fixedParams = fixedParams;
        
        [results(i,1).tm, results(i,1).a,results(i,1).tauBar, results(i,1).tauSig, results(i,1).geoMean,...
            results(i,1).geoStd, results(i,1).chiSq, results(i,1).shift,...
            results(i,1).offset, results(i,1).residTrace,results(i,1).SSE,...
            results(i,1).exitFlag] = floptimize3_logNorm2(results(i,1).decays,...
            IRF,cShift,shiftFixed,dcOffset,startParams,fixedParams,stFi,period_ns);
    end
    %note that the order for start params in lognorm2 is A (x2), tauBar
    %(x2), tauSig (x2)
elseif(strcmp(model,'1exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'shift',0,'IRF',0,'chiSq',0,...
        'offset',0,'SSE',0,'residTrace',0,'exitFlag',0,'startParams',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParams = startParams;
        
        [results(i,1).tau, results(i,1).shift, results(i,1).offset,...
            results(i,1).chiSq, results(i,1).residTrace,results(i,1).SSE,...
            results(i,1).exitFlag] = floptimize3_1exp(results(i,1).decays,...
            IRF,cShift,shiftFixed,dcOffset,startParams,stFi,period_ns);
    end
    
elseif(strcmp(model,'2exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'a',0,'tm',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParams',0,'fixedParams',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParams = startParams;
        results(i,1).fixedParams = fixedParams;
        
        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq,results(i,1).residTrace,...
            results(i,1).SSE, results(i,1).exitFlag] = floptimize3_2exp(results(i,1).decays,IRF,cShift,...
            shiftFixed,dcOffset,startParams,fixedParams,stFi,period_ns);
    end        
    
elseif(strcmp(model,'3exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'a',0,'tm',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParams',0,'fixedParams',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParams = startParams;
        results(i,1).fixedParams = fixedParams;

        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq, results(i,1).residTrace, ...
            results(i,1).SSE, results(i,1).exitFlag] = floptimize3_3exp(results(i,1).decays,...
            IRF,cShift,shiftFixed,dcOffset,startParams,fixedParams,stFi,period_ns);
    end    
    
elseif(strcmp(model,'4exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'a',0,'tm',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParams',0,'fixedParams',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParams = startParams;
        results(i,1).fixedParams = fixedParams;

        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq, results(i,1).residTrace,...
            results(i,1).SSE,results(i,1).exitFlag] = floptimize3_4exp(results(i,1).decays,...
            IRF,cShift,shiftFixed,dcOffset,startParams,fixedParams,stFi,period_ns);
    end   
else
    err('Model string entered is currently unrecognized.');
end

%may want to document start point, etc.

end

