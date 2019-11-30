%Wrapper function for calling one of the floptimize3 fitting functions, as specified by the model.

%This function takes a set of decays and an IRF (as well as a variety of
%other fit parameters) and performs a single fit per decay (i.e. global
%analysis). It returns a structure of the results.

function [results] = batchFloptimizeExp(decays,IRF,model,cShift,shiftFixed,offFixed,startParam,fixedParam,stFi,period_ns,viewDecay)

if(strcmp(model,'1exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'shift',0,'IRF',0,'chiSq',0,...
        'offset',0,'SSE',0,'residTrace',0,'exitFlag',0,'startParam',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParam = startParam;
        results(i,1).fixedParam = fixedParam;
        
        [results(i,1).tau, results(i,1).shift, results(i,1).offset,...
            results(i,1).chiSq, results(i,1).residTrace,results(i,1).SSE,...
            results(i,1).exitFlag] = floptimize3_1exp(results(i,1).decays,...
            IRF,cShift,shiftFixed,offFixed,startParam,stFi,period_ns,viewDecay);
    end
    
elseif(strcmp(model,'2exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'a',0,'tm',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParam',0,'fixedParam',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParam = startParam;
        results(i,1).fixedParam = fixedParam;
        
        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq,results(i,1).residTrace,...
            results(i,1).SSE, results(i,1).exitFlag] = floptimize3_2exp(results(i,1).decays,IRF,cShift,...
            shiftFixed,offFixed,startParam,fixedParam,stFi,period_ns,viewDecay);
    end        
    
elseif(strcmp(model,'3exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'a',0,'tm',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParam',0,'fixedParam',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParam = startParam;
        results(i,1).fixedParam = fixedParam;

        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq, results(i,1).residTrace, ...
            results(i,1).SSE, results(i,1).exitFlag] = floptimize3_3exp(results(i,1).decays,...
            IRF,cShift,shiftFixed,offFixed,startParam,fixedParam,stFi,period_ns,viewDecay);
    end    
    
elseif(strcmp(model,'4exp'))
    %set up a new structure
    results = struct('decays',0,'tau',0,'a',0,'tm',0,'shift',0,'IRF',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParam',0,'fixedParam',0);    

    for i = 1:size(decays,2)
        results(i,1).decays = decays(:,i);
        results(i,1).IRF = IRF;
        results(i,1).startParam = startParam;
        results(i,1).fixedParam = fixedParam;

        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq, results(i,1).residTrace,...
            results(i,1).SSE,results(i,1).exitFlag] = floptimize3_4exp(results(i,1).decays,...
            IRF,cShift,shiftFixed,offFixed,startParam,fixedParam,stFi,period_ns,viewDecay);
    end   
else
    err('Model string entered is currently unrecognized.');
end

%may want to document start point, etc.

end

