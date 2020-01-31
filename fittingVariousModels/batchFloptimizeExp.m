%Wrapper function for repeatedly calling floptimize fitting functions for
%an array of decays that all use the same IRF. The function will call a
%fitting function using models containing 1, 2, 3 or 4 exponential decay
%terms, as specified by the model term of the configuration structure.

%This function takes a set of decays, an IRF, and a configuration structure
% as inputs. It then performs a single fit per decay (i.e. global
%analysis). It returns a structure of the results, as well as some of the
%starting parameters integrated so that the user doesn't always have to go
%back to the saved configuration file to figure out the fit conditions
%used.

function [results] = batchFloptimizeExp(decays,IRF,configS)

%pull the model out of the configuration structure to decide which fit
%function to call and set up the appropriate results structures for it.
model = configS(1,1).model;

%1 exponential fit
if(strcmp(model,'1exp'))
    results = struct('tau',0,'shift',0,'chiSq',0,'offset',0,'SSE',0,'residTrace',0,'exitFlag',0);    

    for i = 1:size(decays,2)
        [results(i,1).tau, results(i,1).shift, results(i,1).offset,...
            results(i,1).chiSq, results(i,1).residTrace,results(i,1).SSE,...
            results(i,1).exitFlag] = floptimize3_1exp(decays(:,i),...
            IRF,configS);
    end

%2 exponential fit   
elseif(strcmp(model,'2exp'))
    results = struct('tau',0,'a',0,'tm',0,'shift',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,'exitFlag',0);    

    for i = 1:size(decays,2)       
        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq,results(i,1).residTrace,...
            results(i,1).SSE, results(i,1).exitFlag] = floptimize3_2exp(decays(:,i),IRF,configS);
    end
    
%3 exponential fit    
elseif(strcmp(model,'3exp'))
    results = struct('tau',0,'a',0,'tm',0,'shift',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParam',0,'fixedParam',0);    

    for i = 1:size(decays,2)
        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq, results(i,1).residTrace, ...
            results(i,1).SSE, results(i,1).exitFlag] = floptimize3_3exp(decays(:,i),...
            IRF,configS);
    end    

%4 exponential fit    
elseif(strcmp(model,'4exp'))
    results = struct('tau',0,'a',0,'tm',0,'shift',0,...
        'chiSq',0,'offset',0,'residTrace',0,'SSE',0,...
        'exitFlag',0,'startParam',0,'fixedParam',0);    

    for i = 1:size(decays,2)
        [results(i,1).tm, results(i,1).a, results(i,1).tau, results(i,1).shift,...
            results(i,1).offset, results(i,1).chiSq, results(i,1).residTrace,...
            results(i,1).SSE,results(i,1).exitFlag] = floptimize3_4exp(decays(:,i),...
            IRF,configS);
    end   
else
    err('Model string entered is currently unrecognized.');
end

end

