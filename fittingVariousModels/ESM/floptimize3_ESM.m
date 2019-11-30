%comment this function here
%TODO - implement removal of non-needed coefficients
%chi squared will be weird given that there are a ton of params - need to
%sort this out

%shiftZero and offsetZero are booleans to decide if the shift and offset
%should be fixed

function [aFinal, chiSq, cF, offset, residTrace, SSE, exitFlag] = floptimize3_ESM(decay,irf,tMax,nSteps,cShift,shiftFixed,offFixed)

%create variables to hold the coefficients for nSteps from 0 ns to tauMax
aSize = 1/nSteps;
tSize = tMax/nSteps;
aList = zeros(nSteps,1);
tList = (tSize:tSize:tMax)';
aList(:) = aSize;

%initialize data parameters correctly for deckard usage
p = 12.5; %time range
n = length(irf);
t = (1:n)'; %time scale in time bin units
dt = p/n;
tp = (0:dt:(p-dt))'; %time scale in ns
st = 20;
fi = 240;
c = cShift; %initialization for shift
residTrace = zeros(size(tp));

if(shiftFixed)
    param = aList;
    nParam = nSteps;
else
    param = aList;
    param(end+1,1) = c; %add shift as the final parameter
    nParam = nSteps +1;
end

%define a figure to show the fit results in
figure(1);

%setup for fmincon
options = optimoptions('fmincon','Display','none','Algorithm','interior-point',...
    'MaxFunctionEvaluations',1e7,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,...
    'MaxIterations',5e5);
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(1,nParam);
ub = zeros(1,nParam);
if (shiftFixed)
    for i=1:nParam
        ub(1,i) = 10;
        lb(1,i) = 0;
    end
else
    for i=1:(nParam-1)
        ub(1,i) = 10;
        lb(1,i) = 0;
    end
    %bound shift between +/- 5 time bin units
    ub(1,nParam) = 5;
    lb(1,nParam) = -5;
end
nonlcon = [];
[pOpt,chiSq,exitFlag] = fmincon(@expScoring,param,A,b,Aeq,beq,lb,ub,nonlcon,options);
 
    function [resid] = expScoring(param)
        %	This nested function calculates the scoring function and has
        %	access to all of the values within floptimize
        
        %isolate the coefficients from the shift and calculate the
        %estimated decay
        if(~shiftFixed)
            %sum up the contribution of each lifetime term to the fit
            %ignore last parameter because it is the shift
            aTotal = sum(param(1:end-1,1));
            aValues = param(1:end-1,1);
            cG = param(end,1); %for re-weighting the coefficients
            x = zeros(n,1);
            for ind=1:size(aValues,1)
                x(:) = x(:) + aValues(ind,1)/aTotal*exp(-tp/tList(ind,1));
            end
            irs = (1-cG+floor(cG))*irf(rem(rem(t-floor(cG)-1, n)+n,n)+1) + (cG-floor(cG))*irf(rem(rem(t-ceil(cG)-1, n)+n,n)+1);
        else
            %sum up the contribution of each lifetime term to the fit
            %all parameters are coefficients
            aTotal = sum(param); %for re-weighting the coefficients
            x = zeros(n,1);
            for ind=1:size(param,1)
                x(:) = x(:) + param(ind,1)/aTotal*exp(-tp/tList(ind,1));
            end
            cG = c;
            irs = (1-cG+floor(cG))*irf(rem(rem(t-floor(cG)-1, n)+n,n)+1) + (cG-floor(cG))*irf(rem(rem(t-ceil(cG)-1, n)+n,n)+1);
        end
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
        resid = sum((z(st:fi) - decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);   
    end

%calculate the final state of the function
if(~shiftFixed)
    %sum up the contribution of each lifetime term to the fit
    %ignore last parameter because it is the shift
    aFinalRaw = pOpt(1:end-1,1);
    aFinalTot = sum(aFinalRaw);
    aFinal = aFinalRaw/aFinalTot;
    cF = pOpt(end,1);
    g = zeros(n,1);
    for j=1:size(aFinal,1)
        g(:) = g(:) + aFinal(j,1)*exp(-tp/tList(j,1));
    end
    irs = (1-cF+floor(cF))*irf(rem(rem(t-floor(cF)-1, n)+n,n)+1) + (cF-floor(cF))*irf(rem(rem(t-ceil(cF)-1, n)+n,n)+1);
else
    %sum up the contribution of each lifetime term to the fit
    %all parameters are coefficients
    aFinalTot = sum(pOpt);
    aFinal = pOpt/aFinalTot;
    g = zeros(n,1);
    for j=1:size(pOpt,1)
        g(:) = g(:) + aFinal(j,1)*exp(-tp/tList(j,1));
    end
    irs = irf;
    cF = c;
end
%re convolve with the correctly shifted IRF
gC = convol(irs, g);

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

%make a plot of the guesses versus the data
hold off
subplot('position',[0.1 0.4 0.8 0.5])
plot(t,log10(decay),t,log10(irs./max(irs).*max(decay)),t,log10(gCW));
v = axis;
v(1) = min(t);
v(2) = max(t);
axis(v);
xlabel('Time in ns');
ylabel('Log Count');
subplot('position',[0.1 0.1 0.8 0.2])
residTrace(st:fi) = (decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi)));
plot(t(st:fi),residTrace(st:fi));
v = axis;
v(1) = min(t);
v(2) = max(t);

axis(v);
xlabel('Time in ns');
ylabel('Residue');
s = sprintf('%3.3f', chiSq);
text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])

%make a plots of the weight at each time bin
figure;
plot(tList,aFinal);
xlabel('Time(ns)');
ylabel('Coefficient Weight');

SSE=0;
for i=st:fi
	%SSE = SSE + (decay(i,1)-gCW(i,1)).^2/abs(gCW(i,1));
    SSE = SSE + (decay(i,1)-gCW(i,1)).^2;
end

end

